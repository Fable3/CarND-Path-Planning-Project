#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable: 4251)
#include <uWS/uWS.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include "helpers.h"
#include "json.hpp"
#include "spline.h"

#include <conio.h>
// for convenience
using nlohmann::json;
using std::string;
using std::vector;

#define NUM_LANES 3

#define EPSILON 1e-5

FILE *fLog = NULL;

void error_condition(const char *filename, int line, const char *text)
{
	std::cerr << "Assert failed in " << filename << " line " << line << ": " << text << std::endl;
}
#undef assert
#define assert(a) if(!(a)) error_condition (__FILE__, __LINE__, #a)

double relaxed_jerk = 5; // m/s^3
double relaxed_acc = 5; // m/s^2
double min_relaxed_acc_while_braking = 4; // m/s^2
double maximum_jerk = 9;
double maximum_acc = 9;

struct Car
{
	// input
	int id;
	double x, y, vx, vy;
	
	// calculated:
	double s, d;
	double acc_x, acc_y;
	double vs, vd;
	int lane;

	double predicted_s(double delta_t)
	{
		return s + vs * delta_t;
	}
	double predicted_d(double delta_t)
	{
		return d + vd * delta_t;
	}
};

struct Map
{
	struct Waypoint
	{
		Point ref; // reference pos
		Point lane_center[NUM_LANES]; // lane centerline position with averaged normalized vectors at joints
		double len;
		double nx, ny; //normal vector (same for ref and all centerlines)
	};
	vector<Waypoint> waypoints;
	double get_lane_center_offset(int lane)
	{
		double lane_width = 4.0;
		return lane_width * (lane + 0.5);
	}
	void Init(const vector<double> &map_waypoints_x,
		const vector<double> &map_waypoints_y)
	{
		waypoints.resize(map_waypoints_x.size());
		// init reference line
		for (int i = 0; i < (int)waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			wp.ref.x = map_waypoints_x[i];
			wp.ref.y = map_waypoints_y[i];
		}
		// init normal vectors
		for (int i = 0; i < (int)waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			auto prev_wp = get_waypoint(i - 1);
			Point delta = wp.ref - prev_wp.ref;
			double delta_len = delta.length();
			wp.nx = delta.y / delta_len;
			wp.ny = -delta.x / delta_len;
		}
		// create centerlines from averaged normals
		for (int i = 0; i < (int)waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			auto next_wp = get_waypoint(i + 1);
			double nx = (wp.nx + next_wp.nx) / 2;
			double ny = (wp.ny + next_wp.ny) / 2;
			// the required lane distance is measured on the perpendicular at the straight part, not at the avgnormal, which is at an angle (alpha)
			double angle_of_normal = atan2(wp.ny, wp.nx);
			double angle_of_avgnormal = atan2(ny, nx);
			double cos_alpha = cos(angle_of_avgnormal - angle_of_normal);
			nx /= cos_alpha;
			ny /= cos_alpha;

			for (int r = 0; r < 3; r++)
			{
				double offset = get_lane_center_offset(r);
				wp.lane_center[r].x = wp.ref.x + nx * offset;
				wp.lane_center[r].y = wp.ref.y + ny * offset;
			}
		}
	}
	int reference_waypoint_id;
	double reference_waypoint_ratio[NUM_LANES]; // 's' calculation zero point for each lane, 0: at prev waypoint, 1: at ref waypoint
	Waypoint &get_waypoint(int idx)
	{
		return waypoints[(idx + waypoints.size()) % waypoints.size()];
	}
	double get_lane_length(int wp, int lane)
	{
		Point diff = get_waypoint(wp).lane_center[lane] - get_waypoint(wp - 1).lane_center[lane];
		return diff.length();
	}
	void init_reference_waypoint(double x, double y)
	{
		int closest_wp = 0;
		Point p(x, y);
		double closest_distsq = (get_waypoint(0).ref - p).lengthsq();
		for (int i = 1; i < (int)waypoints.size(); i++)
		{
			double d = (get_waypoint(i).ref - p).lengthsq();
			if (d < closest_distsq)
			{
				closest_wp = i;
				closest_distsq = d;
			}
		}
		// choose waypoint segment by distance to ref line
		double rnom, snom, rdenom;
		double distsq_ref[2];
		for (int next_or_prev = 0; next_or_prev < 2; next_or_prev++)
		{
			distsq_ref[next_or_prev] = distancesq_pt_seg(p,
				get_waypoint(closest_wp + next_or_prev - 1).ref,
				get_waypoint(closest_wp + next_or_prev).ref,
				rnom, rdenom, snom);
		}
		if (distsq_ref[1] < distsq_ref[0])
		{
			// next segment is closer
			closest_wp++;
		}
		else if (distsq_ref[1] == distsq_ref[0])
		{
			// same distance, so the point is in the triangle around the averaged normal of the 2 segments, decide which side
			Point avg_normal((get_waypoint(closest_wp - 1).nx + get_waypoint(closest_wp).nx) / 2,
				(get_waypoint(closest_wp - 1).ny + get_waypoint(closest_wp).ny) / 2);
			Point delta_p = p - get_waypoint(closest_wp).ref;
			double dotp = avg_normal.x*delta_p.x + avg_normal.y*delta_p.y;
			if (dotp > 0)
			{
				// 'delta_p' vector is after the averaged normal, so it passed the junction of the 2 waypoint segments
				closest_wp++;
			}
		}

		reference_waypoint_id = closest_wp;

		// create reference (0) 's' point on each lane centerline
		for (int lane = 0; lane < NUM_LANES; lane++)
		{
			double distsq = distancesq_pt_seg(p,
				get_waypoint(closest_wp - 1).lane_center[lane],
				get_waypoint(closest_wp).lane_center[lane],
				rnom, rdenom, snom);
			reference_waypoint_ratio[lane] = rnom / rdenom;
		}
	}

	bool lane_matching(double x, double y, double &out_s, double &out_d, int &out_lane, int *p_next_wp_id = NULL, int lane_mask = 0xFFFF) // only works around car
	{
		int search_direction = 0; // any
		Point p(x, y);
		bool stop_search = false;
		int current_wp = reference_waypoint_id;
		double sum_s[NUM_LANES] = { 0, 0, 0 };
		double s_ratio[NUM_LANES]; // ratio shift, reference ratio at start, 1 after going backwards, 0 after forward
		for (int i = 0; i < NUM_LANES; i++) s_ratio[i] = reference_waypoint_ratio[i];

		double best_distsq = 1000 * 1000;
		bool found_any = false;
		for (;;)
		{
			double rnom, snom, rdenom;
			bool search_improved = false;
			for (int lane = 0; lane < NUM_LANES; lane++)
			{
				if (((1 << lane)&lane_mask) == 0) continue;
				double distsq = distancesq_pt_seg(p,
					get_waypoint(current_wp - 1).lane_center[lane],
					get_waypoint(current_wp).lane_center[lane],
					rnom, rdenom, snom);
				if (distsq < best_distsq)
				{
					best_distsq = distsq;
					search_improved = true;
					found_any = true;
					double ratio_from_start = rnom / rdenom;
					double r_mod = ratio_from_start - s_ratio[lane];
					double seg_len = get_lane_length(current_wp, lane); // same as sqrt(rdenom)
					out_s = sum_s[lane] + seg_len * r_mod;
					double d = sqrt(distsq);
					if (snom < 0) d = -d; // signed distance
					out_d = d + get_lane_center_offset(lane);
					out_lane = lane;
					if (p_next_wp_id != NULL) *p_next_wp_id = current_wp;
				}
				if (rnom == 0) // go towards prev
				{
					if (search_direction == 1) stop_search = true; // change direction
					search_direction = -1;
				}
				else if (rnom == rdenom)
				{
					if (search_direction == -1) stop_search = true; // change direction
					search_direction = 1;
				}
				else
				{
					// projection is on the segment, stop search
					stop_search = true;
				}
			}
			if (!search_improved || stop_search) break;
			assert(search_direction != 0);
			if (search_direction > 0)
			{
				for (int lane = 0; lane < NUM_LANES; lane++)
				{
					sum_s[lane] += (1 - s_ratio[lane])*get_lane_length(current_wp, lane);
					s_ratio[lane] = 0; // at start of next waypoint, ratio shift is 0
				}
				current_wp++;
			}
			else
			{
				for (int lane = 0; lane < NUM_LANES; lane++)
				{
					sum_s[lane] -= s_ratio[lane] * get_lane_length(current_wp, lane);
					s_ratio[lane] = 1; // because we're at the end of the previous waypoint segment
				}
				current_wp--;
			}
		}
		return found_any;
	}

	Point get_lane_pos(double s, int lane, int &next_waypoint_id, double &next_waypoint_distance)
	{
		double ratio = reference_waypoint_ratio[lane];
		int wp_id = reference_waypoint_id;
		Point next_pt, prev_pt;
		double dest_ratio = 0; // result is between prev_pt and next_pt, at next_pt if 1.0

		for (;;)
		{
			next_pt = get_waypoint(wp_id).lane_center[lane];
			prev_pt = get_waypoint(wp_id - 1).lane_center[lane];
			double wp_len = (next_pt - prev_pt).length();
			if (s > 0)
			{
				double remaining_length = wp_len * (1 - ratio);
				if (s <= remaining_length)
				{
					// if wp_len is 100m, ratio is 0.4, s is 50:
					//    remaining_length = 60
					//    dest_ratio = 1-(60-50)/100 or 0.4+50/60*0.6
					dest_ratio = 1 - (remaining_length - s) / wp_len;
					next_waypoint_distance = remaining_length - s;
					break;
				}
				s -= remaining_length;
				ratio = 0;
				wp_id++;
			}
			else
			{
				double remaining_length = wp_len * ratio;
				if (-s <= remaining_length)
				{
					// if wp_len is 100m, ratio is 0.4, s is -30:
					//    remaining_length = 40
					//    dest_ratio = (40-30)/100 or (40-30)/40*0.4
					//    next_waypoint_distance = 100*(1-0.4)+30
					dest_ratio = (remaining_length + s) / wp_len;
					next_waypoint_distance = wp_len*(1-ratio) - s;
					break;
				}
				s += remaining_length;
				ratio = 1;
				wp_id--;
			}
		}
		Point ret;
		ret.x = next_pt.x*dest_ratio + prev_pt.x*(1 - dest_ratio);
		ret.y = next_pt.y*dest_ratio + prev_pt.y*(1 - dest_ratio);
		next_waypoint_id = wp_id;
		return ret;
	}
};


struct SpeedController
{
	double start_speed;
	double target_speed;
	double target_time;
	double get_speed(double current_t)
	{
		double speed;
		if (current_t > target_time) speed = target_speed;
		else speed = start_speed + (target_speed - start_speed)*current_t / target_time;
		return speed;
	}
	void maximize_acc(double acc)
	{
		double delta_v = target_speed - start_speed;
		double min_time = fabs(delta_v) / acc;
		if (target_time < min_time)
		{
			target_time = min_time;
		}
	}
};


struct TrajectoryBuilder
{
	double total_dist = 0;
	vector<Point > waypoints;
	void add_waypoint(Point p)
	{
		if (waypoints.size() > 0)
		{
			auto prev_p = waypoints.back();
			total_dist += distance(prev_p.x, prev_p.y, p.x, p.y);
		}
		waypoints.push_back(p);
	}
	
	vector<Point> build(vector<Point> &prev_trajectory,
		double ego_x, double ego_y, double ego_yaw, // fallback for init if not enough points in prev_trajectory
		int ego_lane, int target_lane,
		Map &map,
		SpeedController &speed_controller)
	{
		/**
		 *  define a path made up of (x,y) points that the car will visit
		 *   sequentially every .02 seconds
		 */
		double pos_x;
		double pos_y;
		double angle;
		vector<Point> result_points = prev_trajectory;

		int path_size = result_points.size();

		if (path_size == 0) {
			pos_x = ego_x;
			pos_y = ego_y;
			angle = deg2rad(ego_yaw);
		}
		else 
		{
			pos_x = result_points[path_size - 1].x;
			pos_y = result_points[path_size - 1].y;
			if (path_size == 1)
			{
				angle = deg2rad(ego_yaw);
			}
			else
			{
				double pos_x2 = result_points[path_size - 2].x;;
				double pos_y2 = result_points[path_size - 2].y;
				angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
			}
		}
		add_waypoint(Point(pos_x, pos_y)); // startup

		/*
		int nNextWaypoint = map.reference_waypoint_id;
		//ego_lane = 2;
		// skip next waypoint if it's between car pos and end of fixed trajectory head
		{
			double s, d;
			int lane;
			map.lane_matching(pos_x, pos_y, s, d, lane, &nNextWaypoint, 1 << ego_lane);
		}
		//if (map.reference_waypoint_ratio[ego_lane] > 0.9) nNextWaypoint++; // prevent flukes
		
		for (int i = 0; i < 50 - path_size; ++i) {
			Point p = map.get_waypoint(nNextWaypoint).lane_center[ego_lane];
			add_waypoint(p.x, p.y);
			if (total_dist > 50 && i>2) break;
			nNextWaypoint = (nNextWaypoint + 1) % map.waypoints.size();
		}*/
		double min_control_point_dist = speed_controller.start_speed * 1;
		min_control_point_dist = std::max(min_control_point_dist, 5.0); // at least 5 m
		double contol_point_start_s = 0;
		if (ego_lane != target_lane)
		{
			double lane_switch_time = 2.0;
			double min_lane_switch_dist = 5.0;
			double speed = speed_controller.start_speed;
			double dist = speed * lane_switch_time;
			if (dist < min_lane_switch_dist) dist = min_lane_switch_dist;
			double waypoint_dist;
			int next_wp;
			add_waypoint(map.get_lane_pos(dist, target_lane, next_wp, waypoint_dist));
			contol_point_start_s = dist;
		}
		/*else
		{
			// check if next waypoint is close, keep it as a waypoint
			double waypoint_dist;
			int next_wp;
			map.get_lane_pos(contol_point_start_s, target_lane, next_wp, waypoint_dist);
			if (waypoint_dist < min_control_point_dist && waypoint_dist>5)
			{
				contol_point_start_s += waypoint_dist - min_control_point_dist;
			}
		}*/
		for (int i = 0; i < 5; i++)
		{
			double waypoint_dist;
			int next_wp;
			contol_point_start_s += min_control_point_dist;
			Point next_pt = map.get_lane_pos(contol_point_start_s, target_lane, next_wp, waypoint_dist);
			/*if (waypoint_dist < min_control_point_dist)
			{
				// if the next waypoint is closer than double min_control_point_dist, we skip the next control point,
				// and add the waypoint instead to avoid cutting corners, or adding too dense control points which would make it jerk
				double tmp_dist;
				int tmp_wp;
				//next_pt = map.get_lane_pos(contol_point_start_s- min_control_point_dist/2+waypoint_dist/2, target_lane, tmp_wp, tmp_dist);
				//add_waypoint(next_pt);
				add_waypoint(map.get_waypoint(next_wp).lane_center[target_lane]);
				contol_point_start_s += waypoint_dist;
			}
			else*/
			{
				add_waypoint(next_pt);
			}
			if (total_dist > 50 && waypoints.size() > 2) break;
		}

		if (fLog)
		{
			fprintf(fLog, "first control dist %.2f\n", (Point(pos_x, pos_y) - waypoints[1]).length());
			fprintf(fLog, "prev_trajectory=[");
			for (int i = 0; i < prev_trajectory.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", prev_trajectory[i].x, prev_trajectory[i].y);
			fprintf(fLog, "]\n");
			fprintf(fLog, "control_points=[");
			for (int i = 0; i < waypoints.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", waypoints[i].x, waypoints[i].y);
			fprintf(fLog, "]\n");
		}


		// transform waypoints for spline
		double cos_angle = cos(-angle);
		double sin_angle = sin(-angle);
		Point transform_center(pos_x, pos_y);
		for (int i = 0; i < (int)waypoints.size(); i++)
		{
			Point &p = waypoints[i];
			p = p - transform_center;
			double tx = p.x*cos_angle - p.y*sin_angle;
			double ty = p.x*sin_angle + p.y*cos_angle;
			p.x = tx;
			p.y = ty;
		}
		vector<double> wp_x, wp_y;
		for (int i = 0; i < int(prev_trajectory.size()) - 1; i++)
		{
			Point &p = prev_trajectory[i];
			Point t = p - transform_center;
			double tx = t.x*cos_angle - t.y*sin_angle;
			double ty = t.x*sin_angle + t.y*cos_angle;

			wp_x.push_back(tx);
			wp_y.push_back(ty);
		}

		pos_x = 0;
		pos_y = 0;
		cos_angle = cos(angle);
		sin_angle = sin(angle);

		tk::spline spline;
		
		for (int i = 0; i < (int)waypoints.size(); i++)
		{
			wp_x.push_back(waypoints[i].x);
			wp_y.push_back(waypoints[i].y);
		}

		for (int i = 1; i < (int)wp_x.size(); i++)
		{
			if (wp_x[i] < wp_x[i - 1])
			{
				printf("\nspline input error\n");
				if (fLog) fprintf(fLog, "spline input error at idx %d x values %.4f<%.4f\n", i, wp_x[i], wp_x[i - 1]);
				wp_x.resize(i);
				wp_y.resize(i);
				break;
			}
		}
		if (wp_x.size() < 3)
		{
			return prev_trajectory;
		}
		
		spline.set_points(wp_x, wp_y);

		double dist_total = 0;
		double current_t = 0.02;
		//double speed = max_speed;
		//for (int i = 1; i < (int)waypoints.size(); i++)
		double current_sp_arg = 0;
		while (current_sp_arg<50) // 50 meters max, will stop at 50 points anyway
		{
			/*if (!in_lane_jmt.empty() && current_t <= in_lane_jmt_max_time)
			{
				speed = get_jmt_speed(in_lane_jmt, current_t);
			}*/
			// do linear speed, the only important point is the first one, we'll drop the rest in the next cycle,
			// so can't do smooth stuff
			double speed = speed_controller.get_speed(current_t);
			double dist_step = speed / 50;
			double sp_result = spline(current_sp_arg + dist_step);
			double x = current_sp_arg + dist_step;
			double y = sp_result;
			double dist = distance(pos_x, pos_y, x, y);
			if (dist+EPSILON < dist_step)
			{
				printf("\nspline warning, dist: %.4f, step: %.4f, pos %.2f,%.2f to %.2f,%.2f\n", dist, dist_step, pos_x, pos_y, x, y);
				if (fLog) fprintf(fLog, "spline warning, dist: %.4f, step: %.4f, pos %.2f,%.2f to %.2f,%.2f\n", dist, dist_step, pos_x, pos_y, x, y);
			}
			if (current_t == 0.02) printf(" next sp %.2f", speed);
			current_t += 0.02;
			double sp_step = (x - pos_x)*dist_step / dist;
			current_sp_arg += sp_step;
			pos_x += sp_step;
			pos_y += (y - pos_y)*dist_step / dist;
			// transform back:
			Point tp;
			tp.x = pos_x * cos_angle - pos_y * sin_angle;
			tp.y = pos_x * sin_angle + pos_y * cos_angle;
			tp = tp + transform_center;				
			result_points.push_back(tp);
			if (result_points.size() >= 50) break;
		}
		if (fLog)
		{
			fprintf(fLog, "result=[");
			for (int i = 0; i < result_points.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i==0?"":",", result_points[i].x, result_points[i].y);
			fprintf(fLog, "]\n");
		}
		return result_points;
	}
};



//int jump_to_waypoint = 30; // traffic block
int jump_to_waypoint = 0;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    std::istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  Map map;
  map.Init(map_waypoints_x, map_waypoints_y);
  double prev_car_speed = 0;
  double prev_car_acc = 0;
  bool prev_speed_valid = false;
  bool prev_acc_valid = false;
  std::map<int, Car> sensor_fusion_cars;
  int target_lane = 1;

  bool need_log = true;
  if (need_log)
  {
	  fLog = fopen("trajectory.log", "wt");
	  fprintf(fLog, "wpmap=[");
	  for (int i = 0; i < map_waypoints_x.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", map_waypoints_x[i], map_waypoints_y[i]);
	  fprintf(fLog, "]\n");
  }

  h.onMessage([&map, &prev_car_speed, &prev_speed_valid, &prev_car_acc, &prev_acc_valid, &sensor_fusion_cars, &target_lane]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
               uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double ego_x = j[1]["x"];
          double ego_y = j[1]["y"];
          double ego_s_orig = j[1]["s"];
          double ego_d_orig = j[1]["d"];
          double ego_yaw = j[1]["yaw"];
          double ego_speed = j[1]["speed"] / 2.237;
		  double ego_acc=0, ego_jerk=0;

		  if (prev_speed_valid)
		  {
			  ego_acc = (ego_speed - prev_car_speed) * 50;
			  if (prev_acc_valid)
			  {
				  ego_jerk = (ego_acc - prev_car_acc) * 50;
			  }
			  prev_car_acc = ego_acc;
			  prev_acc_valid = true;
		  }
		  prev_car_speed = ego_speed;
		  prev_speed_valid = true;

		  // doesn't work from car speed, end of prev path would be better
		  ego_acc = 0;
		  ego_jerk = 0;
		  
		  // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

		  vector<Point> prev_trajectory;
		  double delta_t0 = 0; // time in the future where trajectory planning starts, at the end of kept previous path

		  int prev_trajectory_length = 4; // need 1 for pos, 2 for speed, 3 for acc, 4 for jerk
		  if ((int)previous_path_x.size() >= prev_trajectory_length)
		  {
			  for (int i = 0; i < prev_trajectory_length; i++)
			  {
				  Point p(previous_path_x[i], previous_path_y[i]);
				  prev_trajectory.push_back(p);
			  }
			  double v1 = (prev_trajectory[1] - prev_trajectory[0]).length();
			  double v2 = (prev_trajectory[2] - prev_trajectory[1]).length();
			  double v3 = (prev_trajectory[3] - prev_trajectory[2]).length();
			  ego_acc = (v3 - v2) * 50;
			  double prev_acc = (v2 - v1) * 50;
			  ego_jerk = (ego_acc - prev_acc) * 50;
			  ego_speed = v3 * 50;
			  ego_x = prev_trajectory[3].x;
			  ego_y = prev_trajectory[3].y;
			  delta_t0 = prev_trajectory_length / 50.0;
		  }
		  else
		  {
			  if (jump_to_waypoint)
			  {
				  Point p = map.get_waypoint(jump_to_waypoint).lane_center[1];
				  ego_x = p.x;
				  ego_y = p.y;
				  //jump_to_waypoint = 0;
			  }
		  }

		  //printf("v % 4.2f acc % 4.2f jerk % 4.2f", ego_speed, ego_acc, ego_jerk);
		  // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

		  map.init_reference_waypoint(ego_x, ego_y);
		  int ego_lane;
		  double ego_s, ego_d;
		  if (!map.lane_matching(ego_x, ego_y, ego_s, ego_d, ego_lane))
		  {
			  printf("Warning! can't lane match ego\n");
			  ego_s = ego_d = 0;
			  ego_lane = 0;
		  }
		  printf("v% 4.2f a% 4.2f d% 4.2f", ego_speed, ego_acc, ego_d);
		  if (fLog) fprintf(fLog, "v%.2f a%.2f d%.2f\n", ego_speed, ego_acc, ego_d);

		  if (ego_acc > maximum_acc) ego_acc = maximum_acc; // don't try to limit jerk if acc is too high anyway
		  if (ego_acc < -maximum_acc) ego_acc = -maximum_acc;
		  if (ego_jerk > maximum_jerk) ego_jerk = maximum_jerk;
		  if (ego_jerk < -maximum_jerk) ego_jerk = -maximum_jerk;


		  for (int i = 0; i < (int)sensor_fusion.size(); i++)
		  {
			  auto car_data = sensor_fusion[i];
			  int id = car_data[0];
			  auto &car = sensor_fusion_cars[id];
			  car.id = id;
			  car.x = car_data[1];
			  car.y= car_data[2];
			  car.vx = car_data[3];
			  car.vy = car_data[4];
			  int next_wp_id = 0;
			  if (!map.lane_matching(car.x, car.y, car.s, car.d, car.lane, &next_wp_id))
			  {
				  printf("Warning! can't lane match car %d at %.2f %.2f\n", id, car.x, car.y);
				  sensor_fusion_cars.erase(sensor_fusion_cars.find(id));
			  }
			  else
			  {
				  //car.s = car_data[5];
				  //car.d = car_data[6];
				  double rnom, snom, rdenom;
				  Point wp_vector = map.get_waypoint(next_wp_id).ref - map.get_waypoint(next_wp_id-1).ref;
				  Point speed_vector = Point(car.vx, car.vy);
				  double speed_vector_len = speed_vector.length();
				  if (speed_vector_len < EPSILON)
				  {
					  car.vs = speed_vector_len;
					  car.vd = 0;
				  }
				  else
				  {
					  double wp_vector_len = wp_vector.length();
					  wp_vector.x *= speed_vector_len / wp_vector_len;
					  wp_vector.y *= speed_vector_len / wp_vector_len;
					  double speed_sign = 1.0;

					  if (wp_vector.dotp(speed_vector) < 0)
					  {
						  speed_vector.x *= -1;
						  speed_vector.y *= -1;
						  speed_sign = -1;
					  }
					  distancesq_pt_seg(Point(car.vx, car.vy), Point(0, 0), wp_vector, rnom, rdenom, snom);
					  car.vs = (rnom / rdenom)*speed_vector_len * speed_sign;
					  car.vd = (snom / rdenom)*speed_vector_len * speed_sign;
				  }
			  }
		  }

		  int next_car_id = -1;
		  double next_s = 0;
		  for (auto &p : sensor_fusion_cars)
		  {
			  auto &other = p.second;
			  double s0 = other.predicted_s(delta_t0);
			  double d0 = other.predicted_d(delta_t0);
			  if (s0 > ego_s && fabs(d0 - ego_d) < 3)
			  {
				  if (next_car_id == -1 || next_s > s0)
				  {
					  next_car_id = other.id;
					  next_s = s0;
				  }
			  }
		  }
		  double car_length = 4.5;
		  double safety_distance = 2;
		  double keep_distance = 5;
		  double keep_distance_leeway = 0.5; // when following, match car in front speed if in this range
		  double max_speed = 22.2; // 50 mph
		  if (next_car_id != -1)
		  {
			  auto &car = sensor_fusion_cars[next_car_id];
			  printf(" next (%.2f, %.2f) v (%.2f, %.2f) s_dist: %.2f lane %d",
				  car.s, car.d, car.vs, car.vd, next_s-ego_s, car.lane);
		  }
		  //vector<double> in_lane_jmt;
		  //double in_lane_jmt_max_time = 0;
		  SpeedController speed_controller;
		  speed_controller.start_speed = ego_speed;
		  speed_controller.target_speed = max_speed;
		  double delta_to_max_speed = fabs(ego_speed - max_speed);
		  speed_controller.target_time = delta_to_max_speed / relaxed_acc;
		  if (next_car_id != -1)
		  {
			  if (next_s - ego_s < 6)
			  {
				  int alma = 0;
			  }
			  vector<double> start = { ego_s, ego_speed, ego_acc };

			  auto &follow_car = sensor_fusion_cars[next_car_id];
			  double follow_car_distance = next_s - ego_s - car_length;
			  if (follow_car_distance < 0)
			  {
				  printf("detected collision");
				  follow_car_distance = 0; // collision already
			  }
			  double follow_car_speed = sqrt(follow_car.vx*follow_car.vx + follow_car.vy*follow_car.vy);
			  bool can_accelerate = true;
			  // check emergency breaking, safe distance
			  if (ego_speed > follow_car_speed)
			  {
				  double acc = relaxed_acc;
				  if (ego_acc < 0) acc = min_relaxed_acc_while_braking;
				  double delta_v = ego_speed - follow_car_speed;
				  double delta_t = delta_v / acc;
				  double delta_d = ego_speed * delta_t - delta_v / 2 * delta_t;
				  double max_dist = follow_car_distance - safety_distance;
				  /*if (delta_d > max_dist)
				  {
					  // need to start breaking. The most efficient way to maintain safety distance is to pretend that the car ahead breaks at maximum deceleration and stops
					  double time_for_follow_car_to_stop = follow_car_speed / maximum_acc;
					  double extra_distance_while_stopping = time_for_follow_car_to_stop * follow_car_speed / 2;
					  double follow_car_predicted_pos = next_s + extra_distance_while_stopping;
					  // we reach the car pos in delta_t, plus the time to stop
					  vector<double> end = { follow_car_predicted_pos - car_length - safety_distance, 0, 0 };
					  in_lane_jmt_max_time = time_for_follow_car_to_stop + delta_t;
					  in_lane_jmt = JMT(start, end, in_lane_jmt_max_time);
					  can_accelerate = false;

				  }*/
				  if (delta_d > max_dist)
				  {
					  speed_controller.target_speed = follow_car_speed;
					  // we have max_dist to decelerate, so from max_dist = (ego_speed - delta_v / 2)*target_speed_time:
					  speed_controller.target_time = max_dist / (ego_speed - delta_v / 2);
					  // if max accel is reached, use up the safety distance:
					  if (speed_controller.target_time < EPSILON || delta_v / speed_controller.target_time>maximum_acc)
					  {
						  printf(" MAXBRAKE");
						  speed_controller.target_time = delta_v / maximum_acc;
					  }
					  else
					  {
						  printf(" BRAKE");
					  }
					  can_accelerate = false;
				  }
			  }
			  if (can_accelerate) // maintain safe distance
			  {
				  double t_to_reach_optimal_distance = 1; // sec
				  double follow_car_predicted_pos = follow_car_speed * t_to_reach_optimal_distance + next_s;
				  double ego_car_predicted_pos = ego_speed * t_to_reach_optimal_distance + ego_s;

				  if (ego_s + car_length + keep_distance > next_s)
				  {
					  double excess_distance = ego_s + car_length + keep_distance - next_s; // negative when approaching

					  speed_controller.target_speed = follow_car_speed - excess_distance / t_to_reach_optimal_distance;
					  speed_controller.target_time = t_to_reach_optimal_distance;
					  speed_controller.maximize_acc(relaxed_acc);
					  can_accelerate = false;
					  printf(" ADJUST");

					  /*
					  vector<double> start = { ego_s, ego_speed, ego_acc };
					  vector<double> end = { follow_car_predicted_pos - car_length, follow_car_speed, 0 };
					  in_lane_jmt = JMT(start, end, t_to_reach_optimal_distance);*/
				  }
				  else if (ego_s + car_length + keep_distance + keep_distance_leeway > next_s)
				  {
					  // we maintain current speed, don't try to exacly maintain safety distance
					  speed_controller.target_speed = follow_car_speed;
					  speed_controller.target_time = 1.0;
					  speed_controller.maximize_acc(relaxed_acc);
					  printf(" KEEP");
				  }
			  }
			  /*double speed_diff = fabs(follow_car_speed - ego_speed);
			  if (ego_speed > follow_car_speed && can_accelerate)
			  {
				  double collision_time = follow_car_distance / (ego_speed - follow_car_speed);
				  if (collision_time > 10)
				  {
					  // follow car is far, can accelerate
				  }
				  else
				  {
					  // need to break, follow mode
					  double follow_car_predicted_pos = follow_car_speed * collision_time + next_s;
					  vector<double> end = { follow_car_predicted_pos - car_length - safety_distance, follow_car_speed, 0 };
					  in_lane_jmt_max_time = collision_time;
					  in_lane_jmt = JMT(start, end, in_lane_jmt_max_time);
					  
					  can_accelerate = false;
				  }
			  }
			  double follow_sync_time = 3; // if speed is no problem, match distance if it's in 3 seconds range
			  if (can_accelerate && follow_car_distance<ego_speed*follow_sync_time) 
			  {
				  double follow_car_predicted_pos = follow_car_speed * follow_sync_time + next_s;
				  vector<double> end = { follow_car_predicted_pos - car_length - safety_distance, follow_car_speed, 0 };
				  in_lane_jmt_max_time = follow_sync_time;
				  in_lane_jmt = JMT(start, end, in_lane_jmt_max_time);
				  can_accelerate = false;
			  }*/
			  /*
			  double t_to_reach_optimal_distance = 1; // sec
			  double follow_car_predicted_pos = follow_car_speed * t_to_reach_optimal_distance + next_s;
			  double ego_car_predicted_pos = ego_speed * t_to_reach_optimal_distance + ego_s;

			  if (ego_s + car_length + safety_distance > next_s)
			  {
				  

				  double excess_distance = ego_s + car_length - next_s; // negative when approaching
				  
				  
				  max_speed = std::min(max_speed, follow_car_speed) - excess_distance / t_to_reach_optimal_distance;

				  // should start from pos_x
				  vector<double> start = { ego_s, ego_speed, ego_acc };
				  vector<double> end = { follow_car_predicted_pos - car_length, follow_car_speed, 0 };
				  in_lane_jmt = JMT(start, end, t_to_reach_optimal_distance);
			  }*/
		  }
		  printf(" target %.2f @ %.2f", speed_controller.target_speed, speed_controller.target_time);
		  //if (in_lane_jmt.empty()) 
		  //else printf(" JMT start speed %.2f end speed %.2f (%.2fs)", get_jmt_speed(in_lane_jmt, 0), get_jmt_speed(in_lane_jmt, in_lane_jmt_max_time), in_lane_jmt_max_time);

          json msgJson;

		  TrajectoryBuilder trajectory;
		  
		  if (_kbhit())
		  {
			  char ch = _getch();
			  if (ch == 'A' || ch == 'a') target_lane = std::max(0, target_lane - 1);
			  if (ch == 'D' || ch == 'd') target_lane = std::min(2, target_lane + 1);
			  
		  }
		  vector<Point> refined_trajectory = trajectory.build(prev_trajectory, ego_x, ego_y, ego_yaw, ego_lane, target_lane, map, speed_controller);

		  vector<double> next_x_vals(refined_trajectory.size());
		  vector<double> next_y_vals(refined_trajectory.size());

		  for (int i = 0; i < (int)refined_trajectory.size(); i++)
		  {
			  next_x_vals[i] = refined_trajectory[i].x;
			  next_y_vals[i] = refined_trajectory[i].y;
		  }

		  /*if (fLog != NULL)
		  {
			  double prev_heading = 0;
			  for (int i = 0; i < (int)next_x_vals.size(); i++)
			  {
				  double heading = 0;
				  double dist = 0;
				  if (i != 0)
				  {
					  heading = atan2(next_y_vals[i] - next_y_vals[i - 1], next_x_vals[i] - next_x_vals[i - 1]);
					  dist = distance(next_x_vals[i], next_y_vals[i], next_x_vals[i - 1], next_y_vals[i - 1]);
				  }
				  fprintf(fLog, "%.2f, %.2f, %.2f, %2f\n", next_x_vals[i], next_y_vals[i], dist, heading * 180 / pi());
				  if (i >1 && fabs(heading - prev_heading) > 0.2)
				  {
					  int alma = 0;
				  }
				  prev_heading = heading;
			  }
			  fclose(fLog);
		  }*/
		  if (fLog) fflush(fLog);
		  printf("\n");
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }  // end "telemetry" if
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket if
  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen("127.0.0.1", port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  
  h.run();
}