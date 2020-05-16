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

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

#define NUM_LANES 3

void error_condition(const char *filename, int line, const char *text)
{
	std::cerr << "Assert failed in " << filename << " line " << line << ": " << text << std::endl;
}
#undef assert
#define assert(a) if(!(a)) error_condition (__FILE__, __LINE__, #a)

struct Car
{
	int id;
	double x, y, vx, vy, s, d;
};
std::map<int, Car> SensorFusionCars;

struct TrajectoryBuilder
{
	double total_dist = 0;
	vector<Point > waypoints;
	void add_point(double x, double y)
	{
		if (waypoints.size() > 0)
		{
			auto prev_p = waypoints.back();
			total_dist += distance(prev_p.x, prev_p.y, x, y);
		}
		waypoints.push_back(Point(x, y));
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
		for (int i = 0; (int)i < waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			wp.ref.x = map_waypoints_x[i];
			wp.ref.y = map_waypoints_y[i];
		}
		// init normal vectors
		for (int i = 0; (int)i < waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			auto prev_wp = get_waypoint(i - 1);
			Point delta = wp.ref - prev_wp.ref;
			double delta_len = delta.length();
			wp.nx = delta.y/delta_len;
			wp.ny = -delta.x/ delta_len;
		}
		// create centerlines from averaged normals
		for (int i = 0; (int)i < waypoints.size(); i++)
		{
			auto &wp = waypoints[i];
			auto next_wp = get_waypoint(i + 1);
			double nx = (wp.nx + next_wp.nx) / 2;
			double ny = (wp.ny + next_wp.ny) / 2;
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
				get_waypoint(closest_wp- 1).lane_center[lane],
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
					sum_s[lane] -= s_ratio[lane]*get_lane_length(current_wp, lane);
					s_ratio[lane] = 1; // because we're at the end of the previous waypoint segment
				}
				current_wp--;
			}
		}
		return found_any;
	}
};

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

  h.onMessage([&map]
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
          double ego_s = j[1]["s"];
          double ego_d = j[1]["d"];
          double ego_yaw = j[1]["yaw"];
          double ego_speed = j[1]["speed"];

		  // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

		  map.init_reference_waypoint(ego_x, ego_y);

		  for (int i = 0; i < (int)sensor_fusion.size(); i++)
		  {
			  auto car_data = sensor_fusion[i];
			  int id = car_data[0];
			  auto &car = SensorFusionCars[id];
			  car.id = id;
			  car.x = car_data[1];
			  car.y= car_data[2];
			  car.vx = car_data[3];
			  car.vy = car_data[4];
			  car.s = car_data[5];
			  car.d = car_data[6];
		  }

		  int next_car_id = -1;
		  double next_s = 0;
		  for (auto &p : SensorFusionCars)
		  {
			  auto &car = p.second;
			  if (car.s > ego_s && fabs(car.d - ego_d) < 2)
			  {
				  if (next_car_id == -1 || next_s > car.s)
				  {
					  next_car_id = car.id;
					  next_s = car.s;
				  }
			  }
		  }
		  double car_length = 7;
		  double safety_distance = 5;
		  double max_speed = 44.4;//22.2;
		  if (next_car_id != -1)
		  {
			  auto &car = SensorFusionCars[next_car_id];
			  double s, d;
			  int lane;
			  map.lane_matching(car.x, car.y, s, d, lane);
			  double cartesian_dist = distance(car.x, car.y, ego_x, ego_y);
			  printf("next car telemetry (%.2f, %.2f) calculated (%.2f, %.2f) diff %.2f xy_dist: %.2f dist: %.2f in lane %d\n",
				  car.s, car.d, s + ego_s, d, fabs(s + ego_s - car.s), s, cartesian_dist, lane);
			  ego_s = 0;
			  next_s = s;
		  }
		  if (next_car_id != -1 && ego_s+car_length+safety_distance>next_s)
		  {
			  auto &car = SensorFusionCars[next_car_id];
			  double t_to_reach_optimal_distance = 1; // sec
			  double excess_distance = ego_s + car_length - next_s; // negative when approaching
			  max_speed = std::min(max_speed, sqrt(car.vx*car.vx + car.vy*car.vy)) - excess_distance / t_to_reach_optimal_distance;
		  }


          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          /**
           *  define a path made up of (x,y) points that the car will visit
           *   sequentially every .02 seconds
           */

		  
		  double pos_x;
		  double pos_y;
		  double angle;
		  int path_size = std::min(3, (int)previous_path_x.size());

		  for (int i = 0; i < path_size; ++i) {
			  next_x_vals.push_back(previous_path_x[i]);
			  next_y_vals.push_back(previous_path_y[i]);
		  }

		  TrajectoryBuilder trajectory;
		  if (path_size == 0) {
			  pos_x = ego_x;
			  pos_y = ego_y;
			  angle = deg2rad(ego_yaw);
		  }
		  else {
			  pos_x = previous_path_x[path_size - 1];
			  pos_y = previous_path_y[path_size - 1];

			  double pos_x2 = previous_path_x[path_size - 2];
			  double pos_y2 = previous_path_y[path_size - 2];
			  angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
		  }
		  trajectory.add_point(pos_x, pos_y);
		  
		  int nNextWaypoint = map.reference_waypoint_id;
		  int ego_lane = 2;
		  // skip next waypoint if it's between car pos and end of fixed trajectory head
		  {
			  double s, d;
			  int lane;
			  map.lane_matching(pos_x, pos_y, s, d, lane, &nNextWaypoint, 1<<ego_lane);
		  }
		  //if (map.reference_waypoint_ratio[ego_lane] > 0.9) nNextWaypoint++; // prevent flukes

		  for (int i = 0; i < 50 - path_size; ++i) {
			  Point p = map.get_waypoint(nNextWaypoint).lane_center[ego_lane];
			  trajectory.add_point(p.x, p.y);
			  if (trajectory.total_dist > 50) break;
			  nNextWaypoint = (nNextWaypoint + 1) % map.waypoints.size();
		  }

		  double dist_total = 0;
		  double dist_step = max_speed/50;
		  for (int i = 1; i < trajectory.waypoints.size(); i++)
		  {
			  double x = trajectory.waypoints[i].x;
			  double y = trajectory.waypoints[i].y;
			  for (;;)
			  {
				  double dist = distance(pos_x, pos_y, x, y);
				  if (dist < dist_step) break;
				  pos_x += (x-pos_x)*dist_step / dist;
				  pos_y += (y-pos_y)*dist_step / dist;
				  next_x_vals.push_back(pos_x);
				  next_y_vals.push_back(pos_y);
				  if (next_x_vals.size() >= 50) break;
			  }
			  if (next_x_vals.size() >= 50) break;
		  }

		  FILE *fLog = fopen("trajectory.log", "wt");
		  if (fLog != NULL)
		  {
			  double prev_heading = 0;
			  for (int i = 0; i < next_x_vals.size(); i++)
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
		  }


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