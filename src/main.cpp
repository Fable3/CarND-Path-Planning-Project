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

#define EPSILON 1e-5

FILE *fLog = NULL;
bool need_log = true;
bool console_log = true;
bool sparse_console_log = false;
bool test_fast_lane_change = false;

void error_condition(const char *filename, int line, const char *text)
{
	std::cerr << "Assert failed in " << filename << " line " << line << ": " << text << std::endl;
}
#undef assert
#define assert(a) if(!(a)) error_condition (__FILE__, __LINE__, #a)

double relaxed_acc = 5; // m/s^2
double min_relaxed_acc_while_braking = 4; // m/s^2
double maximum_jerk = 9;
double maximum_acc = 7;

double max_speed = 22.2; // 50 mph
double car_length = 4.5;
double safety_distance = 2;
double keep_distance = 10;
double keep_distance_leeway = 0.5; // when following, match car in front speed if in this range

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

class Map
{
public:
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

	void project_speed(Point speed_vector, int next_wp_id, double *p_vs, double *p_vd)
	{
		double rnom, snom, rdenom;
		Point wp_vector = get_waypoint(next_wp_id).ref - get_waypoint(next_wp_id - 1).ref;
		
		double speed_vector_len = speed_vector.length();
		if (speed_vector_len < EPSILON)
		{
			*p_vs = speed_vector_len;
			*p_vd = 0;
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
			distancesq_pt_seg(speed_vector, Point(0, 0), wp_vector, rnom, rdenom, snom);
			*p_vs = (rnom / rdenom)*speed_vector_len * speed_sign;
			*p_vd = (snom / rdenom)*speed_vector_len * speed_sign;
		}
	}
};

class LaneChangePlanner
{
public:
	int calculate_target_lane(std::map<int, Car> &sensor_fusion_cars, int ego_lane, int target_lane, double ego_s, double ego_vs, double delta_t0)
	{
		double lane_speed[NUM_LANES];
		double next_s_in_lane[NUM_LANES];
		bool lane_open[NUM_LANES];
		for (int i = 0; i < NUM_LANES; i++)
		{
			lane_speed[i] = max_speed;
			next_s_in_lane[i] = 1000;
			lane_open[i] = true;
		}

		if (fLog) fprintf(fLog, "ego lane %d target lane %d\n", ego_lane, target_lane);
		for (auto &p : sensor_fusion_cars)
		{
			auto &other = p.second;
			int lane = other.lane;
			double s = other.predicted_s(delta_t0);
			if (s > ego_s)
			{
				if (s < next_s_in_lane[lane])
				{
					next_s_in_lane[lane] = s;
					double max_dist_for_speed = 200; // m
					if (s - ego_s < max_dist_for_speed)
					{
						double cutoff_dist_for_speed = 100; // m
						// over 200 meters it's unreliable, but we want to avoid sudden jumps in scores just because a car appears 200 m ahead
						// so between 100 and 200 meters, the speed is handed gradually
						
						int speed = other.vs;
						if (speed > max_speed) speed = max_speed;
						if (s - ego_s > cutoff_dist_for_speed)
						{
							// linear fit in the [100m,200m] range for [speed, max_speed]:
							speed = speed + (max_speed - speed)*(s - ego_s - cutoff_dist_for_speed) / (max_dist_for_speed - cutoff_dist_for_speed);
						}
						lane_speed[lane] = speed;
					}
				}
			}
			double additional_safety_distance_for_lane_change = 2; // added treshold to avoid aborted lane change
			if (target_lane == lane) additional_safety_distance_for_lane_change = 0; // already decided to change lane, work with normal safety
			double min_dist = car_length + safety_distance + additional_safety_distance_for_lane_change;
			if (fabs(ego_s - s) < min_dist)
			{
				if (fLog && lane_open[lane]) if (fLog) fprintf(fLog, "lane %d closed, car %d at %.2f inside range %.2f\n", lane, other.id, s, min_dist);
				lane_open[lane] = false;
			}
			// check if we could potentially break if it's ahead, or the other car won't crash into us if it's behind:
			if (s > ego_s && other.vs < ego_vs)
			{
				double car_dist = s - ego_s - car_length - safety_distance - additional_safety_distance_for_lane_change;
				double speed_diff = ego_vs - other.vs;
				double decelerate_time = speed_diff / relaxed_acc;
				double decelerate_dist = ego_vs * decelerate_time - speed_diff / 2 * decelerate_time;
				if (fLog) fprintf(fLog, "checking car %d ahead in lane %d: dist %.2f speed %.2f diff %.2f decelerate time %.2f dist %.2f\n",
					other.id, other.lane, car_dist, other.vs, speed_diff, decelerate_time, decelerate_dist);

				if (car_dist < decelerate_dist)
				{
					if (fLog && lane_open[lane]) if (fLog) fprintf(fLog, "lane %d closed, can't decelerate to %.2f in %.2fm, dist %.2f\n", lane, other.vs, decelerate_dist, car_dist);
					lane_open[lane] = false;
				}
			}
			if (s<ego_s && other.vs>ego_vs && s+50>ego_s) // cars far away can have funny speed values
			{
				double car_dist = ego_s - s - car_length - safety_distance - additional_safety_distance_for_lane_change;
				double speed_diff = other.vs - ego_vs;
				double decelerate_time = speed_diff / relaxed_acc;
				if (target_lane == ego_lane) // currently not switching lane, or returning
					decelerate_time += 2; // we give it more space in case it won't break in time, but once we start lane change, try to finish it
				double min_dist = speed_diff * decelerate_time;
				if (fLog) fprintf(fLog, "checking car %d behind in lane %d: dist %.2f speed %.2f diff %.2f decelerate time %.2f dist %.2f\n",
					other.id, other.lane, car_dist, other.vs, speed_diff, decelerate_time, min_dist);
				if (car_dist < min_dist)
				{
					if (fLog && lane_open[lane]) if (fLog) fprintf(fLog, "lane %d closed, approaching car speed %.2f need %.2fm, dist %.2f\n", lane, other.vs, min_dist, car_dist);
					lane_open[lane] = false;
				}
			}
		}

		int best_target_lane = ego_lane;
		double best_score = 0;
		double lane_score[NUM_LANES];
		for (int lane = 0; lane < NUM_LANES; lane++)
		{
			lane_score[lane] = 0;
			if (lane != ego_lane && !lane_open[lane]) continue;
			double speed_score = std::min(lane_speed[lane] / max_speed, 1.0);
			// too much noise in this:
			//double ego_d_future = ego_d+ego_vd * 2;
			//double lane_center_offset = map.get_lane_center_offset(lane);
			//double distance_score = std::max(0.0, 1 - fabs(ego_d_future - lane_center_offset)/8); // roughly 0 for distant, 0.5 for next, 1 for same
			double distance_score = 1 - fabs(target_lane - lane) / 2;
			//if (target_lane != ego_lane && lane != target_lane) distance_score = 0; // if we already started lane change, prefer the target to minimize change of mind (have to finish maneuver in 3 sec)
			double free_score = std::min(1.0, next_s_in_lane[lane] / 100);
			if (test_fast_lane_change) distance_score = 0;
			double total_score = speed_score + distance_score / 2 + free_score;
			if (fLog) fprintf(fLog, "lane %d: %.2f,%.2f,%.2f -> %.2f\n", lane, speed_score, distance_score, free_score, total_score);
			lane_score[lane] = total_score;
			if (total_score > best_score)
			{
				best_score = total_score;
				best_target_lane = lane;
			}
		}
		if (console_log) printf(" lane_scores %.2f %.2f %.2f best %d->%d", lane_score[0], lane_score[1], lane_score[2], ego_lane, best_target_lane);
		if (abs(ego_lane - best_target_lane) > 1)
		{
			// check next lane instead
			int next_lane = best_target_lane > ego_lane ? ego_lane + 1 : ego_lane - 1;
			if (!lane_open[next_lane]) target_lane = ego_lane;
			else target_lane = next_lane;
		}
		else
		{
			target_lane = best_target_lane;
		}
		return target_lane;
	}
};

class SpeedController
{
public:
	double start_speed;
	double target_speed;
	double target_time;
	double time_shift; // when acceleration is maximized, shift the linear function to match the calculated speed value
	SpeedController(double ego_speed)
	{
		time_shift = 0;
		start_speed = ego_speed;
		target_speed = max_speed;
		double delta_to_max_speed = fabs(ego_speed - max_speed);
		target_time = delta_to_max_speed / relaxed_acc;
	}
	double get_speed(double current_t)
	{
		current_t -= time_shift;
		if (current_t < 0)
			current_t = 0;
		double speed;
		if (current_t > target_time) speed = target_speed;
		else speed = start_speed + (target_speed - start_speed)*current_t / target_time;
		return speed;
	}
	void add_limit_breakpoint(double new_target_speed, double new_target_time)
	{
		// make sure we don't divide by 0:
		double time_mod = std::max(target_time, 0.02);
		double new_time_mod = std::max(new_target_time, 0.02);
		double current_grade = (target_speed - start_speed) / time_mod;
		double new_grade = (new_target_speed - start_speed) / new_time_mod;
		if (fLog)
		{
			if (target_speed != max_speed && new_target_speed != max_speed)
			{
				fprintf(fLog, "SpeedController start %.2f prev target %.2f at %.2f grade %.2f, new target at %.2f at %.2f grade %.2f\n",
					start_speed, target_speed, target_time, current_grade, new_target_speed, new_target_time, new_grade);
			}
		}
		if (new_grade < current_grade)
		{
			target_speed =new_target_speed;
			target_time = new_target_time;
		}
	}
	void override_speed(double current_t, double speed)
	{
		// modify time_shift to keep the same slope and target speed, but with speed at current_t
		if (current_t > target_time) return; // zero grade from next step, and also prevents divide by zero
		//double current_grade = (target_speed - start_speed) / target_time;
		// get modified current_t from from speed = start_speed + (target_speed - start_speed)*current_t / target_time
		if (fabs(target_speed - start_speed) < EPSILON) return;
		// speed - start_speed = (target_speed - start_speed)*current_t / target_time
		// target_time*(speed-start_speed)/(target_speed-start_speed) = current_t
		double mod_current_t = target_time * (speed - start_speed) / (target_speed - start_speed);
		// we reported speed at mod_current_t (like at 1 sec), but we'll report that value at current_t from now on (like 1.4 sec)
		// then we need to have time_shift = 0.4 sec, so that get_speed(1.4) would return the value for (1.4-0.4) sec
		time_shift = current_t - mod_current_t;
	}
};

class TrajectoryBuilder
{
public:
	double total_dist = 0;
	vector<Point > control_points;
	void add_control_point(Point p)
	{
		if (control_points.size() > 0)
		{
			auto prev_p = control_points.back();
			total_dist += distance(prev_p.x, prev_p.y, p.x, p.y);
		}
		control_points.push_back(p);
	}
	
	vector<Point> build(vector<Point> &prev_trajectory,
		double ego_x, double ego_y, double ego_yaw, // fallback for init if not enough points in prev_trajectory
		int ego_lane, int target_lane, double ego_d, double ego_vd,
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
		vector<Point> pre_trajectory_spline_control_points;
		Point pre_trajectory_start_pos(ego_x, ego_y);


		int path_size = prev_trajectory.size();
		if (path_size == 0) {
			pos_x = ego_x;
			pos_y = ego_y;
			angle = deg2rad(ego_yaw);
		}
		else 
		{
			pre_trajectory_start_pos = prev_trajectory[0];
			pos_x = prev_trajectory[path_size - 1].x;
			pos_y = prev_trajectory[path_size - 1].y;
			if (path_size == 1)
			{
				angle = deg2rad(ego_yaw);
			}
			else
			{
				double pos_x2 = prev_trajectory[path_size - 2].x;;
				double pos_y2 = prev_trajectory[path_size - 2].y;
				Point v(pos_x - pos_x2, pos_y - pos_y2);
				if (v.lengthsq() < EPSILON)
					angle = deg2rad(ego_yaw);
				else
				{
					angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
				}
			}
		}
		bool pre_trajectory_from_angle = true;
		double pre_trajectory_point_dist = 1.0;
		Point v_pre; // normalized vector pointing backwards
		if (prev_trajectory.size() >= 2)
		{
			v_pre = prev_trajectory[0] - prev_trajectory[1];
			double v_pre_len = v_pre.length();
			if (v_pre_len > EPSILON)
			{
				pre_trajectory_from_angle = false;
				v_pre.x /= v_pre_len;
				v_pre.y /= v_pre_len;
			}
		}
		if (pre_trajectory_from_angle)
		{
			v_pre.x = -cos(angle);
			v_pre.y = -sin(angle);
		}
		
		for (int i = 3; i >=1; i--)
		{
			pre_trajectory_spline_control_points.push_back(
				Point(pre_trajectory_start_pos.x + v_pre.x*i*pre_trajectory_point_dist, 
					pre_trajectory_start_pos.y + v_pre.y*i*pre_trajectory_point_dist));
		}
		//printf(" ego yaw %.2f angle %.2f", deg2rad(ego_yaw), angle);
		add_control_point(Point(pos_x, pos_y)); // startup

		double min_control_point_dist = speed_controller.start_speed * 1.5;
		min_control_point_dist = std::max(min_control_point_dist, 5.0); // at least 5 m
		//double first_control_point_dist = min_control_point_dist; // first control point is to switch lane or return to center line
		double contol_point_start_s = 0;
		
		//if (ego_lane != target_lane)
		{
			// we need more distance if heading the opposite way, and less if car is already aligned
			// about half the distance if aligned, because we need to steer back to straight
			double d_diff = map.get_lane_center_offset(target_lane) - ego_d;

			// we have d-wise speed and constant d-wise acceleration
			// first, we accelerate up, then slow down
			// a normal lane change is 2 seconds, so we accelerate for 1 seconds, then slow down, and the result distance is 4 meters
			// first half or second half: 2 = a/2*1*1 => d-wise acceleration is 4 m/s^2
			// depending on d_diff and ego_vd, we determine if we need to still accelerate or just slow down is enough:
			// time to stop: t = ego_v/a
			// d_diff_max = ego_v*t-a/2*t*t and 0=ego_v-a*t
			// or simply: d_diff_max = ego_v*t/2 = ego_v*ego_v/a/2
			double d_acc = 4; // d-wise acceleration, m/s^2
			bool slow_down_phase = false;
			double lane_switch_time = 2.0; // calculate below
			if ((ego_vd<0) == (d_diff<0))
			{
				double d_diff_max = ego_vd * ego_vd / d_acc / 2;
				if (d_diff_max > abs(d_diff))
				{
					slow_down_phase = true;
					lane_switch_time = fabs(ego_vd) / d_acc;
				}
			}
			// if we're not slowing down, we accelerate first, then decellerate
			if (!slow_down_phase)
			{
				// we'll reach peak_vd after t1 and decelerate from there back to 0:
				// decelerate_phase_dist = peak_vd*peak_vd/a/2
				// accelerate_phase_dist = ego_vd*t+a/2*t*t
				// acceleration time: t=(peak_vd-ego_vd)/a
				// accelerate_phase_dist + decelerate_phase_dist = d_diff
				// so:
				// ego_vd*t+a/2*t*t + peak_vd*peak_vd/a/2 = d_diff
				// ego_vd*(peak_vd-ego_vd)/a+a/2*(peak_vd-ego_vd)/a*(peak_vd-ego_vd)/a + peak_vd*peak_vd/a/2  = d_diff
				// ego_vd*(peak_vd-ego_vd)/a+(peak_vd-ego_vd)*(peak_vd-ego_vd)/a/2 + peak_vd*peak_vd/a/2  = d_diff
				// 2*ego_vd*(peak_vd-ego_vd)+(peak_vd-ego_vd)*(peak_vd-ego_vd) + peak_vd*peak_vd  = d_diff * a * 2
				// 2ep-2ee+pp-2pe+ee+pp=d_diff
				// 2pp-ee = d_diff*a*2
				// peak_vd = sqrt(d_diff*a+ego_vd*ego_vd/2)
				// checking with normal change, ego_vd=0, d_diff=4, a=4:
				// peak_vd = sqrt(4*4) = 4 m/s, average speed is 2 m/s, lane change completes in 2 seconds
				// checking halfway lane change, ego_vd=4, d_diff=2, a=4:
				// peak_vd = sqrt(2*4+4*4/2) = 4 m/s lane change completes in 1 second

				// maybe easier another calculation:
				// accelerate phase from 0: peak_vd*peak_vd/a/2
				// if ego_vd has the same sign:
				// accelerate phase from ego_vd: peak_vd*peak_vd/a/2-ego_vd*ego_vd/a/2
				// total dist: d_diff = peak_vd*peak_vd/a-ego_vd*ego_vd/a/2
				// peak_vd = sqrt(d_diff*a+ego_vd*ego_vd/2)
				// need the time for the lane change: t_acc=(peak_vd-ego_vd)/a
				// t_decelerate = peak_vd/a
				// total time: (peak_vd*2-ego_vd)/a
				// if ego_v is negative, we first accelerate back to 0, we'll calculate with relative ego_v, it's relative to d_diff direction
				double rel_ego_vd = ego_vd;
				if (d_diff < 0) rel_ego_vd *= -1;
				double abs_d_diff = fabs(d_diff);

				double peak_vd = sqrt(abs_d_diff * d_acc + rel_ego_vd * rel_ego_vd / 2);
				lane_switch_time = (peak_vd * 2 - rel_ego_vd) / d_acc;
				if (lane_switch_time < 0)
				{
					// should be handled in slow_down_phase
					printf("lane_switch_time(%.2f) < 0\n", lane_switch_time);
				}
			}

			// we need 0 speed and acceleration when we arrive at the center of target lane
			// 0 = ego_v + t*acc
			// d_diff = ego_v*t + acc/2*t^2

			//double vd_diff = fabs(ego_vd); // need to straighten up
			//double disalignment = fabs(ego_vd - d_diff); // comparing d in 1 sec into future to see disalignment
			// if disalignment is 0, like middle of lane switch, we need about 1 sec for 2 meter,
			// if 4 meter, we just start the lane switch, need about 2 seconds
			// if 6 meter, we might in the middle of an opposite lane change, need about 3 sec
			//double lane_switch_time = 2.0 * d_diff/4.0+disalignment; // 2 seconds / 4 meter
			double min_lane_switch_dist = 15.0;
			double speed = speed_controller.start_speed;
			double dist = speed * lane_switch_time;
			if (dist < min_lane_switch_dist) dist = min_lane_switch_dist;
			if (dist > 50) dist = 50; // for offroad
			contol_point_start_s = dist;
			if (fLog) fprintf(fLog, "lane switch time = %.2f raw %.2f dist=%.2f increment %.2f\n", lane_switch_time, speed * lane_switch_time, dist, min_control_point_dist);
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
			
			Point next_pt = map.get_lane_pos(contol_point_start_s, target_lane, next_wp, waypoint_dist);
			/* This didn't work out:
			if (waypoint_dist < min_control_point_dist)
			{
				// if the next waypoint is closer than double min_control_point_dist, we skip the next control point,
				// and add the waypoint instead to avoid cutting corners, or adding too dense control points which would make it jerk
				double tmp_dist;
				int tmp_wp;
				//next_pt = map.get_lane_pos(contol_point_start_s- min_control_point_dist/2+waypoint_dist/2, target_lane, tmp_wp, tmp_dist);
				//add_control_point(next_pt);
				add_control_point(map.get_waypoint(next_wp).lane_center[target_lane]);
				contol_point_start_s += waypoint_dist;
			}
			else*/
			{
				add_control_point(next_pt);
			}
			if (total_dist > 50 && control_points.size() > 2) break;
			contol_point_start_s += min_control_point_dist;
		}

		if (fLog)
		{
			fprintf(fLog, "first control dist %.2f\n", (Point(pos_x, pos_y) - control_points[1]).length());
			fprintf(fLog, "pre_control=[");
			for (int i = 0; i < pre_trajectory_spline_control_points.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", pre_trajectory_spline_control_points[i].x, pre_trajectory_spline_control_points[i].y);
			fprintf(fLog, "]\n");
			fprintf(fLog, "prev_trajectory=[");
			for (int i = 0; i < prev_trajectory.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", prev_trajectory[i].x, prev_trajectory[i].y);
			fprintf(fLog, "]\n");
			fprintf(fLog, "control_points=[");
			for (int i = 0; i < control_points.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", control_points[i].x, control_points[i].y);
			fprintf(fLog, "]\n");
		}


		// transform control_points for spline
		double cos_angle = cos(-angle);
		double sin_angle = sin(-angle);
		Point transform_center(pos_x, pos_y);
		for (int i = 0; i < (int)control_points.size(); i++)
		{
			Point &p = control_points[i];
			p = p - transform_center;
			double tx = p.x*cos_angle - p.y*sin_angle;
			double ty = p.x*sin_angle + p.y*cos_angle;
			p.x = tx;
			p.y = ty;
		}
		vector<double> wp_x, wp_y;
		for (int i = 0; i < (int)pre_trajectory_spline_control_points.size(); i++)
		{
			Point &p = pre_trajectory_spline_control_points[i];
			p = p - transform_center;
			double tx = p.x*cos_angle - p.y*sin_angle;
			double ty = p.x*sin_angle + p.y*cos_angle;
			wp_x.push_back(tx);
			wp_y.push_back(ty);
		}
		for (int i = 0; i < int(prev_trajectory.size()) - 1; i++)
		{
			Point &p = prev_trajectory[i];
			Point t = p - transform_center;
			double tx = t.x*cos_angle - t.y*sin_angle;
			double ty = t.x*sin_angle + t.y*cos_angle;

			wp_x.push_back(tx);
			wp_y.push_back(ty);
		}
		int min_control_point_count = wp_x.size(); // so far we only added previous points for shaping
		pos_x = 0;
		pos_y = 0;
		double transform_angle = angle;
		cos_angle = cos(transform_angle);
		sin_angle = sin(transform_angle);

		tk::spline spline;
		
		for (int i = 0; i < (int)control_points.size(); i++)
		{
			wp_x.push_back(control_points[i].x);
			wp_y.push_back(control_points[i].y);
		}

		for (int i = 1; i < (int)wp_x.size(); i++)
		{
			if (wp_x[i] <= wp_x[i - 1])
			{
				printf("\nspline input error\n");
				if (fLog) fprintf(fLog, "spline input error at idx %d x values %.4f<%.4f\n", i, wp_x[i], wp_x[i - 1]);
				wp_x.resize(i);
				wp_y.resize(i);
				break;
			}
		}

		double dist_total = 0;
		double current_t = 0.02;

		if (wp_x.size() < 3 || wp_x.size()<= min_control_point_count || fabs(ego_d)>20) // offroad is not good with splines
		{
			// not enough control points, can't generate spline, generate trajectory based on angle
			double speed = speed_controller.get_speed(current_t);
			double current_angle = 0; // working in transformed space
			int next_control_point = 1;
			for (; result_points.size() < 50 && next_control_point< control_points.size();)
			{
				double dist_step = speed / 50;
				Point next_pt = control_points[next_control_point];
				Point current_p = Point(pos_x, pos_y);
				Point next_delta = (next_pt - current_p);
				double control_point_dist = next_delta.length();
				if (control_point_dist < 5) // approach 5 meter radius, we may not reach the control point exactly due to angle constraints
				{
					next_control_point++;
					continue;
				}
				current_t += 0.02;
				//if (speed > EPSILON)
				{
					double next_control_point_angle = atan2(next_delta.y, next_delta.x);
					double angle_diff = fmod(next_control_point_angle - current_angle + 3 * pi(), 2 * pi()) - pi();
					// v=sqrt(radius)*const, so radius = v^2*K. Also, turn radius is minimum 10 meters.
					// Tried in simulator, K is about 1/centrifugal force.
					double max_acceleration = 4; // m/s^2
					double min_radius = std::max(10.0, speed*speed / max_acceleration); // limit to 4 m/s^2
					// with this speed, we complete the 2pi circle t= min_radius*2*pi/speed
					// so in sec/rad: min_radius/speed
					// we need rad/sec: speed/min_radius
					double rad_per_sec = speed / min_radius;
					double max_angle_step =  rad_per_sec / 50;
					if (fabs(angle_diff) > max_angle_step)
					{
						if (angle_diff > 0) current_angle += max_angle_step;
						else current_angle -= max_angle_step;
					}
					else
					{
						current_angle += angle_diff;
					}
				}
				//pos_x += next_delta.x*dist_step / control_point_dist;
				//pos_y += next_delta.y*dist_step / control_point_dist;
				pos_x += cos(current_angle)*dist_step;
				pos_y += sin(current_angle)*dist_step;
				// transform back:
				Point tp;
				tp.x = pos_x * cos_angle - pos_y * sin_angle;
				tp.y = pos_x * sin_angle + pos_y * cos_angle;
				tp = tp + transform_center;
				result_points.push_back(tp);
			}
		}
		else
		{
			spline.set_points(wp_x, wp_y);

			//double speed = max_speed;
			//for (int i = 1; i < (int)waypoints.size(); i++)
			double current_sp_arg = 0;
			double prev_speed = speed_controller.start_speed;
			double prev_angle = 0;
			while (current_sp_arg < 50) // 50 meters max, will stop at 50 points anyway
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
				if (dist + EPSILON < dist_step)
				{
					printf("\nspline warning, dist: %.4f, step: %.4f, pos %.2f,%.2f to %.2f,%.2f\n", dist, dist_step, pos_x, pos_y, x, y);
					if (fLog) fprintf(fLog, "spline warning, dist: %.4f, step: %.4f, pos %.2f,%.2f to %.2f,%.2f\n", dist, dist_step, pos_x, pos_y, x, y);
				}

				// we now have a general direction, just have to update speed or decrease angle change to keep below target acceleration
				double acceleration = fabs(speed - prev_speed)*50;
				double angle_step = atan2(y - pos_y, x - pos_x);
				double angle_diff = fmod(angle_step - prev_angle + 3 * pi(), 2 * pi()) - pi();
				// from: rad_per_sec = speed / min_radius and angle_step = rad_per_sec / 50
				//double turn_radius = speed / 50 / fabs(angle_diff);
				// from radius = speed^2/acc -> centrifugal_acceleration = speed * speed / turn_radius
				// centrifugal_acceleration = speed * speed / speed * 50 * angle_diff;
				double centrifugal_acceleration = speed * 50 * fabs(angle_diff);
				if (acceleration + centrifugal_acceleration > maximum_acc)
				{
					// need to intervent, either by smoothing the angle or reducing acceleration
					// reduced acceleration is more simple:
					if (speed > prev_speed) // limit to acceleration and not braking
					{
						double new_acceleration = maximum_acc - centrifugal_acceleration;
						if (new_acceleration < 0)
						{
							if (fLog) fprintf(fLog, "accT too high, accN=%.2f accT=%.2f\n", centrifugal_acceleration, acceleration);
							new_acceleration = 0;
						}

						double new_speed;
						if (speed > prev_speed)
						{
							new_speed = prev_speed + new_acceleration / 50;
						}
						else 
						{
							// if still adjusting braking
							new_speed = prev_speed - new_acceleration / 50;
						}
						if (fLog) fprintf(fLog, "acceleration override, accN=%.2f accT=%.2f -> %.2f speed %.2f -> %.2f\n", centrifugal_acceleration, acceleration, new_acceleration, speed, new_speed);
						// reset speed control to new value: (todo)
						speed_controller.override_speed(current_t, new_speed);
						speed = new_speed;
						speed_controller.target_time += 0.02;
						dist_step = speed / 50;
						acceleration = new_acceleration;
					}
					if (acceleration + centrifugal_acceleration > maximum_acc)
					{
						double new_centrifugal_acceleration = maximum_acc - acceleration;
						
						if (new_centrifugal_acceleration < 0)
						{
							if (fLog) fprintf(fLog, "accN too high (overbraking?), accN=%.2f accT=%.2f\n", centrifugal_acceleration, acceleration);
							new_centrifugal_acceleration = 0;
						}
						if (fLog) fprintf(fLog, "adjusting curvature for pt %d, accN=%.2f -> %.2f accT=%.2f\n", result_points.size(), centrifugal_acceleration, new_centrifugal_acceleration, acceleration);
						// from double centrifugal_acceleration = speed * 50 * fabs(angle_diff);
						double new_angle_diff = new_centrifugal_acceleration / speed / 50;
						if (angle_diff < 0) new_angle_diff *= -1;
						double new_angle = angle_step + angle_diff;
						// need to recalculate transform center, 
						// rotate around pos_x and pos_y with new_angle_diff, 
						// and calculate new transform center with transform_angle+new_angle_diff going through pos_x, pos_y
						// for that, we need to rotate transform_center around transformed(pos_x,posy) with new_angle_diff
						// transform back pos:
						Point tp;
						tp.x = pos_x * cos_angle - pos_y * sin_angle;
						tp.y = pos_x * sin_angle + pos_y * cos_angle;
						tp = tp + transform_center;
						// rotate transform center around tp:
						Point vect = transform_center - tp;
						Point rotated_vect;
						rotated_vect.x = vect.x*cos(new_angle_diff) - vect.y*sin(new_angle_diff);
						rotated_vect.y = vect.x*sin(new_angle_diff) + vect.y*cos(new_angle_diff);
						transform_center = tp + rotated_vect;
						transform_angle += new_angle_diff;
						cos_angle = cos(transform_angle);
						sin_angle = sin(transform_angle);
						// just testing, pos_x, posy should be the same:
						Point tp_save = tp;
						tp.x = pos_x * cos_angle - pos_y * sin_angle;
						tp.y = pos_x * sin_angle + pos_y * cos_angle;
						tp = tp + transform_center;
						if ((tp_save - tp).lengthsq() > EPSILON)
						{
							if (fLog) fprintf(fLog, "transform calculation error, angle: %.2f around pos %.2f, %.2f got back %.2f,%.2f\n",
								new_angle_diff, tp_save.x, tp_save.y, tp.x, tp.y);
						}
					}
				}
				
				//if (current_t == 0.02) printf(" next sp %.2f", speed);
				current_t += 0.02;
				prev_speed = speed;
				prev_angle = angle_step;
				double sp_step = (x - pos_x)*dist_step / dist;
				pos_y += (y - pos_y)*dist_step / dist;
				current_sp_arg += sp_step;
				//extra precise spline:
				//y = spline(current_sp_arg);
				//pos_y = y;
				pos_x += sp_step;
				
				// transform back:
				Point tp;
				tp.x = pos_x * cos_angle - pos_y * sin_angle;
				tp.y = pos_x * sin_angle + pos_y * cos_angle;
				tp = tp + transform_center;
				result_points.push_back(tp);
				if (result_points.size() >= 50) break;
			}
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

class LimitSpeed // limits the speed of the ego vehicle based when following single car
{
public:
	double target_speed;
	double target_time;
	bool can_accelerate = true;

	void maximize_acc(double ego_speed, double acc)
	{
		double delta_v = target_speed - ego_speed;
		double min_time = fabs(delta_v) / acc;
		if (target_time < min_time)
		{
			target_time = min_time;
		}
	}
	void calculate(Car &follow_car, double next_s, double ego_s, double ego_speed, double ego_acc)
	{
		target_speed = max_speed;
		double delta_to_max_speed = fabs(ego_speed - max_speed);
		target_time = delta_to_max_speed / relaxed_acc;

		double follow_car_distance = next_s - ego_s - car_length;
		if (follow_car_distance < 0)
		{
			printf("detected collision\n");
			follow_car_distance = 0; // collision already
		}
		double follow_car_speed = sqrt(follow_car.vx*follow_car.vx + follow_car.vy*follow_car.vy);
		
		const char *code = "FREEFLOW";
		// check emergency breaking, safe distance
		if (ego_speed > follow_car_speed)
		{
			double acc = relaxed_acc;
			if (ego_acc < 0) acc = min_relaxed_acc_while_braking;
			double delta_v = ego_speed - follow_car_speed;
			double delta_t = delta_v / acc;
			double delta_d = ego_speed * delta_t - delta_v / 2 * delta_t;
			double max_dist = follow_car_distance - safety_distance;

			if (delta_d > max_dist)
			{
				target_speed = follow_car_speed;
				// we have max_dist to decelerate, so from max_dist = (ego_speed - delta_v / 2)*target_speed_time:
				target_time = max_dist / (ego_speed - delta_v / 2);
				// if max accel is reached, use up the safety distance:
				if (target_time < EPSILON || delta_v / target_time>maximum_acc)
				{
					code = "MAXBRAKE";
					if (console_log) printf(" MAXBRAKE");
					if (fLog) fprintf(fLog, "maxbrake: need %.2f have %.2f delta v %.2f, target time %.2f acc %.2f\n",
						delta_d, max_dist, delta_v, target_time, delta_v / target_time);
					target_time = delta_v / maximum_acc;
				}
				else
				{
					code = "BRAKE";
					if (console_log) printf(" BRAKE");
					if (fLog) fprintf(fLog, "normalbrake: need %.2f have %.2f delta v %.2f, target time %.2f acc %.2f\n",
						delta_d, max_dist, delta_v, target_time, delta_v / target_time);
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

				target_speed = follow_car_speed - excess_distance / t_to_reach_optimal_distance;
				target_time = t_to_reach_optimal_distance;
				maximize_acc(ego_speed, relaxed_acc);
				can_accelerate = false;
				code = "ADJUST";
				if (console_log) printf(" ADJUST");

				/*
				vector<double> start = { ego_s, ego_speed, ego_acc };
				vector<double> end = { follow_car_predicted_pos - car_length, follow_car_speed, 0 };
				in_lane_jmt = JMT(start, end, t_to_reach_optimal_distance);*/
			}
			else if (ego_s + car_length + keep_distance + keep_distance_leeway > next_s)
			{
				// we maintain current speed, don't try to exacly maintain safety distance
				target_speed = follow_car_speed;
				target_time = 1.0;
				maximize_acc(ego_speed, relaxed_acc);
				if (console_log) printf(" KEEP");
				code = "KEEP";
			}
		}
		if (fLog) fprintf(fLog, "limitspeed %s for car %d in lane %d dist %.2f: %.2f in %.2f\n", code, follow_car.id, follow_car.lane, next_s - ego_s, target_speed, target_time);
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

  if (need_log)
  {
	  fLog = fopen("trajectory.log", "wt");
	  fprintf(fLog, "wpmap=[");
	  for (int i = 0; i < map_waypoints_x.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", map_waypoints_x[i], map_waypoints_y[i]);
	  fprintf(fLog, "]\n");
	  for (int lane = 0; lane < 3; lane++)
	  {
		  fprintf(fLog, "lane%d=[", lane);
		  for (int i = 0; i < map.waypoints.size(); i++) fprintf(fLog, "%s[%.4f,%.4f]", i == 0 ? "" : ",", map.waypoints[i].lane_center[lane].x, map.waypoints[i].lane_center[lane].y);
		  fprintf(fLog, "]\n");
	  }
  }

  int frame_count = 0;
  bool lane_change_complete = true;

  h.onMessage([&map, &prev_car_speed, &prev_speed_valid, &prev_car_acc, &prev_acc_valid, &sensor_fusion_cars, &target_lane, &lane_change_complete, &frame_count]
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
          double ego_speed = j[1]["speed"];
		  ego_speed /= 2.237;
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

		  Point ego_speed_vector;
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
			  ego_speed_vector = prev_trajectory[3] - prev_trajectory[2];
			  double v3 = ego_speed_vector.length();
			  ego_acc = (v3 - v2) * 50;
			  double prev_acc = (v2 - v1) * 50;
			  ego_jerk = (ego_acc - prev_acc) * 50;
			  ego_speed = v3 * 50;
			  ego_speed_vector.x *= 50;
			  ego_speed_vector.y *= 50;
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
		  if (sparse_console_log)
		  {
			  console_log = !(frame_count % 50);
		  }
		  double ego_vs, ego_vd;
		  map.project_speed(ego_speed_vector, map.reference_waypoint_id, &ego_vs, &ego_vd);
		  if (console_log) printf("%04d v% 4.2f a% 4.2f d% 4.2f vd % 4.2f", frame_count, ego_speed, ego_acc, ego_d, ego_vd);
		  
		  if (fLog) fprintf(fLog, "#%04d wp%d v%.2f a%.2f d%.2f\n", frame_count, map.reference_waypoint_id, ego_speed, ego_acc, ego_d);
		  frame_count++;

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
				  map.project_speed(Point(car.vx, car.vy), next_wp_id, &car.vs, &car.vd);
				  //car.s = car_data[5];
				  //car.d = car_data[6];
				  if (fLog) fprintf(fLog, "car %d at (%.2f,%.2f) speed (%.2f,%.2f) lane %d s %.2f d %.2f vs %.2f vd %.2f\n",
					  car.id, car.x, car.y, car.vx, car.vy, car.lane, car.s, car.d, car.vs, car.vd);

			  }
		  }

		  LaneChangePlanner lane_change_planner;
		  target_lane = lane_change_planner.calculate_target_lane(sensor_fusion_cars, ego_lane, target_lane, ego_s, ego_vs, delta_t0);
		  // check if we're off center to the other direction:
		  if (target_lane != ego_lane)
		  {
			  double d_of_target_lane = map.get_lane_center_offset(target_lane);
			  double predict_t = 1.0; //sec, see where the car would be
			  double lane_d_diff = fabs(ego_vd*predict_t + ego_d - d_of_target_lane);
			  if (fLog) fprintf(fLog, "checking lane switch to %d: diff %.2f ego d %.2f vd %.2f\n", target_lane, lane_d_diff, ego_d, ego_vd);
			  if (lane_d_diff > 6.0)
			  {
				  if (fLog) fprintf(fLog, "target lane too far: %.2f ego d %.2f vd %.2f, keep in lane\n", lane_d_diff, ego_d, ego_vd);
				  target_lane = ego_lane;
			  }
		  }
		  if (target_lane != ego_lane) lane_change_complete = false;

		  if (target_lane == ego_lane && !lane_change_complete)
		  {
			  double d_of_target_lane = map.get_lane_center_offset(target_lane);
			  double lane_d_diff = fabs(ego_d - d_of_target_lane);
			  if (lane_d_diff < 0.5)
			  {
				  if (fLog) fprintf(fLog, "lane change complete to lane %d, d=%.2f\n", ego_lane, ego_d);
				  lane_change_complete = true;
			  }
		  }

		  int next_car_id = -1;
		  double next_s = 0;
		  int next_car_in_target_lane = -1;
		  double next_s_in_target_lane = 0;
		  double d_of_target_lane = map.get_lane_center_offset(target_lane);
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

			  if (s0 >= ego_s - car_length - safety_distance && fabs(d0 - d_of_target_lane) < 3)
			  {
				  if (next_car_in_target_lane == -1 || next_s_in_target_lane > s0)
				  {
					  next_car_in_target_lane = other.id;
					  next_s_in_target_lane = s0;
				  }
			  }
		  }
		  if (next_car_in_target_lane == next_car_id) next_car_in_target_lane = -1;// only check once
		  
		  /*
		  if (next_car_id != -1)
		  {
			  auto &car = sensor_fusion_cars[next_car_id];
			  printf(" next (%.2f, %.2f) v (%.2f, %.2f) s_dist: %.2f lane %d",
				  car.s, car.d, car.vs, car.vd, next_s-ego_s, car.lane);
		  }*/
		  //vector<double> in_lane_jmt;
		  //double in_lane_jmt_max_time = 0;
		  SpeedController speed_controller(ego_speed);
		  
		  // during lane change, we may have to adjust speed on 2 cars
		  if (next_car_id != -1)
		  {
			  LimitSpeed limit_speed_in_lane;
			  auto &follow_car = sensor_fusion_cars[next_car_id];
			  limit_speed_in_lane.calculate(follow_car, next_s, ego_s, ego_speed, ego_acc);
			  speed_controller.add_limit_breakpoint(limit_speed_in_lane.target_speed, limit_speed_in_lane.target_time);
		  }
		  if (next_car_in_target_lane != -1)
		  {
			  LimitSpeed limit_speed_in_target_lane;
			  auto &follow_car = sensor_fusion_cars[next_car_in_target_lane];
			  limit_speed_in_target_lane.calculate(follow_car, next_s_in_target_lane, ego_s, ego_speed, ego_acc);
			  speed_controller.add_limit_breakpoint(limit_speed_in_target_lane.target_speed, limit_speed_in_target_lane.target_time);
		  }

		  //printf(" target %.2f @ %.2f", speed_controller.target_speed, speed_controller.target_time);
		  //if (in_lane_jmt.empty()) 
		  //else printf(" JMT start speed %.2f end speed %.2f (%.2fs)", get_jmt_speed(in_lane_jmt, 0), get_jmt_speed(in_lane_jmt, in_lane_jmt_max_time), in_lane_jmt_max_time);

          json msgJson;

		  TrajectoryBuilder trajectory;
		  
		  vector<Point> refined_trajectory = trajectory.build(prev_trajectory, ego_x, ego_y, ego_yaw, ego_lane, target_lane, ego_d, ego_vd, map, speed_controller);

		  vector<double> next_x_vals(refined_trajectory.size());
		  vector<double> next_y_vals(refined_trajectory.size());

		  for (int i = 0; i < (int)refined_trajectory.size(); i++)
		  {
			  next_x_vals[i] = refined_trajectory[i].x;
			  next_y_vals[i] = refined_trajectory[i].y;
		  }

		  if (fLog) fflush(fLog);
		  if (console_log) printf("\n");
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