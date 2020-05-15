#define _USE_MATH_DEFINES
#include <uWS/uWS.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "spline.h"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;

struct Point
{
	double x, y;
	Point() { x = y = 0; }
	Point(double _x, double _y) :x(_x), y(_y) {}
};

double distancesq_pt_pt(Point p, Point A)
{
	return (p.x - A.x)*(p.x - A.x) + (p.y - A.y)*(p.y - A.y);
}
double distancesq_pt_seg(Point p, Point A, Point B,
	double &_rnom,
	double &_rdenom,
	double &_snom)
{
	_rnom = 0;
	_rdenom = 1;
	_snom = 0;

	//if start==end, then use pt distance
	if (A.x == B.x && A.y == B.y)
		return distancesq_pt_pt(A, B);

	//otherwise, we use comp.graphics.algorithms Frequently Asked Questions method

	/*(1)     	      AC dot AB
				r = ---------
						||AB||^2
		r has the following meaning:
		r=0 P = A
		r=1 P = B
		r<0 P is on the backward extension of AB
		r>1 P is on the forward extension of AB
		0<r<1 P is interior to AB
	*/

	const double rdenom = distancesq_pt_pt(A, B);
	const double pdx = p.x - A.x, dx = B.x - A.x;
	const double pdy = p.y - A.y, dy = B.y - A.y;
	const double rnom1 = pdx * dx;
	const double rnom2 = pdy * dy;
	const double rnom = rnom1 + rnom2;
	_rdenom = rdenom;

	const double snom1 = pdx * dy;
	const double snom2 = pdy * dx;
	const double snom = snom1 - snom2;
	_snom = snom;

	if (rnom < -1)
	{
		_rnom = 0;
		return distancesq_pt_pt(p, A);
	}
	if (rnom > rdenom)
	{
		_rnom = rdenom;
		return distancesq_pt_pt(p, B);
	}
	_rnom = rnom;

	/*(2)
			(Ay-Cy)(Bx-Ax)-(Ax-Cx)(By-Ay)
		s = -----------------------------
						L^2

		Then the distance from C to P = |s|*L.

	*/

	return snom * snom / rdenom;
}

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


  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
               &map_waypoints_dx,&map_waypoints_dy]
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
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

		  //vector<double> calc_xy = getXY(car_s, car_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
		  //printf("%.2f, %.2f calc %.2f %.2f\n", car_x, car_y, calc_xy[0], calc_xy[1]);
		  //printf("XY calc diff %.2f, %.2f\n", car_x-calc_xy[0], car_y-calc_xy[1]);

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
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

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          /**
           * TODO: define a path made up of (x,y) points that the car will visit
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
			  pos_x = car_x;
			  pos_y = car_y;
			  angle = deg2rad(car_yaw);
		  }
		  else {
			  pos_x = previous_path_x[path_size - 1];
			  pos_y = previous_path_y[path_size - 1];

			  double pos_x2 = previous_path_x[path_size - 2];
			  double pos_y2 = previous_path_y[path_size - 2];
			  angle = atan2(pos_y - pos_y2, pos_x - pos_x2);
		  }
		  trajectory.add_point(pos_x, pos_y);
		  
		  int nNextWaypoint = NextWaypoint(pos_x, pos_y, angle, map_waypoints_x, map_waypoints_y);
		  int nPrevWaypoint = nNextWaypoint == 0 ? map_waypoints_x.size() - 1 : nNextWaypoint - 1;
		  double rnom, rdenom, snom;
		  double distsq = distancesq_pt_seg(
			  Point(pos_x, pos_y),
			  Point(map_waypoints_x[nNextWaypoint], map_waypoints_y[nNextWaypoint]),
			  Point(map_waypoints_x[nPrevWaypoint], map_waypoints_y[nPrevWaypoint]),
			  rnom, rdenom, snom);
		  double wp_len = distance(map_waypoints_x[nPrevWaypoint], map_waypoints_y[nPrevWaypoint],
			  map_waypoints_x[nNextWaypoint], map_waypoints_y[nNextWaypoint]);
		  double wp_dist = rnom < 0 ? 0 : rnom >= rdenom ? wp_len : wp_len * rnom / rdenom;
		  printf("next waypoint dist: s = %.2f d = %.2f\n", wp_dist, sqrt(distsq));
		  if (wp_dist < 2) // fix for very close waypoint to pos
		  {
			  nNextWaypoint = (nNextWaypoint + 1) % map_waypoints_x.size();
		  }
		  
		  for (int i = 0; i < 50 - path_size; ++i) {
			  int nPrevWaypoint = nNextWaypoint == 0 ? map_waypoints_x.size()-1: nNextWaypoint - 1;
			  double heading = atan2((map_waypoints_y[nNextWaypoint] - map_waypoints_y[nPrevWaypoint]),
				  (map_waypoints_x[nNextWaypoint] - map_waypoints_x[nPrevWaypoint]));
			  double perp_heading = heading - pi() / 2;
			  double d = 10;
			  double dx = d * cos(perp_heading);
			  double dy = d * sin(perp_heading);
			  trajectory.add_point(map_waypoints_x[nNextWaypoint]+dx, map_waypoints_y[nNextWaypoint]+dy);
			  if (trajectory.total_dist > 50) break;
			  nNextWaypoint= (nNextWaypoint+1)%map_waypoints_x.size();
		  }
		  double dist_total = 0;
		  double dist_step = 0.9;
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

		  /* Wiggly 
		  static int counter = 0;
		  counter++;
		  double dist_inc = 0.5;
		  double angle_inc = (pi() / 100)*sin(counter/100.0);
		  for (int i = 0; i < 50 - path_size; ++i) {
			  next_x_vals.push_back(pos_x + (dist_inc)*cos(angle + (i + 1)*angle_inc));
			  next_y_vals.push_back(pos_y + (dist_inc)*sin(angle + (i + 1)*angle_inc));
			  pos_x += (dist_inc)*cos(angle + (i + 1)*angle_inc);
			  pos_y += (dist_inc)*sin(angle + (i + 1)*angle_inc);
		  }*/

		  


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