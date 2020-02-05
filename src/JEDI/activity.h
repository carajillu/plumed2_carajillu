#include<vector>
#include "tools/PDB.h"
#include "kernel.h"

using namespace std;

#ifndef activity_h_
#define activity_h_

class activity
{
 public:
   activity();
   vector<double> activity_grid;
   vector<double> S_on_mindist;
   vector<double> S_off_mindist;
   vector<vector<double>> d_activity_dx;
   vector<vector<double>> d_activity_dy;
   vector<vector<double>> d_activity_dz;
   double sum_activity;
   vector<double> d_sum_activity_dx;
   vector<double> d_sum_activity_dy;
   vector<double> d_sum_activity_dz;
   void compute_activities(vector<double> mindist, 
                           vector<vector<double>> d_mindist_dx,
                           vector<vector<double>> d_mindist_dy,
                           vector<vector<double>> d_mindist_dz,
                           double CCmin, double deltaCC,
                           double r_hydro, double deltar_hydro);
  void filter_activities(vector<unsigned> cluster);
  void print_activities();

};

#endif