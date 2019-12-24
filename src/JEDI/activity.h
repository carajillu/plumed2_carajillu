#include<vector>
#include "tools/PDB.h"
#include "kernel.h"

using namespace std;

#ifndef activity_h_
#define activity_h_

class Activity
{
 private:
   vector<vector<int>> neighbours; //move to private after debug
   void compute_neighbours(vector<PLMD::Vector> &positions,double &GPmin, double &GPmax);

   vector<double> depth_sum;
   vector<vector<double>> d_depthsum_dx;
   vector<vector<double>> d_depthsum_dy;
   vector<vector<double>> d_depthsum_dz;
   void compute_farawayness_sum(vector<vector<int>> &neighbours,S_off_grid farawayness);

   vector<double> depth;
   vector<vector<double>> d_depth_dx;
   vector<vector<double>> d_depth_dy;
   vector<vector<double>> d_depth_dz;
   void compute_depth();

 public:
   Activity();
   void Activity_init(vector<PLMD::Vector> &positions,
            double &CC_min, double deltaCC,
            double &GPmin, double &GPmax,
            double &CC2_min, double &deltaCC2);
   vector<double> activity;
   vector<vector<double>> d_activity_dx;
   vector<vector<double>> d_activity_dy;
   vector<vector<double>> d_activity_dz;
   void compute_activities(vector<double> &mindist, 
                           vector<vector<double>> &d_mindist_dx,
                           vector<vector<double>> &d_mindist_dy,
                           vector<vector<double>> &d_mindist_dz,
                           double &CC_min, double &deltaCC,
                           double &GPmin, double &GPmax,
                           double &CC2_min, double &deltaCC2);
   
   
   
   
};

#endif