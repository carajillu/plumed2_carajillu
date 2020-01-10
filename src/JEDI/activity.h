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
                           vector<double> apolar_contacts,
                           vector<vector<double>> d_apolar_dx,
                           vector<vector<double>> d_apolar_dy,
                           vector<vector<double>> d_apolar_dz,
                           vector<double> polar_contacts,
                           vector<vector<double>> d_polar_dx,
                           vector<vector<double>> d_polar_dy,
                           vector<vector<double>> d_polar_dz,
                           double Emin, double deltaE);

};

#endif