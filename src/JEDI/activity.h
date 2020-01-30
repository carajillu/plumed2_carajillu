#include<vector>
#include "tools/PDB.h"
#include "kernel.h"

using namespace std;

#ifndef activity_h_
#define activity_h_

class activity
{
 private:
   vector<double> S_on_mindist;
   vector<double> S_on_contacts;
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
   vector<unsigned> active_indices;
   void compute_activities(vector<double> mindist, 
                           vector<vector<double>> d_mindist_dx,
                           vector<vector<double>> d_mindist_dy,
                           vector<vector<double>> d_mindist_dz,
                           double CCmin, double deltaCC,
                           vector<double> total_contacts,
                           vector<vector<double>> d_contacts_total_dx,
                           vector<vector<double>> d_contacts_total_dy,
                           vector<vector<double>> d_contacts_total_dz,
                           double Emin, double deltaE);
  void filter_activities();
  void print_activities();

};

#endif