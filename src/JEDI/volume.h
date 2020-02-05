#include "activity.h"

using namespace std;

#ifndef volume_h_
#define volume_h_

class Volume
{
 public:
  Volume();
  vector<double> hydroactivity_grid;
  vector<vector<double>> d_hydroactivity_dx;
  vector<vector<double>> d_hydroactivity_dy;
  vector<vector<double>> d_hydroactivity_dz;
  
  double volume;
  vector<double> d_volume_dx;
  vector<double> d_volume_dy;
  vector<double> d_volume_dz;
  
  void compute_hydroactivity(vector<double> activity_grid,
                             vector<vector<double>> d_activity_dx,vector<vector<double>> d_activity_dy,vector<vector<double>> d_activity_dz,
                             vector<double> hydrophobicity_grid,
                             vector<vector<double>> d_hydrophobicity_dx,vector<vector<double>> d_hydrophobicity_dy,vector<vector<double>> d_hydrophobicity_dz);
  
  void compute_volume(vector<double> hydroactivity_grid,
                      vector<vector<double>> d_hydroactivity_dx,vector<vector<double>> d_hydroactivity_dy,vector<vector<double>> d_hydroactivity_dz,
                      double resolution);
};
#endif