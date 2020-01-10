#include "tools/PDB.h"
#include "distances.h"
#include "kernel.h"
#include "activity.h"

using namespace std;

#ifndef hydrophobicity_h_
#define hydrophobicity_h_

class hydrophobicity
{
 public:
   hydrophobicity();
   double Ha;
   vector<double> d_Ha_dx;
   vector<double> d_Ha_dy;
   vector<double> d_Ha_dz;
   void compute_hydrophobicity(vector<double> contacts_apolar,
                               vector<vector<double>> d_apolar_dx,vector<vector<double>> d_apolar_dy,vector<vector<double>> d_apolar_dz,
                               vector<double> contacts_polar,
                               vector<vector<double>> d_polar_dx,vector<vector<double>> d_polar_dy,vector<vector<double>> d_polar_dz,
                               vector<double> activity,
                               vector<vector<double>> d_activity_dx,vector<vector<double>> d_activity_dy,vector<vector<double>> d_activity_dz,
                               double sum_activity,
                               vector<double> d_sum_activity_dx, vector<double> d_sum_activity_dy,vector<double> d_sum_activity_dz);
   
   

};
#endif