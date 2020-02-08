#include <vector>
#include "tools/PDB.h"

using namespace std;

#ifndef jedi_distances_h_
#define jedi_distances_h_

class distances
{
 public:
   distances();
   vector<vector<double>> r_matrix;
   vector<vector<double>> dr_matrix_dx;
   vector<vector<double>> dr_matrix_dy;
   vector<vector<double>> dr_matrix_dz;
   void compute_distance_matrix(vector<PLMD::Vector> &protein, vector<PLMD::Vector> &grid);
   void filter_distance_matrix(vector<unsigned> &indices);
};

class mindist
{
 public:
  mindist();
  vector<double> min_dist;
  vector<vector<double>> d_mindist_dx;
  vector<vector<double>> d_mindist_dy;
  vector<vector<double>> d_mindist_dz;
  void compute_mindist(vector<vector<double>> &r_matrix, 
                       vector<vector<double>> &dr_matrix_dx, 
                       vector<vector<double>> &dr_matrix_dy, 
                       vector<vector<double>> &dr_matrix_dz,
                       double &theta);
};

#endif