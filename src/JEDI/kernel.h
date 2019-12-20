#include <vector>
#include "tools/PDB.h"

using namespace std;

#ifndef jedi_kernel_h_
#define jedi_kernel_h_

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