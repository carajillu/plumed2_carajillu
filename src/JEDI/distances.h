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
};
#endif