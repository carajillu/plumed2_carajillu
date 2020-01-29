#include<vector>
#include "tools/PDB.h"
using namespace std;

#ifndef clustering_h_
#define clustering_h_
class clustering
{
 private:
  struct vector5d
  {
   unsigned idx;
   double density;
   int nnhd;
   double delta;
   int cluster;
  };
  bool sort_density(const vector5d a, const vector5d b);
 public:
  clustering();
  vector<vector<unsigned>> clusters;
  void cluster_grid(vector<double> activity, vector<vector<double>> r_matrix,vector<vector<unsigned>> neighbours);
};
#endif
