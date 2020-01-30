#include<vector>
#include "tools/PDB.h"
using namespace std;

#ifndef clustering_h_
#define clustering_h_
class clustering
{
 private:
  struct laio
  {
   double density;
   unsigned idx;
   double activity;
   int nnhd;
   double delta;
   int cluster;
  };
  static bool sort_grid(const laio &a, const laio &b);
  
 public:
  clustering();
  vector<laio> grid_stats;
  vector<vector<unsigned>> clusters;
  void cluster_grid(vector<double> activity, vector<vector<double>> r_matrix,vector<vector<unsigned>> neighbours,double GP_max, double grid_resolution);
  void print_clusters(vector<PLMD::Vector> grid);
};
#endif
