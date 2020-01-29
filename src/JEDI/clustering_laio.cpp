#include "clustering_laio.h"
#include <algorithm>
using namespace std;

clustering::clustering()
{

}

bool clustering::sort_density(const vector5d a, const vector5d b)
{
  return a.density>b.density;
}

void clustering::cluster_grid(vector<double> activity, vector<vector<double>> r_matrix,vector<vector<unsigned>> neighbours)
{
  unsigned size_grid = activity.size();
  vector<vector5d> grid_5d(grid_size);

  #pragma omp parallel for
  for (unsigned i=0; i<size_grid; i++)
  {
    if (activity[i]==0) continue;

    vector5d p;
    p.idx=i;
    p.density=0;
    p.nnhd=-1;
    p.delta=-1;
    p.cluster=-1;

    for (unsigned k=0; k<neighbours[i].size();k++)
    {
     p.density+=activity[neighbours[i][k]];
    }
    grid_5d[i]=p;
  }
   
  sort(grid_5d.begin(),grid_5d.end(),sort_density);
  
  unsigned cluster_id=0
  grid_5d[0].nnhd=cluster
  grid_5d[0].delta=99999999; // This is the point with the highest density
  grid_5d[0].cluster=0;
  for (unsigned i=1; i<size_grid; i++)
  {
   unsigned idx_i=grid_5d[i].idx;
   unsigned idx_k=grid_5d[i-1].idx;
   grid_5d[i].nnhd=idx_k;
   grid_5d[i].delta=r_matrix[idx_i][idx_k];
   if (grid_5d[i].delta>GPmax)
   {
    cluster_id++
    grid_5d[i].cluster=cluster_id
   }
   else
   {
    grid_5d[i].cluster=grid_5d[i-1].cluster;
   }
  }

  unsigned nclust=cluster_id+1;
  vector<unsigned> cluster_l;
  clusters=vector<vector<unsigned>>()
}