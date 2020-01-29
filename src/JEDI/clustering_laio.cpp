#include "clustering_laio.h"
#include <algorithm>
#include <iostream>
using namespace std;

clustering::clustering()
{

}

bool clustering::sort_density(const vector5d a, const vector5d b)
{
  return a.density>b.density;
}

void clustering::cluster_grid(vector<double> activity, vector<vector<double>> r_matrix,vector<vector<unsigned>> neighbours, double GP_max)
{
  unsigned size_grid = activity.size();
  vector<double> density(size_grid,0);
  vector<double> delta(size_grid,9999999); //not pretty! any way to use infinite?
  vector<int> nnhd(size_grid,-1);
  vector<int> cluster_assign(size_grid,-1);
  vector<unsigned> clusters_l;
  vector<pair<double,unsigned>> density_pair(size_grid,make_pair(0,0));

  #pragma omp parallel
  {
    #pragma omp for
    for (unsigned i=0; i<size_grid;i++)
    {
      double rho=0;
      for (unsigned k=0; k<neighbours[i].size();k++)
      {
        unsigned neighbour_idx=neighbours[i][k];
        if (activity[neighbour_idx]==0) continue;
        rho++;
      }
      density[i]=rho;
      density_pair[i]=make_pair(rho,i);
    }
    #pragma omp barrier
    #pragma omp single
    {
     sort(density_pair.rbegin(),density_pair.rend());
     for (unsigned i=0; i<size_grid; i++)
     {
       cout << density_pair[i].first << " " << density_pair[i].second << endl;
     }
    }
    #pragma omp barrier

    #pragma omp single
    for(unsigned i=0; i<size_grid;i++)
    {
      for (unsigned k=0; k<size_grid; k++)
      {
        if (i==k) continue;
        if ((density[k]>density[i]) and (r_matrix[i][k]<delta[i]))
        {
          delta[i]=r_matrix[i][k];
          nnhd[i]=k;
        }
      }
    }
    #pragma omp barrier

    #pragma omp single
    {
      for (unsigned i=0; i<size_grid;i++)
      {
       unsigned grid_idx=density_pair[i].second;
       cout << "processing point " << grid_idx << " with delta = " << delta[grid_idx] << endl;
       if (delta[grid_idx]>GP_max) //If there is no neighbour with a higher density, create a new cluster
       {
        clusters.push_back(clusters_l);
        clusters[clusters.size()-1].push_back(grid_idx);
        cluster_assign[grid_idx]=clusters.size()-1;
       }
       else
       {
        cluster_assign[grid_idx]=cluster_assign[nnhd[grid_idx]];
        clusters[cluster_assign[grid_idx]].push_back(grid_idx);
       }
      }
    }

  }
 
  int count=0;
  for (unsigned i=0; i<clusters.size();i++)
  {
    cout << "Cluster " << i << " Has " << clusters[i].size() << " elements" << endl;
    count+=clusters[i].size();
  }
  cout << "Total points clustered: " << count << endl;
  cout << "Total number of grid points " << size_grid << endl;
  exit(0);
}