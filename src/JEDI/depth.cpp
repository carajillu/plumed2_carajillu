#include <iostream>
#include "activity.h"
#include "kernel.h"

void Activity::compute_farawayness_sum(vector<vector<int>> &neighbours, S_off_grid farawayness)
{
  vector<double> d_depthsum_i(farawayness.d_Soff_grid_dx[0].size(),0);
  for (unsigned i=0; i<neighbours.size();i++)
  {
    depth_sum.push_back(0);
    d_depthsum_dx.push_back(d_depthsum_i);
    d_depthsum_dy.push_back(d_depthsum_i);
    d_depthsum_dz.push_back(d_depthsum_i);
  }

  #pragma omp parallel for
  for (unsigned i=0; i<neighbours.size();i++)
  {
      for (unsigned k=0; k<neighbours[i].size();k++)
      {
          cout << "point " << neighbours[i][k] << " is a neighbour of point " << i << endl;
          int neighbour_idx = neighbours[i][k];
          depth_sum[i]+=farawayness.S_off_grid_vector[neighbour_idx];
          for (unsigned j=0; j<d_depthsum_dx[i].size();j++)
          {  
            d_depthsum_dx[i][j]+=farawayness.d_Soff_grid_dx[neighbour_idx][j];
            d_depthsum_dy[i][j]+=farawayness.d_Soff_grid_dy[neighbour_idx][j];
            d_depthsum_dz[i][j]+=farawayness.d_Soff_grid_dz[neighbour_idx][j];
          }
      }
  }

  //Uncomment the following lines for testing
  cout << "Testing neighbours farawayness" << endl;
  for (unsigned i=0; i<neighbours.size();i++)
  {
   cout << "farawayness point " << i << " = " << farawayness.S_off_grid_vector[i] << endl;
   cout << "depth_sum point " << i << " = " << depth_sum[i] << endl;
   for (unsigned j=0; j<d_depthsum_dx[i].size();j++)
   { 
     cout << "derivatives with respect to atom " << j << ": " << d_depthsum_dx[i][j] << " " 
                                                              << d_depthsum_dy[i][j] << " " 
                                                              << d_depthsum_dz[i][j] << " " << endl;
   }
   
  }

}