#include <iostream>
#include "activity.h"
#include "kernel.h"

void Activity::compute_depth(vector<vector<int>> neighbours, vector<S_off> proximity)
{
  vector<double> d_depth_i(proximity[0].d_Soff_dx.size());
  depth=vector<double>(neighbours.size(),0);
  d_depth_dx=vector<vector<double>>(neighbours.size(),d_depth_i);
  d_depth_dz=vector<vector<double>>(neighbours.size(),d_depth_i);
  d_depth_dy=vector<vector<double>>(neighbours.size(),d_depth_i);

  #pragma omp parallel for
  for (unsigned i=0; i<neighbours.size();i++)
  {
      for (unsigned k=0; k<neighbours[i].size();k++)
      {
          int neighbour_idx = neighbours[i][k];
          depth[i]+=proximity[neighbour_idx].S_off_value/max_neighbours;
          for (unsigned j=0; j<d_depth_dx[i].size();j++)
          {  
            d_depth_dx[i][j]+=proximity[neighbour_idx].d_Soff_dx[j]/max_neighbours;
            d_depth_dy[i][j]+=proximity[neighbour_idx].d_Soff_dy[j]/max_neighbours;
            d_depth_dz[i][j]+=proximity[neighbour_idx].d_Soff_dz[j]/max_neighbours;
          }
      }
  }

  //Uncomment the following lines for testing
  /*
  cout << "Testing neighbours proximity" << endl;
  for (unsigned i=0; i<neighbours.size();i++)
  {
   for (unsigned k=0;k<neighbours[i].size();k++)
   {
     cout << "proximity of point " << neighbours[i][k] <<"=" << proximity[neighbours[i][k]].S_off_value << endl;
   }
   cout << "depth point " << i << " = " << depth[i] << endl;
   for (unsigned j=0; j<d_depth_dx[i].size();j++)
   { 
     cout << "derivatives with respect to atom " << j << ": " << d_depth_dx[i][j] << " " 
                                                              << d_depth_dy[i][j] << " " 
                                                              << d_depth_dz[i][j] << " " << endl;
   }
  }
 */
}