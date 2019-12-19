#include "distances.h"

distances::distances()
{

}

/*
Once the grid is in place, the first step is to calculate the matrices that contain the euclidean distances
between atoms j and grid points i and their derivatives with respect to the coordinates x y and z of each atom.
These matrices have the form of vector<vector<double>> in which the positions of the external vector represent
atoms j and the positions of the internal vectors represent grid points i. 
This is expensive in terms of memory, but since these distances and their derivatives are calculated a high
number of times during the calculation of the JEDI derivatives, storing them in memory makes the code 
significantly more efficient. 
Regarding the calculation of the derivatives, it must be noted that the numerator
is the coordinate of the atom minus the coordinate of the grid. This is important because the derivatives of JEDI
are calculated with respect to the atom coordinates. This order guarantees the correct directionality.
*/
void distances::compute_distance_matrix(vector<PLMD::Vector> &protein, vector<PLMD::Vector> &grid)
{
  vector<double> r_ji_vec(grid.size(),0); //inside vector
  
  //ugly way of initialising the vectors
  for (unsigned j=0; j<protein.size();j++)
  {
   r_matrix.push_back(r_ji_vec);
   dr_matrix_dx.push_back(r_ji_vec);
   dr_matrix_dy.push_back(r_ji_vec);
   dr_matrix_dz.push_back(r_ji_vec);
  }

  #pragma omp parallel for
  for (unsigned j=0; j<protein.size();j++)
  {
      for (unsigned i=0; i<grid.size(); i++)
      {
          r_matrix[j][i]=sqrt(pow((protein[j][0]-grid[i][0]),2)+pow((protein[j][1]-grid[i][1]),2)+pow((protein[j][2]-grid[i][2]),2));
          dr_matrix_dx[j][i]=(protein[j][0]-grid[i][0])/r_matrix[j][i];
          dr_matrix_dy[j][i]=(protein[j][1]-grid[i][1])/r_matrix[j][i];
          dr_matrix_dz[j][i]=(protein[j][2]-grid[i][2])/r_matrix[j][i];
      }
  }
}

