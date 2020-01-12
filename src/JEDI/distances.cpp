#include <iostream>
#include "distances.h"

distances::distances()
{
 
}

/*
Once the grid is in place, the first step is to calculate the matrices that contain the euclidean distances
between atoms j and grid points i and their derivatives with respect to the coordinates x y and z of each atom.
These matrices have the form of vector<vector<double>> in which the positions of the external vector represent
grid points i and the positions of the internal vectors represent atoms j. 
This is expensive in terms of memory, but since these distances and their derivatives are calculated a high
number of times during the calculation of the JEDI derivatives, storing them in memory makes the code 
significantly more efficient. 
Regarding the calculation of the derivatives, it must be noted that the numerator
is the coordinate of the atom minus the coordinate of the grid. This is important because the derivatives of JEDI
are calculated with respect to the atom coordinates. This order guarantees the correct directionality.
*/
void distances::compute_distance_matrix(vector<PLMD::Vector> protein, vector<PLMD::Vector> grid)
{
  vector<double> r_ij_vec(protein.size(),0); //inside vector
  r_matrix=vector<vector<double>>(grid.size(),r_ij_vec);
  dr_matrix_dx=vector<vector<double>>(grid.size(),r_ij_vec);
  dr_matrix_dy=vector<vector<double>>(grid.size(),r_ij_vec);
  dr_matrix_dz=vector<vector<double>>(grid.size(),r_ij_vec);

  #pragma omp parallel for
  for (unsigned i=0; i<grid.size();i++)
  {
      for (unsigned j=0; j<protein.size(); j++)
      {
          r_matrix[i][j]=delta(protein[j],grid[i]).modulo();
          if (r_matrix[i][j]<0.0000001)
          {
            r_matrix[i][j]=0.0000001;
            dr_matrix_dx[i][j]=0;
            dr_matrix_dy[i][j]=0;
            dr_matrix_dz[i][j]=0;
          }
          else
          {
           dr_matrix_dx[i][j]=(protein[j][0]-grid[i][0])/r_matrix[i][j];
           dr_matrix_dy[i][j]=(protein[j][1]-grid[i][1])/r_matrix[i][j];
           dr_matrix_dz[i][j]=(protein[j][2]-grid[i][2])/r_matrix[i][j]; 
          }
      }
  }
  
  //Uncomment the following lines for testing
  /*
  for (unsigned i=0; i<grid.size();i++)
  {
    for (unsigned j=0; j<protein.size(); j++)
    {
      cout << " Distance atom " << j << " - gridpoint " << i << " (derivatives) = " << r_matrix[i][j] << "("
                                                                                    << dr_matrix_dx[i][j] << " "
                                                                                    << dr_matrix_dy[i][j] << " "
                                                                                    << dr_matrix_dz[i][j] << " "
                                                                                    << ")" << endl;
    }
  }
  */
}

