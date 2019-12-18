#include "distances.h"

distances::distances()
{

}

/*
This function returns a vector<vector<double>> in which the rows (external vector)
are the grid points (i) and the columns (internal vector) are the atoms (j). May be
memory expensive but the derivatives calculate that over and over so better to have it
in the memory.
*/
void distances::compute_distance_matrix(vector<PLMD::Vector> &protein, vector<PLMD::Vector> &grid)
{
  vector<double> r_ij_vec(protein.size(),0); //inside vector
  
  for (unsigned i=0; i<grid.size();i++)
  {
      r_matrix.push_back(r_ij_vec);
      for (unsigned j=0; j<protein.size(); j++)
      {
          r_matrix[i][j]=sqrt(pow((grid[i][0]-protein[j][0]),2)+pow((grid[i][0]-protein[j][2]),2)+pow((grid[i][2]-protein[j][2]),2));
      }
  }
}