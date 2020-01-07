#include <iostream>
#include "distances.h"
#include "kernel.h"

using namespace std;

contacts::contacts()
{

}
/*
This will take the distance matrix and its derivatives with respect to every atom
and return the S_off values of the distance for each gridpoint-atom pair.
*/
void contacts::compute_contacts(vector<vector<double>> r_matrix, 
                                 vector<vector<double>> dr_matrix_dx, 
                                 vector<vector<double>> dr_matrix_dy, 
                                 vector<vector<double>> dr_matrix_dz,
                                 double CC2min, double deltaCC2)
{
  unsigned size_grid=r_matrix.size();
  unsigned size_protein=r_matrix[0].size();
  
  vector<double> contacts_i(size_protein,0);
  contacts_matrix=vector<vector<double>>(size_grid,contacts_i);
  d_contacts_dx=vector<vector<double>>(size_grid,contacts_i);
  d_contacts_dy=vector<vector<double>>(size_grid,contacts_i);
  d_contacts_dz=vector<vector<double>>(size_grid,contacts_i);

  #pragma omp parallel for
  for (unsigned i=0; i<size_grid;i++)
  {
   for (unsigned j=0; j<size_protein;j++)
   { 
     vector<double> d_contactsij_dx(1,dr_matrix_dx[i][j]);
     vector<double> d_contactsij_dy(1,dr_matrix_dy[i][j]);
     vector<double> d_contactsij_dz(1,dr_matrix_dz[i][j]);
     S_off contacts_ij(r_matrix[i][j],CC2min,deltaCC2,d_contactsij_dx,d_contactsij_dy,d_contactsij_dz);
     contacts_ij.compute_S_off();
     contacts_matrix[i][j]=contacts_ij.S_off_value;
     d_contacts_dx[i][j]=contacts_ij.d_Soff_dx[0];
     d_contacts_dy[i][j]=contacts_ij.d_Soff_dy[0];
     d_contacts_dz[i][j]=contacts_ij.d_Soff_dz[0];
   }
 }

  //Uncomment the following lines for testing
  
  for (unsigned i=0; i<size_grid;i++)
  {
    for (unsigned j=0; j<size_protein; j++)
    {
      cout << " S_off(r,cc2min,deltacc2) atom " << j << " - gridpoint " << i << " (derivatives) = " << contacts_matrix[i][j] << "("
                                                                                                    << d_contacts_dx[i][j] << " "
                                                                                                    << d_contacts_dy[i][j] << " "
                                                                                                    << d_contacts_dz[i][j] << ")"
                                                                                                    << ")" << endl;
    }
  }
  

}

