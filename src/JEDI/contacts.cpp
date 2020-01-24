#include <iostream>
#include "contacts.h"
#include "kernel.h"

using namespace std;

contacts_Soff::contacts_Soff()
{

}
/*
This will take the distance matrix and its derivatives with respect to every atom
and return the S_off values of the distance for each gridpoint-atom pair.
*/
void contacts_Soff::compute_contacts_S_off(vector<vector<double>> r_matrix, 
                                           vector<vector<double>> dr_matrix_dx, 
                                           vector<vector<double>> dr_matrix_dy, 
                                           vector<vector<double>> dr_matrix_dz,
                                           double v0, double delta_v)
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
     S_off contacts_ij;
     contacts_ij.compute_S_off(r_matrix[i][j],v0,delta_v,d_contactsij_dx,d_contactsij_dy,d_contactsij_dz);
     contacts_matrix[i][j]=contacts_ij.S_off_value;
     d_contacts_dx[i][j]=contacts_ij.d_Soff_dx[0];
     d_contacts_dy[i][j]=contacts_ij.d_Soff_dy[0];
     d_contacts_dz[i][j]=contacts_ij.d_Soff_dz[0];
   }
 }

  //Uncomment the following lines for testing
  /*
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
  */

}

contacts_sum::contacts_sum()
{

}

void contacts_sum::compute_contacts_sum(vector<vector<double>> contacts_matrix,
                                        vector<vector<double>> d_contacts_dx,
                                        vector<vector<double>> d_contacts_dy,
                                        vector<vector<double>> d_contacts_dz,
                                        vector<string> atomnames)
{
 unsigned size_grid=contacts_matrix.size();
 unsigned size_protein=contacts_matrix[0].size();
 vector<double> d_contacts(size_protein,0);

 contacts_apolar=vector<double>(size_grid,0);
 d_contacts_apolar_dx=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_apolar_dz=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_apolar_dy=vector<vector<double>>(size_grid,d_contacts);
 
 contacts_polar=vector<double>(size_grid,0);
 d_contacts_polar_dx=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_polar_dz=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_polar_dy=vector<vector<double>>(size_grid,d_contacts);

 contacts_total=vector<double>(size_grid,0);
 d_contacts_total_dx=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_total_dz=vector<vector<double>>(size_grid,d_contacts);
 d_contacts_total_dy=vector<vector<double>>(size_grid,d_contacts);
 
 //#pragma omp parallel for
 for (unsigned i=0; i<size_grid;i++)
 {
   for(unsigned j=0;j<size_protein;j++)
   {
    if (atomnames[j][0]=='C' or atomnames[j][0]=='S')
    {
     contacts_apolar[i]+=contacts_matrix[i][j];
     d_contacts_apolar_dx[i][j]=d_contacts_dx[i][j];
     d_contacts_apolar_dy[i][j]=d_contacts_dy[i][j];
     d_contacts_apolar_dz[i][j]=d_contacts_dz[i][j];
    }
    else if (atomnames[j][0]=='N' or atomnames[j][0]=='O')
    {
     contacts_polar[i]+=contacts_matrix[i][j];
     d_contacts_polar_dx[i][j]=d_contacts_dx[i][j];
     d_contacts_polar_dy[i][j]=d_contacts_dy[i][j];
     d_contacts_polar_dz[i][j]=d_contacts_dz[i][j];
    }
    else
    {
      cout << "Unknown atom detected at position: " << j << ": " << atomnames[j] << ". Exiting." << endl;
      exit(0);
    }
    d_contacts_total_dx[i][j]=d_contacts_apolar_dx[i][j]+d_contacts_polar_dx[i][j];
    d_contacts_total_dy[i][j]=d_contacts_apolar_dy[i][j]+d_contacts_polar_dy[i][j];
    d_contacts_total_dz[i][j]=d_contacts_apolar_dz[i][j]+d_contacts_polar_dz[i][j];
   }
   contacts_total[i]=contacts_apolar[i]+contacts_polar[i];
 }
 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<size_grid; i++)
 {
   cout << "Apolar contacts sum grid point " << i << " = " << contacts_apolar[i] << endl;
   for (unsigned j=0; j<size_protein; j++)
   {
     cout << "Derivatives with respect to atom " << j << " with name " << atomnames[j] << ": " 
     << d_contacts_apolar_dx[i][j] << " " << d_contacts_apolar_dy[i][j] << " " << d_contacts_apolar_dz[i][j] << endl;
   }

   cout << "Polar contacts sum grid point " << i << " = " << contacts_polar[i] << endl;
   for (unsigned j=0; j<size_protein; j++)
   {
     cout << "Derivatives with respect to atom " << j << " with name " << atomnames[j] << ": " 
     << d_contacts_polar_dx[i][j] << " " << d_contacts_polar_dy[i][j] << " " << d_contacts_polar_dz[i][j] << endl;
   }
 }
 */
}
