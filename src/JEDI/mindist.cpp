#include <iostream>
#include "distances.h"

/*
This function returns the following vectors:
1) vector<double> mindist(grid_size): contains the mindist according to Eq.5 JCTC 2015 of every grid point.
                                      This funtion is necessary because we need a minimum distance
                                      parameter that ensures that JEDI is continuous and differentiable
                                      with respect to the coordinates of all atoms and not only the 
                                      one closest to a given grid point.
2) vector<vector<double>> d_mindist_dx(protein_size,grid_size): for every protein atom (outside vector),
                          contains the derivatives of each grid point i's (mindist) with respect of the x
                          coordinate of each atom j.
3) vector<vector<double>> d_mindist_dy(protein_size,grid_size): same as 2) but for the y coordinate
3) vector<vector<double>> d_mindist_dz(protein_size,grid_size): same as 2) but for the z coordinate
*/
void mindist::compute_mindist(vector<vector<double>> &r_matrix, 
                              vector<vector<double>> &dr_matrix_dx, 
                              vector<vector<double>> &dr_matrix_dy, 
                              vector<vector<double>> &dr_matrix_dz,
                              double &theta)
{
  unsigned grid_size=r_matrix.size();
  unsigned protein_size=r_matrix[0].size();

  vector<vector<double>> exp_r; //exponential term of Eq 5 JCTC 2015
  
  vector<double> protein_zeros(protein_size,0);
  vector<double> grid_zeros(grid_size,0);
  
  // Fill vectors with zeros
  min_dist=grid_zeros;
  for (unsigned j=0; j<grid_size;j++)
  {
      d_mindist_dx.push_back(protein_zeros);
      d_mindist_dy.push_back(protein_zeros);
      d_mindist_dz.push_back(protein_zeros);
      exp_r.push_back(protein_zeros);
  }
  
  //Calculate mindist and its derivatives
  #pragma omp parallel for
  for (unsigned i=0; i<grid_size;i++)
  {
    //calculate mindist
    double sum_exp_r=0;
    for (unsigned j=0; j<protein_size;j++)
    {
      exp_r[i][j]= exp(theta/r_matrix[i][j]);
      sum_exp_r+=exp_r[i][j];
    }
    min_dist[i]=theta/log(sum_exp_r);
    //calculate mindist derivatives
    for (unsigned j=0; j<protein_size;j++)
    {
     d_mindist_dx[i][j]=(pow(min_dist[i],2)*exp_r[i][j])/(sum_exp_r*pow(r_matrix[i][j],2))*dr_matrix_dx[i][j];
     d_mindist_dy[i][j]=(pow(min_dist[i],2)*exp_r[i][j])/(sum_exp_r*pow(r_matrix[i][j],2))*dr_matrix_dy[i][j];
     d_mindist_dz[i][j]=(pow(min_dist[i],2)*exp_r[i][j])/(sum_exp_r*pow(r_matrix[i][j],2))*dr_matrix_dz[i][j];
    }
  }

//Uncomment the following lines for testing

for (unsigned i=0; i< min_dist.size(); i++)
{
  cout << "Mindist grid point " <<  i << ": " << min_dist[i] << endl; 
  for (unsigned j=0; j<d_mindist_dx[i].size();j++)
  {
   cout << "Derivatives with respect to atom " << j << ": " << d_mindist_dx[i][j]<< " " 
                                                            << d_mindist_dy[i][j]<< " " 
                                                            << d_mindist_dz[i][j] << endl;
  }
}

}