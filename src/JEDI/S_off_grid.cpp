#include <iostream>
#include "kernel.h"

S_off_grid::S_off_grid()
{

}

void S_off_grid::compute_S_off_grid(vector<double> &mindist, 
                           vector<vector<double>> &d_mindist_dx,
                           vector<vector<double>> &d_mindist_dy,
                           vector<vector<double>> &d_mindist_dz,
                           double &v0, double &delta_v)
{
 //Empty vector in case they have something
 //Should not be necessary, but let's stay safe
 S_off_grid_vector.empty();
 d_Soff_grid_dx.empty();
 d_Soff_grid_dy.empty();
 d_Soff_grid_dz.empty();

 //Ugly way to fill vectors with zeros
 for (unsigned i=0; i< mindist.size(); i++)
 {
    S_off_grid_vector.push_back(0);
 }

 vector<double> d_Soff_grid_i(0,d_mindist_dx.size());
 for (unsigned i=0; i<mindist.size();i++)
 {
  d_Soff_grid_dx.push_back(d_Soff_grid_i);
  d_Soff_grid_dy.push_back(d_Soff_grid_i);
  d_Soff_grid_dz.push_back(d_Soff_grid_i);
 }

 #pragma omp parallel for
 for (unsigned i=0; i<mindist.size();i++)
 {
   S_off S_off_grid_i(mindist[i],v0,delta_v,d_mindist_dx[i],d_mindist_dy[i],d_mindist_dz[i]);
   S_off_grid_i.compute_S_off();
   S_off_grid_vector[i]=S_off_grid_i.S_off_value;
   d_Soff_grid_dx[i]=S_off_grid_i.d_Soff_dx;
   d_Soff_grid_dy[i]=S_off_grid_i.d_Soff_dy;
   d_Soff_grid_dz[i]=S_off_grid_i.d_Soff_dz;
 }

 // Uncomment the following lines for testing
 
 for (unsigned i=0; i<S_off_grid_vector.size();i++)
 {
     cout << "Mindist point " << i << ": " << mindist[i] << endl;
     cout << "S_off_grid point "<< i << ": " << S_off_grid_vector[i] << endl;
     for (unsigned j=0; j<d_Soff_grid_dx[i].size();j++)
     {
        cout << "Iteration " << j << ": "<< "derivatives with respect to atom " << j << ": " << d_Soff_grid_dx[i][j] << " " <<
                                                                    d_Soff_grid_dy[i][j] << " " <<
                                                                    d_Soff_grid_dz[i][j] << endl;
     }
 }
 
}