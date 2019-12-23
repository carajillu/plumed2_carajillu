#include <iostream>
#include "kernel.h"

S_on_grid::S_on_grid()
{

}

void S_on_grid::compute_S_on_grid(vector<double> &mindist, 
                           vector<vector<double>> &d_mindist_dx,
                           vector<vector<double>> &d_mindist_dy,
                           vector<vector<double>> &d_mindist_dz,
                           double &v0, double &delta_v)
{
 //Ugly way to fill vectors with zeros
 vector<double> d_Son_grid_i(0,d_mindist_dx.size());
 for (unsigned i=0; i< mindist.size(); i++)
 {
    S_on_grid_vector.push_back(0);
    d_Son_grid_dx.push_back(d_Son_grid_i);
    d_Son_grid_dy.push_back(d_Son_grid_i);
    d_Son_grid_dz.push_back(d_Son_grid_i);
 }

 #pragma omp parallel for
 for (unsigned i=0; i<mindist.size();i++)
 {
   S_on S_on_grid_i(mindist[i],v0,delta_v,d_mindist_dx[i],d_mindist_dy[i],d_mindist_dz[i]);
   S_on_grid_i.compute_S_on();
   S_on_grid_vector[i]=S_on_grid_i.S_on_value;
   d_Son_grid_dx[i]=S_on_grid_i.d_Son_dx;
   d_Son_grid_dy[i]=S_on_grid_i.d_Son_dy;
   d_Son_grid_dz[i]=S_on_grid_i.d_Son_dz;
 }

 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<S_on_grid_vector.size();i++)
 {
     cout << "Mindist point " << i << ": " << mindist[i] << endl;
     cout << "S_on_grid point "<< i << ": " << S_on_grid_vector[i] << endl;
     for (unsigned j=0; j<d_Son_grid_dx[i].size();j++)
     {
        cout << "Iteration " << j << ": "<< "derivatives with respect to atom " << j << ": " << d_Son_grid_dx[i][j] << " " <<
                                                                    d_Son_grid_dy[i][j] << " " <<
                                                                    d_Son_grid_dz[i][j] << endl;
     }
 }
 */
 
}