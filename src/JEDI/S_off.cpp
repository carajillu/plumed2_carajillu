#include <iostream>
#include "kernel.h"

/*
The S_off object is initialised and the values of all its derivatives are set to 0.
This helps speed up the code by simply not calculating them when m equals 0 or 1.
*/
S_off::S_off()
{
  pi=3.14159265359;
}

/*
This function computes the value of S_off as a function of m_value. 
If m is different than 1 or 0, the derivatives of S_off with respect to the coordinates
of the atoms in the binding site are calculated as well.
*/
void S_off::compute_S_off(double v, double v0, double delta_v, vector<double> dv_dx, vector<double> dv_dy, vector<double> dv_dz)
{
  d_Soff_dx=vector<double>(dv_dx.size(),0);
  d_Soff_dy=vector<double>(dv_dy.size(),0);
  d_Soff_dz=vector<double>(dv_dz.size(),0);

  if (v<v0)
  {
    S_off_value = 1;
  }
  else
  {
    //S_off_value = (1/(delta_v*sqrt(2*pi)))*exp((-pow((v-v0),2)/(2*pow(delta_v,2))));
    S_off_value = exp((-pow((v-v0),2)/(2*pow(delta_v,2))));
    double d_Soff_dv=-S_off_value*((v-v0)/pow(delta_v,2)); // First part of the derivative
    #pragma omp parallel for
    for (unsigned j=0; j<dv_dx.size();j++)
    {
        d_Soff_dx[j]=d_Soff_dv*dv_dx[j];
        d_Soff_dy[j]=d_Soff_dv*dv_dy[j];
        d_Soff_dz[j]=d_Soff_dv*dv_dz[j];
    } 
  }
  // Uncomment the following lines for testing
  /*
  cout << "S_on = " << S_off_value << endl;
  for (unsigned j=0; j<d_Soff_dx.size();j++)
  {
   cout << "Derivatives with respect to atom " << j << ": " << d_Soff_dx[j] << " " 
                                                            << d_Soff_dy[j] << " " 
                                                            << d_Soff_dz[j] << endl;
  }
  */
}