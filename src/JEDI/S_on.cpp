#include <iostream>
#include "kernel.h"

/*
The S_on object is initialised and the values of all its derivatives are set to 0.
This helps speed up the code by simply not calculating them when m equals 0 or 1.
*/
S_on::S_on(double v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz)
{
 v;
 k=1; //hardcoded because it's always the same
 vector<double> derivatives(dv_dx.size(),0);
 d_Son_dx=derivatives;
 d_Son_dy=derivatives;
 d_Son_dz=derivatives;
 dm_dx=derivatives;
 dm_dy=derivatives;
 dm_dz=derivatives;
 compute_m(v,v0,delta_v,dv_dx,dv_dy,dv_dz);
}

/* 
This function computes m_value and its derivatives, which are then passed on to compute_S_on.
Since m_value is not used anywhere else in the code, this function is private to this class.
*/
void S_on::compute_m(double &v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz)
{
    m_value=(v-v0)/delta_v;

    if (m_value <= 0 or 1 <= m_value) return; // don't calculate derivatives if S_on is flat


    #pragma omp parallel for
    for (unsigned j=0; j<dv_dx.size();j++)
    {
        dm_dx[j]=1/delta_v*dv_dx[j];
        dm_dy[j]=1/delta_v*dv_dy[j];
        dm_dz[j]=1/delta_v*dv_dz[j];
    }
    // Uncomment the following lines for testing
    /*
    cout << "m_value = " << m_value << endl;
    for (unsigned j=0; j<dm_dx.size();j++)
    {
     cout << "Derivatives with respect to atom " << j << ": " << dm_dx[j] << " " 
                                                              << dm_dy[j] << " " 
                                                              << dm_dz[j] << endl;
    }
    */
}

/*
This function computes the value of S_on as a function of m_value. 
If m is different than 1 or 0, the derivatives of S_on with respect to the coordinates
of the atoms in the binding site are calculated as well.
*/
void S_on::compute_S_on()
{
  if (m_value<=0)
  {
    S_on_value = 0;
  }
  else if (0<m_value and m_value<1) 
  {
    S_on_value = k*(3*pow(m_value,4)-2*(pow(m_value,6)));
    
    double d_Son_dm = 12*k*(pow(m_value,3)-pow(m_value,5)); // First part of the derivative
    #pragma omp parallel for
    for (unsigned j=0; j<dm_dx.size();j++)
    {
        d_Son_dx[j]=d_Son_dm*dm_dx[j];
        d_Son_dy[j]=d_Son_dm*dm_dy[j];
        d_Son_dz[j]=d_Son_dm*dm_dz[j];
    }
  }
  else 
  {
    S_on_value = k;
  }


  // Uncomment the following lines for testing
  /*
  cout << "S_on = " << S_on_value << endl;
  for (unsigned j=0; j<d_Son_dx.size();j++)
  {
   cout << "Derivatives with respect to atom " << j << ": " << d_Son_dx[j] << " " 
                                                            << d_Son_dy[j] << " " 
                                                            << d_Son_dz[j] << endl;
  }
  */
}