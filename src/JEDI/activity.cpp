#include "activity.h"
#include "kernel.h"
#include <iostream>
#include <algorithm>
#include <functional>

using namespace std;

/*
This function is blank so that the activity object can be created and used
in different classes. 
*/
activity::activity()
{
 
}

void activity::compute_activities(vector<double> mindist, 
                                  vector<vector<double>> d_mindist_dx,
                                  vector<vector<double>> d_mindist_dy,
                                  vector<vector<double>> d_mindist_dz,
                                  double CC_min, double deltaCC,
                                  vector<double> apolar_contacts,
                                  vector<vector<double>> d_apolar_dx,
                                  vector<vector<double>> d_apolar_dy,
                                  vector<vector<double>> d_apolar_dz,
                                  vector<double> polar_contacts,
                                  vector<vector<double>> d_polar_dx,
                                  vector<vector<double>> d_polar_dy,
                                  vector<vector<double>> d_polar_dz,
                                  double Emin, double deltaE)
{
  
  activity_grid=vector<double>(mindist.size(),0);
  vector<double> d_activity_i(d_mindist_dx[0].size(),0);
  d_activity_dx=vector<vector<double>>(mindist.size(),d_activity_i);
  d_activity_dy=vector<vector<double>>(mindist.size(),d_activity_i);
  d_activity_dz=vector<vector<double>>(mindist.size(),d_activity_i);
  
  
  //Calculate close_contact, depth and then activity and its derivatives
  //#pragma omp parallel for
  for (unsigned i=0; i<mindist.size();i++)
  {
   S_on close_contact_i(mindist[i],CC_min,deltaCC,d_mindist_dx[i],d_mindist_dy[i],d_mindist_dz[i]);
   close_contact_i.compute_S_on();

   /* Notice here that d_apolar values are added to d_polar. 
      This is using a local copy of the object so it should 
      not affect the values of the original contacts object
   */
   double total_contacts=apolar_contacts[i]+polar_contacts[i];
   transform(d_apolar_dx[i].begin(),d_apolar_dx[i].end(),d_polar_dx[i].begin(),d_apolar_dx[i].begin(),plus<double>());
   transform(d_apolar_dy[i].begin(),d_apolar_dy[i].end(),d_polar_dy[i].begin(),d_apolar_dy[i].begin(),plus<double>());
   transform(d_apolar_dz[i].begin(),d_apolar_dz[i].end(),d_polar_dz[i].begin(),d_apolar_dz[i].begin(),plus<double>());
   S_on depth_i(total_contacts,Emin,deltaE,d_polar_dx[i],d_polar_dy[i],d_polar_dz[i]);
   depth_i.compute_S_on();

   activity_grid[i]=close_contact_i.S_on_value*depth_i.S_on_value;
   for (unsigned j=0; j<d_activity_dx[i].size(); j++)
     {
       d_activity_dx[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dx[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dx[j];
       d_activity_dy[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dy[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dy[j];
       d_activity_dz[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dz[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dz[j];
     }

  }

  //Uncomment the following lines for testing
  
  for (unsigned i=0; i<activity_grid.size();i++)
  {
    cout << "Activity point " << i << " = " << activity_grid[i] << endl;
    for (unsigned j=0; j<d_activity_dx[i].size();j++)
    {
      cout << "derivative with respect to atom: " << j << ": " << d_activity_dx[i][j] << " "
                                                               << d_activity_dy[i][j] << " "
                                                               << d_activity_dz[i][j] << " " << endl;
    }
  }
  
}
