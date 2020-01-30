#include "activity.h"
#include "kernel.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <iomanip> // std::setprecision
#include <fstream> // std::ofstream 

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
                                  vector<double> total_contacts,
                                  vector<vector<double>> d_contacts_total_dx,
                                  vector<vector<double>> d_contacts_total_dy,
                                  vector<vector<double>> d_contacts_total_dz,
                                  double Emin, double deltaE)
{
  S_on_contacts=vector<double>(mindist.size(),0);
  S_on_mindist=vector<double>(mindist.size(),0);
  activity_grid=vector<double>(mindist.size(),0);
  vector<double> d_activity_i(d_mindist_dx[0].size(),0);
  d_activity_dx=vector<vector<double>>(mindist.size(),d_activity_i);
  d_activity_dy=vector<vector<double>>(mindist.size(),d_activity_i);
  d_activity_dz=vector<vector<double>>(mindist.size(),d_activity_i);
  
  sum_activity=0;
  d_sum_activity_dx=vector<double>(d_mindist_dx[0].size(),0);
  d_sum_activity_dz=vector<double>(d_mindist_dx[0].size(),0);
  d_sum_activity_dy=vector<double>(d_mindist_dx[0].size(),0);
  
  //Calculate close_contact, depth and then activity and its derivatives
  //#pragma omp parallel for
  for (unsigned i=0; i<mindist.size();i++)
  {
   S_on close_contact_i(mindist[i],CC_min,deltaCC,d_mindist_dx[i],d_mindist_dy[i],d_mindist_dz[i]);
   close_contact_i.compute_S_on();
   S_on_mindist[i]=close_contact_i.S_on_value;

   S_on depth_i(total_contacts[i],Emin,deltaE,d_contacts_total_dx[i],d_contacts_total_dy[i],d_contacts_total_dz[i]);
   depth_i.compute_S_on();
   S_on_contacts[i]=depth_i.S_on_value;

   activity_grid[i]=close_contact_i.S_on_value*depth_i.S_on_value;
   sum_activity+=activity_grid[i];
   for (unsigned j=0; j<d_activity_dx[i].size(); j++)
     {
       d_activity_dx[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dx[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dx[j];
       d_sum_activity_dx[j]+=d_activity_dx[i][j];

       d_activity_dy[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dy[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dy[j];
       d_sum_activity_dy[j]+=d_activity_dy[i][j];

       d_activity_dz[i][j]=depth_i.S_on_value*close_contact_i.d_Son_dz[j]+
                           close_contact_i.S_on_value*depth_i.d_Son_dz[j];
       d_sum_activity_dz[j]+=d_activity_dz[i][j];                    
     }
  }


  //Uncomment the following lines for testing
  /*
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
  */
}

/*
This removes the grid points that have activity=0, but not used at the moment
*/
void activity::filter_activities()
{
   vector<double> activity_grid_filtered;
   vector<vector<double>> d_activity_dx_filtered;
   vector<vector<double>> d_activity_dy_filtered;
   vector<vector<double>> d_activity_dz_filtered;
   for (unsigned i=0; i<activity_grid.size();i++)
   {
     if (activity_grid[i]==0) continue;
     activity_grid_filtered.push_back(activity_grid[i]);
     d_activity_dx_filtered.push_back(d_activity_dx[i]);
     d_activity_dy_filtered.push_back(d_activity_dy[i]);
     d_activity_dz_filtered.push_back(d_activity_dz[i]);
     active_indices.push_back(i);
   }
   activity_grid=activity_grid_filtered;
   d_activity_dx=d_activity_dx_filtered;
   d_activity_dy=d_activity_dy_filtered;
   d_activity_dz=d_activity_dz_filtered;
}

void activity::print_activities()
{
      string filename = "activity";
      stringstream num;
      filename.append("-step-");
      filename.append(".xyz");
      ofstream wfile;
      wfile.open(filename.c_str());
      wfile << "Point Activity Close_contact Depth" << endl;
      for (unsigned i=0; i<activity_grid.size();i++)
        {
          wfile << std::fixed << std::setprecision(5) << i << " " << activity_grid[i] << " " << S_on_mindist[i] << " " << S_on_contacts[i] << endl;
       }
      wfile.close();
}