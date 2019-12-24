#include "activity.h"
#include "kernel.h"
#include <iostream>

using namespace std;

/*
This function is blank so that the activity object can be created and used
in different classes. 
*/
Activity::Activity()
{
 
}

/*
 This function initialises the activity object and takes in the coordinates of
 the grid, the parameters CCmin and deltaCC (to compute close contact), 
 GPmin and GPmax (to compute neighbours) and CC2min and deltaCC2 (to compute depth)
*/
void Activity::Activity_init(vector<PLMD::Vector> &positions,
                   double &CC_min, double deltaCC,
                   double &GPmin, double &GPmax,
                   double &CC2_min, double &deltaCC2)
{
  cout << "Computing grid neighbours" << endl;
  compute_neighbours(positions, GPmin, GPmax);
}

/*
 This function computes the neighbours of each gridpoint i as the points k
 whose distance with i belongs to the [GPmin,GPmax] range. Note that we don't
 use (distance>GPmin & distance<GPmax) as a condition. This is because the c++
 roundoff error leaves some points out. This has been tested with a small grid
 that only as 7 points, all adjacent: 1 in the middle, 2 on the sides, 1 on top,
 1 at the bottom, one in front, one at the back. The point in the middle returns 
 6 neighbours (correct) as opposed to 3 for the simpler condition form.
 */
void Activity::compute_neighbours(vector<PLMD::Vector> &positions,double &GPmin, double &GPmax)
{
 vector<int> neighbours_i;
 neighbours.assign(positions.size(),neighbours_i);
 
 for (unsigned i=0; i<positions.size();i++)
   {
     for (unsigned k=i+1; k<positions.size();k++)
     {  
      double r=delta(positions[k],positions[i]).modulo();
      if ( (r-GPmin)>0.000001 and (r-GPmax)<0.000001)
      {
        neighbours[i].push_back(k);
        neighbours[k].push_back(i);
      }
     }
   }
   int max_neighbours=0;
  for (unsigned i=0; i<positions.size();i++)
  {
   if (neighbours[i].size() > max_neighbours) max_neighbours = neighbours[i].size();
  }
  cout << "Maximum number of grid point neighbours: " << max_neighbours << endl;
}

void Activity::compute_activities(vector<double> &mindist, 
                                 vector<vector<double>> &d_mindist_dx,
                                 vector<vector<double>> &d_mindist_dy,
                                 vector<vector<double>> &d_mindist_dz,
                                double &CC_min, double &deltaCC,
                                double &GPmin, double &GPmax,
                                double &CC2_min, double &deltaCC2,
                                double &Emin, double &deltaE)
{
  
  //compute close contact
  S_on_grid close_contact;
  close_contact.compute_S_on_grid(mindist,d_mindist_dx,d_mindist_dy,d_mindist_dz,CC_min,deltaCC);
  
  
  //compute depth
  S_off_grid farawayness;
  farawayness.compute_S_off_grid(mindist,d_mindist_dx,d_mindist_dy,d_mindist_dz,CC2_min,deltaCC2);
  Activity::compute_farawayness_sum(neighbours,farawayness);
  
  S_on_grid depth;
  depth.compute_S_on_grid(depth_sum,d_depthsum_dx,d_depthsum_dy,d_depthsum_dz,Emin,deltaE);
  
  
  //compute activity and its derivatives
  
  //Ugly way to fill vectors with zeros
  vector<double> d_activity_i(d_mindist_dx[0].size(),0);
  for (unsigned i=0; i<mindist.size();i++)
  {
    activity.push_back(0);
    d_activity_dx.push_back(d_activity_i);
    d_activity_dy.push_back(d_activity_i);
    d_activity_dz.push_back(d_activity_i);
  }

  #pragma omp parallel for
  for (unsigned i=0; i<activity.size();i++)
  {
   activity[i]=close_contact.S_on_grid_vector[i]*depth.S_on_grid_vector[i];
   for (unsigned j=0; j<d_activity_dx[i].size(); j++)
     {
       d_activity_dx[i][j]=depth.S_on_grid_vector[i]*close_contact.d_Son_grid_dx[i][j]+
                           close_contact.S_on_grid_vector[i]*depth.d_Son_grid_dx[i][j];
       d_activity_dy[i][j]=depth.S_on_grid_vector[i]*close_contact.d_Son_grid_dy[i][j]+
                           close_contact.S_on_grid_vector[i]*depth.d_Son_grid_dy[i][j];
       d_activity_dz[i][j]=depth.S_on_grid_vector[i]*close_contact.d_Son_grid_dz[i][j]+
                           close_contact.S_on_grid_vector[i]*depth.d_Son_grid_dz[i][j];
     }
  }

  //Uncomment the following lines for testing
  /*
  for (unsigned i=0; i<activity.size();i++)
  {
    cout << "Activity point " << i << " = " << activity[i] << endl;
    for (unsigned j=0; j<d_activity_dx[i].size();j++)
    {
      cout << "derivative with respect to atom: " << j << ": " << d_activity_dx[i][j] << " "
                                                               << d_activity_dy[i][j] << " "
                                                               << d_activity_dz[i][j] << " " << endl;
    }
  }
  */
}
