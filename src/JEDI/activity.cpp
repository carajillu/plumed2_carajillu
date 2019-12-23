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
                                double &CC2_min, double &deltaCC2)
{
  cout << "computing close contact"<<endl;
  compute_close_contact(mindist,d_mindist_dx,d_mindist_dy,d_mindist_dz,CC_min,deltaCC);
  exit(0);
          
}
