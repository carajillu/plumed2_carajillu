#include "activity.h"
#include <iostream>

using namespace std;

Activity::Activity()
{
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
     //cout << "Grid point " << i << " looking at :" << endl;
     for (unsigned k=i+1; k<positions.size();k++)
     {
      //cout << " " << k << " distance: ";
      double r=delta(positions[k],positions[i]).modulo();
      //cout  << r << ", GPmin: " << GPmin << ", GPmax: " << GPmax;
      if ( (r-GPmin)>0.000001 and (r-GPmax)<0.000001)
      {
        //cout << "--adding point ";
        neighbours[i].push_back(k);
        neighbours[k].push_back(i);
      }
      //cout << endl;
     }
   }
}
