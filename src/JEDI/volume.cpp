#include <iostream>
#include "volume.h"
#include "activity.h"
using namespace std;

Volume::Volume()
{
 volume=0;
}

void Volume::compute_volume(double sum_activity, double volume_element,
                            vector<double> d_sum_activity_dx, vector<double> d_sum_activity_dy, vector<double> d_sum_activity_dz)
{
  
  d_volume_dx=d_sum_activity_dx;
  d_volume_dy=d_sum_activity_dy;
  d_volume_dz=d_sum_activity_dz;
  
  volume=volume_element*sum_activity;
  
  for (unsigned j=0;j<d_sum_activity_dx.size();j++)
  {
    d_volume_dx[j]*=volume_element;
    d_volume_dy[j]*=volume_element;
    d_volume_dz[j]*=volume_element;
  }
  //Uncomment the following lines for testing
  /*
  cout << "Volume: " << volume << endl;
  for (unsigned j=0;j<d_volume_dx.size();j++)
  {
      cout << "derivative with respect to atom " << j << ": "  << d_volume_dx[j] 
                                                        << " " << d_volume_dy[j] 
                                                        << " " << d_volume_dz[j] << endl;
  }
  */
}