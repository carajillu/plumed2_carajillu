#include <iostream>
#include "volume.h"
#include "activity.h"
using namespace std;

Volume::Volume()
{

}

void Volume::compute_volume(Activity &activity, double &volume_element)
{
  //Ugly way to fill vectors with zeros
  for (unsigned j=0; j<activity.d_activity_dx[0].size();j++)
  {
    d_volume_dx.push_back(0);
    d_volume_dy.push_back(0);
    d_volume_dz.push_back(0);
  }

  //computing volume
  #pragma omp parallel for 
  for (unsigned i=0; i<activity.activity.size();i++)
  {
      volume+=activity.activity[i]*volume_element;
      for (unsigned j=0; j<activity.d_activity_dx[i].size();j++)
      {
        d_volume_dx[j]+=activity.d_activity_dx[i][j]*volume_element;
        d_volume_dy[j]+=activity.d_activity_dy[i][j]*volume_element;
        d_volume_dz[j]+=activity.d_activity_dz[i][j]*volume_element;
      }
  }

  

  //Uncomment the following lines for testing
  cout << "Volume: " << volume << endl;
  for (unsigned j=0;j<d_volume_dx.size();j++)
  {
      cout << "derivative with respect to atom " << j << ": "  << d_volume_dx[j] 
                                                        << " " << d_volume_dy[j] 
                                                        << " " << d_volume_dz[j] << endl;
  }
}