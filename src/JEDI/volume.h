#include "activity.h"

using namespace std;

#ifndef volume_h_
#define volume_h_

class Volume
{
 public:
  Volume();
  double volume;
  vector<double> d_volume_dx;
  vector<double> d_volume_dy;
  vector<double> d_volume_dz;
  void compute_volume(Activity &activity, double &volume_element);
};
#endif