#include <vector>

using namespace std;

#ifndef corefunctions_h_
#define corefunctions_h_

class corefunctions
{
 private:
 public:
  corefunctions();

  double m_v(double v, double v0, double delta);
  double dm_dv(double delta);

  double Son_m(double m);
  double dSon_dm(double m);

  double Soff_m(double m);
  double dSoff_dm(double m);
};
#endif