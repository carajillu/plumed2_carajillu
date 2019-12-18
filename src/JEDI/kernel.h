#include <vector>

using namespace std;

#ifndef jedi_kernel_h_
#define jedi_kernel_h_

class kernel
{
 public:
   kernel();
   double mindist();
   vector<double> d_mindist;
   
   double m;
   vector<double> d_m;

   double S_on;
   vector<double> d_S_on;

   double S_off;
   vector<double> d_S_off;
};
#endif