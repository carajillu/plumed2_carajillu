#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

using namespace std;
typedef double real;

#include <vector>
#include <cmath>
#include <iostream>

#ifndef distance_h

#define distance_h

double r(vector <double> & vec1, vector <double> & vec2)
{
    if (vec1.size()!=vec2.size())
    {
        cout << "Vectors are not the same size. Can't calculate distance.";
        exit(0);
    }
    double r2=0;
    #pragma omp parallel for reduction(+:sum)
    for(unsigned i=0;i<vec1.size();i++)
    {
        r2+=(vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }
    double r=sqrt(r2);
    
    return r;
}

double dr_di(double r, unsigned i, vector<double> & vec1, vector<double> & vec2)
{
  /* 
  The distance is an absolute value but its derivative with respect to one component is a vector
  Therefore its sign depends on which point is considered the one of reference
  */
  double dr=(vec1[i]-vec2[i])/r;
  return dr;
}
#endif