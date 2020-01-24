#include <vector>
#include "tools/PDB.h"

using namespace std;

#ifndef jedi_kernel_h_
#define jedi_kernel_h_

class S_on
{
 private:
    double pi;
 public:
    S_on();
    double k;
    double S_on_value;
    vector<double> d_Son_dx;
    vector<double> d_Son_dy;
    vector<double> d_Son_dz;
    void compute_S_on(double v, double v0, double delta_v, vector<double> dv_dx, vector<double> dv_dy, vector<double> dv_dz);
};

class S_off
{
 private:
    double pi;
 public:
    S_off();
    double k;
    double S_off_value;
    vector<double> d_Soff_dx;
    vector<double> d_Soff_dy;
    vector<double> d_Soff_dz;
    void compute_S_off(double v, double v0, double delta_v, vector<double> dv_dx, vector<double> dv_dy, vector<double> dv_dz);
};

#endif