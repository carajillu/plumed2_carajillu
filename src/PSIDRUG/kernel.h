/*
Code that reads the relevant PDB file(s) and loads the atom numbers, 
coordinates rescaled by a factor of 10 (so converted to nm)
and the atom name, will be used later on for apolarpolar purposes.
*/

#include <vector>

#include "tools/PDB.h"

using namespace std;

#ifndef kernel_h_
#define kernel_h_

class kernel
{
 private:
   // corefunctions
   double m_v(double v, double v0, double delta);
   double dm_dv(double delta);
   double Son_m(double m, double k);
   double dSon_dm(double m, double k);
   double Soff_m(double m, double k);
   double dSoff_dm(double m, double k);
   // parameters
    double kCCmin=pow(10,6); //hardcoded
    double CCmin;
    double deltaCC;
    double SRmin;
    double deltaSR;

    vector<double> r;
    vector<double> dr_dx;
    vector<double> dr_dy;
    vector<double> dr_dz;
    void calculate_distance(vector<double> &gridpoint_crd, 
                            vector<PLMD::Vector> &atom_crd, 
                            vector<unsigned> &bsite_bin);

    vector<double> exp_r;
    double sum_exp;
    double CC;
    double dCC_dr;
    vector<double> dCC_dx;
    vector<double> dCC_dy;
    vector<double> dCC_dz;
    void calculate_CC(vector<unsigned> &bsite_bin);

    double HC;
    vector<double> d_HC_dx;
    vector<double> d_HC_dy;
    vector<double> d_HC_dz;

 public:
    kernel(unsigned &n_atoms,
           double CCmin, double deltaCC, 
           double SRmin, double deltaSR);
    void reset();
    void calculate_activity(vector<double> &gridpoint_crd, 
                            vector<PLMD::Vector> &atom_crd,
                            vector<unsigned> &bsite_bin);
    double activity;
    vector<double> d_activity_dx;
    vector<double> d_activity_dy;
    vector<double> d_activity_dz;
};

#endif