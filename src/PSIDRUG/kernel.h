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
    vector<double> r;
    vector<double> dr_dx;
    vector<double> dr_dy;
    vector<double> dr_dz;
    void calculate_distance(vector<double> &gridpoint_crd, 
                            vector<PLMD::Vector> &atom_crd, 
                            vector<unsigned> &bsite_bin);

    vector<double> exp_r;
    double CC;
    vector<double> d_CC_dx;
    vector<double> d_CC_dy;
    vector<double> d_CC_dz;

    double HC;
    vector<double> d_HC_dx;
    vector<double> d_HC_dy;
    vector<double> d_HC_dz;

 public:
    kernel(unsigned &n_atoms);
    void calculate_activity(vector<double> &gridpoint_crd, 
                            vector<PLMD::Vector> &atom_crd,
                            vector<unsigned> &bsite_bin);
    double activity;
    vector<double> d_activity_dx;
    vector<double> d_activity_dy;
    vector<double> d_activity_dz;
};

#endif