/*
Code that reads the relevant PDB file(s) and loads the atom numbers, 
coordinates rescaled by a factor of 10 (so converted to nm)
and the atom name, will be used later on for apolarpolar purposes.
*/

#include <vector>

#include "tools/PDB.h"

using namespace std;

#ifndef grid_h_
#define grid_h_

class grid
{
 private:
    vector< vector<double>> positions;
    vector<double> centre;
    double PsiGrid;
    vector<double> d_Psigrid_dx;
    vector<double> d_Psigrid_dy;
    vector<double> d_Psigrid_dz;
 public:
    grid(double &radius, double &spacing);
    void place_random(vector<PLMD::Vector> &atoms);
    void print_grid(int id, int step);
};

#endif
