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
    vector< vector<double>> positions; // Positions of grid points
    vector<double> centre; // Centre of the grid
    //bsite_bin[j]=1 if atom j is closer than or at RSITE nm from the center of the grid
    //biste_bool[j]=0 if atom j is further away than RSITE nm from the center of the grid
    vector<unsigned> bsite_bin; 
    double PsiGrid;
    vector<double> d_Psigrid_dx;
    vector<double> d_Psigrid_dy;
    vector<double> d_Psigrid_dz;
 public:
    grid(double &radius, double &spacing, unsigned &n_atoms);
    void place_random(vector<PLMD::Vector> &atom_crd, double &rtol);
    void assign_bsite_bin(vector<PLMD::Vector> &atom_crd, double &rsite);
    void center_grid(vector<PLMD::Vector> &atom_crd);
    void print_grid(int id, int step);
};

#endif
