#include <vector>

#include "tools/PDB.h"

using namespace std;

#ifndef grid_h_
#define grid_h_

class grid
{
 private:
    
     
 public:
    grid(unsigned n_atoms);
    vector<double> centre; // Centre of the grid
    void grid_read(string grid_file);
    void grid_setup(double &radius, double &spacing, unsigned &n_atoms);
    unsigned size_grid;
    //bsite_bin[j]=1 if atom j is closer than or at RSITE nm from the center of the grid
    //biste_bool[j]=0 if atom j is further away than RSITE nm from the center of the grid
    vector< vector<double>> positions; // Positions of grid points
    vector<unsigned> bsite_bin;
    void place_random(vector<PLMD::Vector> &atom_crd, double &rtol);
    void assign_bsite_bin(vector<PLMD::Vector> &atom_crd, double &rsite);
    void center_grid(vector<PLMD::Vector> &atom_crd);
    void print_grid(int id, int step);
    double PsiGrid;
    vector<double> d_Psigrid_dx;
    vector<double> d_Psigrid_dy;
    vector<double> d_Psigrid_dz;
    void init_psigrid(unsigned &n_atoms);
    void add_activity(double &activity,
                      vector<double> &d_activity_dx,
                      vector<double> &d_activity_dy,
                      vector<double> &d_activity_dz);
};

#endif
