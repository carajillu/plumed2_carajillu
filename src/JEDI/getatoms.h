/*
Code that reads the relevant PDB file(s) and loads the atom numbers, 
coordinates rescaled by a factor of 10 (so converted to nm)
and the atom name, will be used later on for apolarpolar purposes.
*/

#include <vector>

#include "tools/PDB.h"

using namespace std;

#ifndef atomlist_h_
#define atomlist_h_

class getatoms
{
 private:
    vector<PLMD::AtomNumber> atomlist;
 public:
    getatoms();
    void readAtoms(string pdb_file);
    void center_atoms(vector<PLMD::Vector> &positions, vector<double> cog_in);
    void select_atoms(vector<PLMD::Vector> &positions, vector<PLMD::Vector> grid_positions, vector<string> &atomnames, double r_max);
    vector<PLMD::AtomNumber> atomnumbers;
    vector<PLMD::Vector> positions;
    vector<string> atomnames;
    vector<double> cog;
    vector<double> cog0;
    vector<unsigned> atoms_jedi;
};

#endif
