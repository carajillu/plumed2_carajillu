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
 public:
    getatoms();
    bool readAtoms(string &pdb_file);
    vector<PLMD::AtomNumber> atomnumbers;
    vector<PLMD::Vector> positions;
    vector<string> atomnames;
};

#endif
