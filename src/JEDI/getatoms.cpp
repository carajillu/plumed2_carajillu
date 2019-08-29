/*
Code that reads the relevant PDB file(s) and loads the atom numbers, 
coordinates rescaled by a factor of 10 (so converted to nm)
and the atom name, will be used later on for apolarpolar purposes.
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "getatoms.h"

using namespace std;

getatoms::getatoms()
{   
}

bool getatoms::readAtoms(string &pdb_file)
  {
    PLMD::PDB pdb_handle;
    if( !pdb_handle.read(pdb_file,false,0.1)) //PDB files are always in angstroms, so we need to rescale the coordinates by a factor of 10
      {
        cout << "File: " << pdb_file.c_str() << " not found. Exiting." << endl; //JCN Aug2019: Needs to crash if the file is malformed. Also needs to do it with a plumed error.
        exit(0);
      }
    atomnumbers = pdb_handle.getAtomNumbers();
    positions = pdb_handle.getPositions();
    for (int j=0; j<pdb_handle.getAtomNumbers().size(); j++)
      {
       atomnames.push_back(pdb_handle.getAtomName(atomnumbers[j]));
      }
    
    cout << "Loaded file " << pdb_file << " and found " << atomnumbers.size() << " elements." << endl;
    return true;
  }