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

void getatoms::readAtoms(string &pdb_file)
  {
    PLMD::PDB pdb_handle;
    if( !pdb_handle.read(pdb_file,false,0.1)) //PDB files are always in angstroms, so we need to rescale the coordinates by a factor of 10
      {
        cout << "File: " << pdb_file.c_str() << " not found. Exiting." << endl; //JCN Aug2019: Needs to crash if the file is malformed. Also needs to do it with a plumed error.
        exit(0);
      }
    atomlist = pdb_handle.getAtomNumbers();
    
    cout << "Entering atomnames loop " << endl;
    for (int j=0; j<atomlist.size(); j++)
      {
       if (pdb_handle.getAtomName(atomlist[j])[0]=='C' or pdb_handle.getAtomName(atomlist[j])[0]=='S' or
           pdb_handle.getAtomName(atomlist[j])[0]=='N' or pdb_handle.getAtomName(atomlist[j])[0]=='O' or
           pdb_handle.getAtomName(atomlist[j])=="GRI" )
           {
              atomnames.push_back(pdb_handle.getAtomName(atomlist[j]));
              atomnumbers.push_back(atomlist[j]);
              positions.push_back(pdb_handle.getPositions()[j]);
           }
      }
      cout << "Exiting atomnames loop " << endl;
  }

/*
This function removes a center of geometry from an atom group (not necessary its own).

If this center of geometry is not provided (vector cog has size 0), the function will calculate the cog
of the atom group and remove it. The original cog will be stored in the vector cog0.

If the center of geometry is provided (cog has size 3). The function will remove it from the coordinates of
the atom group. The original cog of the atom group will still be calculated and stored in cog0.

This is done because we want to center the binding site on the origin of coordinates, but want to preserve the 
relative location of the grid and the ligand (if any) with respect to the binding site, so grid and ligand are
potentially centered on a point slightly displaced from the origin.
*/
void getatoms::center_atoms(vector<PLMD::Vector> &positions, vector<double> &cog_in)
  {
    
    /* Depending on the size of cog_in, we decide if we have to calculate a
       center of geometry to remove or if we are going to remove an existing one.
       If it is the existing one, we'll save the calculated one as cog0
    */
    if (cog_in.size()==0)
    {
     cog_in.push_back(0);
     cog_in.push_back(0);
     cog_in.push_back(0);
     for (unsigned j=0; j<positions.size();j++)
      {
       cog_in[0] += positions[j][0]/positions.size();
       cog_in[1] += positions[j][1]/positions.size();
       cog_in[2] += positions[j][2]/positions.size();
      }
     cout << "Removing own center of geometry: "<< cog_in[0] << " " << cog_in[1] << " " << cog_in[2] << endl;
    }
    else if (cog_in.size() == 3)
    {
     cout << "Removing alien of geometry: "<< cog_in[0] << " " << cog_in[1] << " " << cog_in[2] << endl;
     cog0.push_back(0);
     cog0.push_back(0);
     cog0.push_back(0);
     for (unsigned j=0; j<positions.size();j++)
      { 
       cog0[0] += positions[j][0]/positions.size();
       cog0[1] += positions[j][1]/positions.size();
       cog0[2] += positions[j][2]/positions.size();
      }
    }
    else
    {
     cout << "Fatal Error in getatoms::center_atoms. cog has dimension different than 0 or 3. Exiting.";
     exit(0); 
    }
    
    // Here we remove the supplied or calculated center of geometry
    for (unsigned j=0; j<positions.size();j++)
     {
      positions[j][0] -= cog_in[0];
      positions[j][1] -= cog_in[1];
      positions[j][2] -= cog_in[2];
     }
    
    // Here, if the center of geometry was calculated, we'll save it as cog0
    if (cog0.size()==0)
     {
      cog0.push_back(cog_in[0]);
      cog0.push_back(cog_in[1]);
      cog0.push_back(cog_in[2]);
     }
    
    
    // Here we calculate the new cog (should be close to 0)
    cog.push_back(0.);
    cog.push_back(0.);
    cog.push_back(0.);
    for (unsigned j=0; j<positions.size();j++)
       {
        cog[0] += positions[j][0]/positions.size();
        cog[1] += positions[j][1]/positions.size();
        cog[2] += positions[j][2]/positions.size();
       }
   }