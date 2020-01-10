/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Matrix.h"
#include "tools/PDB.h"

#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

// Lapack needed for l2-mininum norm solution
#include "../tools/lapack/lapack.h"

// Classes necessary for reading in input files
#include "jediparameters.h" //Class to read JEDI parameters
#include "getatoms.h" // Class to read in atoms and grid

//Classes necessary for initialisation and calculation of the Collective variable
#include "kernel.h"
#include "distances.h"
#include "contacts.h"
#include "activity.h"
#include "volume.h"
#include "hydrophobicity.h"


#include "kabsch.h" // kabsch algorithm for grid update


typedef double real;
using namespace std;


namespace PLMD {
namespace colvar {

//Here you can put some documentation if you want

class jedi : public Colvar {
private:  
  bool pbc;
  int iteration = 0;
  jediparameters params;
  getatoms all_atoms;
  getatoms grid;
  getatoms ligand;
public:
  explicit jedi(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(jedi,"JEDI")

void jedi::registerKeywords(Keywords& keys) 
{
  Colvar::registerKeywords(keys);
  keys.add("compulsory","PARAMETERS","a file listing the parameters of the JEDI estimator.");
  keys.add("compulsory","BINDINGSITE","pdbfile containing the atoms of the binding site");
  keys.add("compulsory","GRID","PDB file containing the JEDI grid");
  keys.add("optional","LIGAND","PDB file containing a known ligand (Optional)");
}

jedi::jedi(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  #pragma omp parallel
  {
   #pragma omp single
    {
     int nthreads=omp_get_num_threads();
     cout << "JEDI is running with " << nthreads << " OMP threads" << endl;
    }
   #pragma omp critical
    {
     int thread_id=omp_get_thread_num();
     cout << "This is thread number " << thread_id << endl;
    }
  }

  /* INITIALISING CV AND DERIVATIVES */
  addValueWithDerivatives(); //Developers' note: this goes before requestAtoms()
  setNotPeriodic();

  /*READING IN INPUT FILES*/

  //READ jedi.parameters here
  string parameters_file;
  parse("PARAMETERS",parameters_file);
  params.readParams(parameters_file);

  //Read binding site file
  cout << "*********************************" << endl;
  string pdb_bsite;
  parse("BINDINGSITE",pdb_bsite);
  all_atoms.readAtoms(pdb_bsite);
  cout << "Loaded file " << pdb_bsite << " and found " << all_atoms.atomnumbers.size() << " elements." << endl;
  //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  //vector<double> cog_atoms;
  //all_atoms.center_atoms(all_atoms.positions,cog_atoms);
  //cout << "Center of geometry went from: " << all_atoms.cog0[0] << "," << all_atoms.cog0[1] << "," << all_atoms.cog0[2] << " to " << \
                                              all_atoms.cog[0] << "," << all_atoms.cog[1] << "," << all_atoms.cog[2] << endl;
  requestAtoms(all_atoms.atomnumbers);

  //Read grid file
  cout << "*********************************" << endl;
  string pdb_grid;
  parse("GRID",pdb_grid);
  grid.readAtoms(pdb_grid);
  cout << "Loaded file " << pdb_grid << " and found " << grid.atomnumbers.size() << " elements." << endl;
  //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  //grid.center_atoms(grid.positions,all_atoms.cog0);
  //cout << "Center of geometry went from: " << grid.cog0[0] << "," << grid.cog0[1] << "," << grid.cog0[2] << " to " << \
                                              grid.cog[0] << "," << grid.cog[1] << "," << grid.cog[2] << endl;

  //Read reference ligand file (to be deprecated)
  string pdb_ligand;
  parse("LIGAND",pdb_ligand);
  if(pdb_ligand.length() > 0)
  {
   cout << "*********************************" << endl;
   ligand.readAtoms(pdb_ligand);
   cout << "Loaded file " << pdb_ligand << " and found " << ligand.atomnumbers.size() << " elements." << endl;
   //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
   //ligand.center_atoms(ligand.positions,all_atoms.cog0);
   //cout << "Center of geometry went from: " << ligand.cog0[0] << "," << ligand.cog0[1] << "," << ligand.cog0[2] << " to " << \
                                              ligand.cog[0] << "," << ligand.cog[1] << "," << ligand.cog[2] << endl;
  }
  cout << "*********************************" << endl;
  checkRead();
   
  /*INITIALISING COLLECTIVE VARIABLE*/
  cout << "INITIALISING JEDI COLLECTIVE VARIABLE" << endl;
}


// calculator
void jedi::calculate() {
  if (pbc) makeWhole();

  /////////////////////////////////////////////////
  // Update Grid coordinates and Binding site    //
  /////////////////////////////////////////////////

  /*
  The cubic grid is updated by simply centering it
  in the COG of a rigid part of the ligand
  */
  
  /*
  Shrink grid and binding site so that only atoms/points that
  are within interacting distance from each other are kept
  */
  cout << "Shrinking binding site" << endl;
  all_atoms.positions.clear();
  for (unsigned j=0; j<all_atoms.atomnumbers.size();j++)
  {
    all_atoms.positions.push_back(getPosition(j));
  }
  all_atoms.select_atoms(all_atoms.positions,grid.positions,all_atoms.atomnames,params.r_max);
  cout << "Kept " << all_atoms.positions.size() << " atoms" << endl;
  //grid.select_atoms(grid.positions,all_atoms.positions,grid.atomnames,params.r_max); //NOT YET
  
  
  /////////////////////////////////////////////////
  //                JEDI score                   //
  /////////////////////////////////////////////////
  cout << "calculating distance matrix" << endl;
  distances distance_matrix;
  distance_matrix.compute_distance_matrix(all_atoms.positions,grid.positions);

  cout << "mindist" << endl;
  mindist min_dist;
  min_dist.compute_mindist(distance_matrix.r_matrix,
                          distance_matrix.dr_matrix_dx,
                          distance_matrix.dr_matrix_dy,
                          distance_matrix.dr_matrix_dz,
                          params.theta);

  cout << "calculating contacts" << endl;                        
  contacts_Soff contacts;
  contacts.compute_contacts_S_off(distance_matrix.r_matrix,
                                  distance_matrix.dr_matrix_dx,
                                  distance_matrix.dr_matrix_dy,
                                  distance_matrix.dr_matrix_dz,
                                  params.r_hydro, params.deltar_hydro);

 cout << "calculating contacts_sum" << endl;
  contacts_sum contacts_sum;
  contacts_sum.compute_contacts_sum(contacts.contacts_matrix,
                                    contacts.d_contacts_dx,
                                    contacts.d_contacts_dy,
                                    contacts.d_contacts_dz,
                                    all_atoms.atomnames);

  cout << "calculating activity" << endl;
  activity activity;
  activity.compute_activities(min_dist.min_dist, 
                              min_dist.d_mindist_dx, min_dist.d_mindist_dy, min_dist.d_mindist_dz,
                              params.CC_mind,params.deltaCC,
                              contacts_sum.contacts_apolar,
                              contacts_sum.d_contacts_apolar_dx,contacts_sum.d_contacts_apolar_dy,contacts_sum.d_contacts_apolar_dz,
                              contacts_sum.contacts_polar,
                              contacts_sum.d_contacts_polar_dx,contacts_sum.d_contacts_polar_dy,contacts_sum.d_contacts_polar_dz,
                              params.Emin, params.deltaE);

 cout << "calculating volume" << endl;
  Volume volume;
  double volume_element=pow(params.resolution,3);
  volume.compute_volume(activity.sum_activity,volume_element,
                        activity.d_sum_activity_dx,
                        activity.d_sum_activity_dy,
                        activity.d_sum_activity_dz);

  cout << "hydrophobicity" << endl;
  hydrophobicity hydrophobicity;
  hydrophobicity.compute_hydrophobicity(contacts_sum.contacts_apolar,
                                        contacts_sum.d_contacts_apolar_dx,contacts_sum.d_contacts_apolar_dy,contacts_sum.d_contacts_apolar_dz,
                                        contacts_sum.contacts_polar,
                                        contacts_sum.d_contacts_polar_dx,contacts_sum.d_contacts_polar_dy,contacts_sum.d_contacts_polar_dz,
                                        activity.activity_grid,
                                        activity.d_activity_dx,activity.d_activity_dy,activity.d_activity_dz,
                                        activity.sum_activity,
                                        activity.d_sum_activity_dx,activity.d_sum_activity_dy,activity.d_sum_activity_dz);
  
  cout << "JEDI" << endl;
  //double Jedi=params.alpha*volume.volume/params.V_max+params.beta*hydrophobicity.Ha+params.gamma;
  double Jedi=hydrophobicity.Ha;

  setValue(Jedi);

  cout << "DERIVATIVES" << endl;
  vector<double> dJedi_dx(all_atoms.atomnumbers.size(),0);
  vector<double> dJedi_dy(all_atoms.atomnumbers.size(),0);
  vector<double> dJedi_dz(all_atoms.atomnumbers.size(),0);

  
  #pragma omp parallel for
  for (unsigned j=0; j<all_atoms.atoms_jedi.size();j++)
  {
    unsigned atom_idx=all_atoms.atoms_jedi[j];
    dJedi_dx[atom_idx]=hydrophobicity.d_Ha_dx[j];
    dJedi_dy[atom_idx]=hydrophobicity.d_Ha_dy[j];
    dJedi_dz[atom_idx]=hydrophobicity.d_Ha_dz[j];
    //dJedi_dx[atom_idx]=params.alpha*volume.d_volume_dx[j]/params.V_max+params.beta*hydrophobicity.d_Ha_dx[j];
    //Jedi_dy[atom_idx]=params.alpha*volume.d_volume_dy[j]/params.V_max+params.beta*hydrophobicity.d_Ha_dy[j];
    //dJedi_dz[atom_idx]=params.alpha*volume.d_volume_dz[j]/params.V_max+params.beta*hydrophobicity.d_Ha_dz[j];
  }

  #pragma omp parallel for
  for (unsigned j=0; j<all_atoms.atomnumbers.size();j++)
  {
    setAtomsDerivatives(j,Vector(dJedi_dx[j],dJedi_dy[j],dJedi_dz[j]));
  }

  //iteration++;
  //cout << "Iteration number: " << iteration;
  //cout << "Volume = " << volume.volume << " Hydrophobicity = " << hydrophobicity.Ha << " JEDI = " << Jedi << endl;
}

}
}