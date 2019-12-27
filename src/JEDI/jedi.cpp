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
#include "distances.h"
#include "kernel.h"
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
  jediparameters params;
  getatoms all_atoms;
  getatoms grid;
  getatoms ligand;
  Activity activity;
  bool print_benchmark;
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
  pbc(true),
  print_benchmark(false)
{
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
  cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  vector<double> cog_atoms;
  all_atoms.center_atoms(all_atoms.positions,cog_atoms);
  cout << "Center of geometry went from: " << all_atoms.cog0[0] << "," << all_atoms.cog0[1] << "," << all_atoms.cog0[2] << " to " << \
                                              all_atoms.cog[0] << "," << all_atoms.cog[1] << "," << all_atoms.cog[2] << endl;
  requestAtoms(all_atoms.atomnumbers);

  //Read grid file
  cout << "*********************************" << endl;
  string pdb_grid;
  parse("GRID",pdb_grid);
  grid.readAtoms(pdb_grid);
  cout << "Loaded file " << pdb_grid << " and found " << grid.atomnumbers.size() << " elements." << endl;
  cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  grid.center_atoms(grid.positions,all_atoms.cog0);
  cout << "Center of geometry went from: " << grid.cog0[0] << "," << grid.cog0[1] << "," << grid.cog0[2] << " to " << \
                                              grid.cog[0] << "," << grid.cog[1] << "," << grid.cog[2] << endl;

  //Read reference ligand file (to be deprecated)
  string pdb_ligand;
  parse("LIGAND",pdb_ligand);
  if(pdb_ligand.length() > 0)
  {
   cout << "*********************************" << endl;
   ligand.readAtoms(pdb_ligand);
   cout << "Loaded file " << pdb_ligand << " and found " << ligand.atomnumbers.size() << " elements." << endl;
   cout << "Moving center of geometry towards origin (0,0,0)" << endl;
   ligand.center_atoms(ligand.positions,all_atoms.cog0);
   cout << "Center of geometry went from: " << ligand.cog0[0] << "," << ligand.cog0[1] << "," << ligand.cog0[2] << " to " << \
                                              ligand.cog[0] << "," << ligand.cog[1] << "," << ligand.cog[2] << endl;
  }
  cout << "*********************************" << endl;
  checkRead();
   
  /*INITIALISING COLLECTIVE VARIABLE*/
  cout << "INITIALISING JEDI COLLECTIVE VARIABLE" << endl;
  
  activity.Activity_init(grid.positions, 
                    params.CC_mind, params.deltaCC,
                    params.GP_min,params.GP_max,
                    params.CC2_min, params.deltaCC2);
  
}


// calculator
void jedi::calculate() {
  
  /////////////////////////////////////////////////
  ///               JEDI score                   //
  /////////////////////////////////////////////////
  distances distance_matrix;
  distance_matrix.compute_distance_matrix(all_atoms.positions,grid.positions);

  mindist min_dist;
  min_dist.compute_mindist(distance_matrix.r_matrix,
                          distance_matrix.dr_matrix_dx,
                          distance_matrix.dr_matrix_dy,
                          distance_matrix.dr_matrix_dz,
                          params.theta);
  
  activity.compute_activities(min_dist.min_dist, min_dist.d_mindist_dx, min_dist.d_mindist_dy, min_dist.d_mindist_dz,
                              params.CC_mind,params.deltaCC,params.GP_min,params.GP_max,params.CC2_min,params.deltaCC2,params.Emin,params.deltaE);
    
  Volume volume;
  double volume_element=pow(params.resolution,3);
  volume.compute_volume(activity,volume_element);

  Hydrophobicity hydrophobicity;
  hydrophobicity.compute_hydrophobicity(all_atoms.atomnames,distance_matrix,activity,params.r_hydro,params.deltar_hydro);

  
  


  double Jedi=12345;

  setValue(Jedi);
  
  vector<double> dJedi_dx(all_atoms.atomnumbers.size(),0);
  vector<double> dJedi_dy(all_atoms.atomnumbers.size(),0);
  vector<double> dJedi_dz(all_atoms.atomnumbers.size(),0);
  #pragma omp parallel for
  for (unsigned j=0; j<all_atoms.atomnumbers.size();j++)
  {
    setAtomsDerivatives(j,Vector(dJedi_dx[j],dJedi_dy[j],dJedi_dz[j]));
  }
  
}

}
}