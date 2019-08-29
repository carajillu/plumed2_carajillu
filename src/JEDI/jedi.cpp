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

#include "jediparameters.h" //Class to read JEDI parameters
#include "getatoms.h"
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
}

jedi::jedi(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  //READ jedi.parameters here
  string parameters_file;
  parse("PARAMETERS",parameters_file);
  params.readParams(parameters_file);

  //Read binding site file
  string pdb_bsite;
  parse("BINDINGSITE",pdb_bsite);
  all_atoms.readAtoms(pdb_bsite);

  //Read grid file
  string pdb_grid;
  parse("GRID",pdb_grid);
  grid.readAtoms(pdb_grid);


  addValue(); // to be replaced by AddValueWithDerivatives()
  setNotPeriodic();
  checkRead();
}


// calculator
void jedi::calculate() {

  double Jedi=12345.0;
  setValue(Jedi);
  double dJedi=1.0;
}

}
}