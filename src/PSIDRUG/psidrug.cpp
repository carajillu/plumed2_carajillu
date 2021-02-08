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
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"

#include <string>
#include <cmath>

//CV modules
#include "grid.h"

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR TEMPLATE
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

\plumedfile
# This should be a sample input.
t: TEMPLATE ATOMS=1,2
PRINT ARG=t STRIDE=100 FILE=COLVAR
\endplumedfile
<!---You should reference here the other actions used in this example--->
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class Psidrug : public Colvar {
  bool pbc;
  bool debug;
  unsigned ngrid=0;
  double rgrid=0;
  double spacing=0;
  double rtol=0;
  double rsite=0;
  vector<grid> grids;

public:
  explicit Psidrug(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Psidrug,"PSIDRUG")

void Psidrug::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.addFlag("DEBUG",false,"Running in debug mode");
  keys.add("atoms","ATOMS","Atoms to include in druggability calculations (start at 1)");
  keys.add("optional","NGRID","Number of quasi-spherical grids to place in the system (default = 1)");
  keys.add("optional","RGRID","Radius of the quasi-spherical grids that will be placed in the system (default = 0.3 nm)");
  keys.add("optional","SPACING","Space between adjacent grid points (default = 0.1 nm)");
  keys.add("optional","RTOL","Allowed maximum distance between the protein and the center of the grid at setup (default = 0.3 nm)");
  keys.add("optional","RSITE","Allowed maximum distance from the edge of the grid for any given atom to be taken into account (default = 0.45 nm)");
}

Psidrug::Psidrug(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  debug(false)
{
  parseFlag("DEBUG",debug);
  if (debug)
     log.printf("RUNNING IN DEBUG MODE\n");
  
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);

  parse("NGRID",ngrid);
  if (!ngrid) ngrid=1;

  parse("RGRID",rgrid);
  if (!rgrid) rgrid=0.3;

  parse("SPACING",spacing);
  if (!spacing) spacing=0.1;

  parse("RTOL",rtol);
  if (!rtol) rtol=0.3;

  parse("RSITE",rsite);
  if (!rsite) rsite=0.45;
  // Apply correction to rsite so that bsite_bin can be measured from the center
  rsite+=rgrid;

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  log.printf("  using atoms %d to %d\n",atoms[0].serial(),atoms[atoms.size()-1].serial());
  log.printf("  using %d quasi-spherical grid(s)\n",ngrid);
  log.printf("  of radius equal to %f nm\n",rgrid);
  log.printf("  and spacing between adjacent points equal to %f nm\n",spacing);
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  requestAtoms(atoms);
  unsigned n_atoms=atoms.size();

  for (unsigned i=0; i<ngrid; i++)
  {
   grids.push_back(grid(rgrid,spacing,n_atoms));
  }
}


// calculator
void Psidrug::calculate() {
  if (pbc) makeWhole();
  vector<Vector> atom_crd=getPositions();

  int step=getStep();
  
  for (unsigned i=0; i<grids.size();i++)
  {
    if (step==0)
    {
     grids[i].place_random(atom_crd,rtol);
     grids[i].assign_bsite_bin(atom_crd,rsite);
    }
    else
    {
     grids[i].center_grid(atom_crd);
     grids[i].assign_bsite_bin(atom_crd,rsite);
    }
    grids[i].print_grid(i,step);
  }

  //setAtomsDerivatives(0,-invvalue*distance);
  //setAtomsDerivatives(1,invvalue*distance);
  //setBoxDerivatives  (-invvalue*Tensor(distance,distance));
  //setValue           (value);
}

}
}



