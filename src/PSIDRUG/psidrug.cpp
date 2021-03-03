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
#include <iostream>
#include <cmath>

//CV modules
#include "grid.h"
#include "kernel.h"

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
  string grid_file;
  unsigned ngrid=0;
  double rgrid=0;
  double spacing=0;
  double rtol=0;
  double rsite=0;
  unsigned n_atoms;
  vector<grid> grids;
  vector<kernel>kernels;
  double PsiDrug;
  vector<double> d_PsiDrug_dx;
  vector<double> d_PsiDrug_dy;
  vector<double> d_PsiDrug_dz;
  int numstep=0;

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
  keys.add("optional","GRID_FILE","XYZ file containing a sample grid. Used for debug purposes.");
}

Psidrug::Psidrug(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  debug(false)
{
  parseFlag("DEBUG",debug);
  if (debug)
  {
     log.printf("RUNNING IN DEBUG MODE\n");
     parse("GRID_FILE",grid_file);
     if(grid_file=="")
        {
          log.printf("Please specify a grid xyz file if you are running in debug mode");
          exit(0);
        }
  }
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

  addValueWithDerivatives(); 
  setNotPeriodic();

  requestAtoms(atoms);
  n_atoms=atoms.size();

  // Setup grid(s)
  for (unsigned i=0; i<ngrid; i++)
  {
   grids.push_back(grid(n_atoms));
   if (debug)
   {
    grids[i].grid_read(grid_file);
   }
   else
   {
   grids[i].grid_setup(rgrid,spacing,n_atoms);
   }
   kernels.push_back(kernel(n_atoms));
  }

  //Initialise Psidrug and its derivatives
  PsiDrug=0;
  d_PsiDrug_dx=vector<double>(n_atoms,0);
  d_PsiDrug_dy=vector<double>(n_atoms,0);
  d_PsiDrug_dz=vector<double>(n_atoms,0);
}


// calculator
void Psidrug::calculate() {
  if (pbc) makeWhole();
  vector<Vector> atom_crd=getPositions();
  int step=getStep();
  numstep++;
  PsiDrug=0;

  for (unsigned k=0; k<grids.size();k++)
  {
    //Update Grid Coordinates//
    if (step==0)
    {
      if (!debug) grids[k].place_random(atom_crd,rtol);
      grids[k].assign_bsite_bin(atom_crd,rsite);
    }
    /*
    This if block is here just to debug derivatives
    The grid update bugs when trying to calculate the derivatives of the box
    Need to ask (Giovanni?) how to get around that and what impact it may
    have on simulations
    */
    if (numstep<=n_atoms*3 and checkNumericalDerivatives())
    {
      grids[k].center_grid(atom_crd);
      grids[k].assign_bsite_bin(atom_crd,rsite);
    }
 
    //if (debug) 
    grids[k].print_grid(k,step);
    grids[k].reset_psigrid();
    for (unsigned i=0; i<grids[k].size_grid;i++)
    {
     kernels[k].calculate_activity(grids[k].positions[i],atom_crd,grids[k].bsite_bin);
     grids[k].add_activity(kernels[k].activity,
                           kernels[k].d_activity_dx,
                           kernels[k].d_activity_dy,
                           kernels[k].d_activity_dz);
    }
    PsiDrug+=grids[k].PsiGrid;
    for (unsigned j=0;j<n_atoms;j++)
    {
      d_PsiDrug_dx[j]+=grids[k].d_Psigrid_dx[j];
      d_PsiDrug_dy[j]+=grids[k].d_Psigrid_dy[j];
      d_PsiDrug_dz[j]+=grids[k].d_Psigrid_dz[j];
    }
  }
  
  for (unsigned j=0; j<n_atoms;j++)
  {
    setAtomsDerivatives(j,Vector(d_PsiDrug_dx[j],d_PsiDrug_dy[j],d_PsiDrug_dz[j]));
  }
  setBoxDerivativesNoPbc();
  setValue(PsiDrug);
  cout << "Step " << numstep << "PsiDrug = " << PsiDrug << endl;
}

}
}



