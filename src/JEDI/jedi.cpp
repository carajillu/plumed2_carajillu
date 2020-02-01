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
#include "clustering_laio.h"


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
  getatoms protein;
  getatoms grid;
  getatoms ligand;
  int iteration=0;
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

  /* OPENING STATS FILE */
  ofstream summary;
  summary.open("jedi_stats.dat");
  summary << "Step sum_activity Va Ha JEDI" << endl;
  summary.close();

  /*READING IN INPUT FILES*/

  //READ jedi.parameters here
  string parameters_file;
  parse("PARAMETERS",parameters_file);
  params.readParams(parameters_file);

  //Read binding site file
  cout << "*********************************" << endl;
  string pdb_bsite;
  parse("BINDINGSITE",pdb_bsite);
  protein.readAtoms(pdb_bsite);
  cout << "Loaded file " << pdb_bsite << " and found " << protein.atomnumbers.size() << " elements." << endl;
  //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  //vector<double> cog_atoms;
  //protein.center_atoms(protein.positions,cog_atoms);
  //cout << "Center of geometry went from: " << protein.cog0[0] << "," << protein.cog0[1] << "," << protein.cog0[2] << " to " << \
                                              protein.cog[0] << "," << protein.cog[1] << "," << protein.cog[2] << endl;

                                                //Read reference ligand file
  string pdb_ligand;
  parse("LIGAND",pdb_ligand);
  if(pdb_ligand.length() > 0)
  {
   cout << "*********************************" << endl;
   ligand.readAtoms(pdb_ligand);
   ligand.compute_cog(ligand.positions);
   cout << "Loaded file " << pdb_ligand << " and found " << ligand.atomnumbers.size() << " elements." << endl;
   //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
   //ligand.center_atoms(ligand.positions,protein.cog0);
   //cout << "Center of geometry went from: " << ligand.cog0[0] << "," << ligand.cog0[1] << "," << ligand.cog0[2] << " to " << \
                                              ligand.cog[0] << "," << ligand.cog[1] << "," << ligand.cog[2] << endl;
  }
  
  vector<PLMD::AtomNumber> atomstorequest;
  for (unsigned j=0; j<protein.atomnumbers.size();j++)
  {
    atomstorequest.push_back(protein.atomnumbers[j]);
  }
  for (unsigned k=0; k<ligand.atomnumbers.size();k++)
  {
    atomstorequest.push_back(ligand.atomnumbers[k]);
  }

  requestAtoms(atomstorequest); //Developer's note: this goes AFTER protein and ligand because we want to request the atoms of the ligand too

  //Read grid file
  cout << "*********************************" << endl;
  string pdb_grid;
  parse("GRID",pdb_grid);
  grid.readAtoms(pdb_grid);
  cout << "Loaded file " << pdb_grid << " and found " << grid.atomnumbers.size() << " elements." << endl;
  grid.compute_neighbours(grid.positions,params.GP_min);
  int max_neighbours=0;
  for (unsigned i=0; i<grid.positions.size();i++)
  {
   if (grid.neighbours[i].size()>max_neighbours) max_neighbours=grid.neighbours[i].size();
  }
  cout << "Maximum number of neighbours found: " << max_neighbours << endl;
  //cout << "Moving center of geometry towards origin (0,0,0)" << endl;
  //grid.center_atoms(grid.positions,protein.cog0);
  //cout << "Center of geometry went from: " << grid.cog0[0] << "," << grid.cog0[1] << "," << grid.cog0[2] << " to " << \
                                              grid.cog[0] << "," << grid.cog[1] << "," << grid.cog[2] << endl;

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
  Update the positions of the ligand and find the displacement of the center of geometry
  */
  ligand.atoms_jedi.clear();
  ligand.atomnames_jedi.clear();
  ligand.positions.clear();
  int protein_natoms=protein.atomnumbers.size(); // This will be the first index of the ligand;
  int total_natoms= protein_natoms+ligand.atomnumbers.size(); //This is the total number of requested atoms
  for (unsigned k=protein_natoms; k<total_natoms;k++)
  {
    ligand.positions.push_back(getPosition(k));
    
  }
  
  ligand.cog0=ligand.cog;
  ligand.compute_cog(ligand.positions);
  vector<double> displacement(3,0);
  displacement[0]=ligand.cog[0]-ligand.cog0[0];
  displacement[1]=ligand.cog[1]-ligand.cog0[1];
  displacement[2]=ligand.cog[2]-ligand.cog0[2];

  /*
   Move the grid by the exact displacement of the ligand center of geometry.
   IMPORTANT: This won't work if the ligand is not rigid
   IMPORTANT: The grid needs to be really BIG and symmetric (cubic, for example) in 
             order to minimise the effects of the lack of rotation
  */
  for (unsigned i=0; i<grid.positions.size();i++)
  {
    grid.positions[i][0]+=displacement[0];
    grid.positions[i][1]+=displacement[1];
    grid.positions[i][2]+=displacement[2];
  }
  
  /*
  Shrink grid and binding site so that only atoms/points that
  are within interacting distance from each other are kept
  */
  protein.atoms_jedi.clear();
  protein.atomnames_jedi.clear();
  protein.positions.clear();
  for (unsigned j=0; j<protein.atomnumbers.size();j++)
  {
    protein.positions.push_back(getPosition(j));
  }
  protein.select_atoms(protein.positions,grid.positions,protein.atomnames,params.r_max);
  
  
  //grid.select_atoms(grid.positions,protein.positions,grid.atomnames,params.r_max); //NOT YET
  
  
  /////////////////////////////////////////////////
  //                JEDI score                   //
  /////////////////////////////////////////////////
  
  distances distance_matrix;
  distance_matrix.compute_distance_matrix(protein.positions,grid.positions);

  
  mindist min_dist;
  min_dist.compute_mindist(distance_matrix.r_matrix,
                          distance_matrix.dr_matrix_dx,
                          distance_matrix.dr_matrix_dy,
                          distance_matrix.dr_matrix_dz,
                          params.theta);

                        
  contacts_Soff contacts;
  contacts.compute_contacts_S_off(distance_matrix.r_matrix,
                                  distance_matrix.dr_matrix_dx,
                                  distance_matrix.dr_matrix_dy,
                                  distance_matrix.dr_matrix_dz,
                                  params.r_hydro, params.deltar_hydro);

 
  contacts_sum contacts_sum;
  contacts_sum.compute_contacts_sum(contacts.contacts_matrix,
                                    contacts.d_contacts_dx,
                                    contacts.d_contacts_dy,
                                    contacts.d_contacts_dz,
                                    protein.atomnames_jedi);

  
  activity activity;
  activity.compute_activities(min_dist.min_dist, 
                              min_dist.d_mindist_dx, min_dist.d_mindist_dy, min_dist.d_mindist_dz,
                              params.CC_mind,params.deltaCC,
                              contacts_sum.contacts_total,
                              contacts_sum.d_contacts_total_dx,contacts_sum.d_contacts_total_dy,contacts_sum.d_contacts_total_dz,
                              params.Emin, params.deltaE);

  clustering clusters;
  clusters.cluster_grid(activity.activity_grid,grid.r_matrix,grid.neighbours,params.GP_max, params.resolution, activity.sum_activity);
  //clusters.print_clusters(grid.positions,activity.activity_grid,activity.S_on_mindist, activity.S_on_contacts);

  vector<unsigned> biggest_cluster=clusters.clusters[clusters.biggest_cluster_idx];
  //cout << "getting biggest cluster: " << clusters.biggest_cluster_idx << " with " <<biggest_cluster.size() << " elements" << endl;
  activity.filter_activities(biggest_cluster);
  contacts_sum.filter_contacts(biggest_cluster);
  
   
  Volume volume;
  double volume_element=pow(params.resolution,3);
  volume.compute_volume(activity.sum_activity,volume_element,
                        activity.d_sum_activity_dx,
                        activity.d_sum_activity_dy,
                        activity.d_sum_activity_dz);

  
  hydrophobicity hydrophobicity;
  hydrophobicity.compute_hydrophobicity(contacts_sum.contacts_apolar,
                                        contacts_sum.d_contacts_apolar_dx,contacts_sum.d_contacts_apolar_dy,contacts_sum.d_contacts_apolar_dz,
                                        contacts_sum.contacts_total,
                                        contacts_sum.d_contacts_total_dx,contacts_sum.d_contacts_total_dy,contacts_sum.d_contacts_total_dz,
                                        activity.activity_grid,
                                        activity.d_activity_dx,activity.d_activity_dy,activity.d_activity_dz,
                                        activity.sum_activity,
                                        activity.d_sum_activity_dx,activity.d_sum_activity_dy,activity.d_sum_activity_dz);
  

  double Jedi=params.alpha*volume.volume+params.beta*hydrophobicity.Ha;
  setValue(Jedi);

  vector<double> dJedi_dx(protein.atomnumbers.size(),0);
  vector<double> dJedi_dy(protein.atomnumbers.size(),0);
  vector<double> dJedi_dz(protein.atomnumbers.size(),0);

  #pragma omp parallel for
  for (unsigned j=0; j<protein.atoms_jedi.size();j++)
  {
    unsigned atom_idx=protein.atoms_jedi[j];
    dJedi_dx[atom_idx]=params.alpha*volume.d_volume_dx[j]+params.beta*hydrophobicity.d_Ha_dx[j];
    dJedi_dy[atom_idx]=params.alpha*volume.d_volume_dy[j]+params.beta*hydrophobicity.d_Ha_dy[j];
    dJedi_dz[atom_idx]=params.alpha*volume.d_volume_dz[j]+params.beta*hydrophobicity.d_Ha_dz[j];
  }
  
  #pragma omp parallel for
  for (unsigned j=0; j<protein.atomnumbers.size();j++)
  {
    setAtomsDerivatives(j,Vector(dJedi_dx[j],dJedi_dy[j],dJedi_dz[j]));
  }

  //PRINT OUTPUT FILE
  iteration++;
  int step=getStep();
  ofstream wfile;
  wfile.open("jedi_stats.dat",std::ios_base::app);
  wfile << step <<" " << activity.sum_activity << " " << volume.volume << " " << hydrophobicity.Ha << " " << Jedi << endl;
  wfile.close();
}

}
}