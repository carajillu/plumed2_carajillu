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
#include "core/ActionRegister.h"
#include "tools/PDB.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <armadillo>
#include <omp.h>
#include "aidefunctions.h"

// CV modules
#include "probe.h"

using namespace std;
using namespace std::chrono;

namespace PLMD
{
  namespace colvar
  {

    /*+PLUMEDOC COLVAR TEMPLATE
    Add CV info
    +ENDPLUMEDOC*/

    class Ghostprobe : public Colvar
    {
      // Execution control variables
      int nthreads=0;     // number of available OMP threads
      int ndev=0;         // number of available OMP accelerators
      //calculation speed
      bool performance;
      time_point<high_resolution_clock> start_psi, end_psi;
      time_point<high_resolution_clock> start_dxfix, end_dxfix;

      // All of these are just for remove_netforcetorque()
      time_point<high_resolution_clock> start_tor, end_tor;
      time_point<high_resolution_clock> start_A, end_A;
      time_point<high_resolution_clock> start_B, end_B;
      time_point<high_resolution_clock> start_Bt, end_Bt;
      time_point<high_resolution_clock> start_c, end_c;
      time_point<high_resolution_clock> start_correction, end_correction;
      time_point<high_resolution_clock> start_test, end_test;

      // MD control variables
      bool pbc;
      // CV control variables
      bool nocvcalc;
      bool nodxfix;
      bool noupdate;
      bool botch_derivatives;
      double kpert=0;
      double kxplor=0;
      unsigned pertstride=0;
      string ref_lig;
      PDB lig_pdb;
      PDB atoms_pdb; //pdb with the involved atoms and their hydrophobicity coefficients as 
      bool restart_probes;
      // Parameters
      double Rmin=0, deltaRmin=0;
      double Rmax=0, deltaRmax=0;
      double Pmin=0, deltaP=0;
      double Cmin=0, deltaC=0;
      double Hmin=0, deltaH=0;
      vector<double> h_coeff;

      // Set up of CV
      vector<PLMD::AtomNumber> atoms; // indices of atoms supplied to the CV (starts at 1)
      unsigned n_atoms=0;               // number of atoms supplied to the CV
      vector<double> atoms_x, atoms_y, atoms_z;
      
      unsigned step=0;

      vector<PLMD::AtomNumber> atoms_init; // Indices of the atoms in which the probes will be initially centered
      unsigned n_init=0;                     // number of atoms used in ATOMS_INIT
      vector<unsigned> init_j;             // Indices of atoms_init in getPositions()

      vector<Probe> probes; // This will contain all the spherical probes
      unsigned nprobes=0;     // number of spherical probes to use

      // Output control variables
      unsigned probestride=0; // stride to print information for post-processing the probe coordinates

      // Calculation of CV and its derivatives

      double Psi=0;
      vector<double> d_Psi_dx, d_Psi_dy, d_Psi_dz;

      // Correction of derivatives
      unsigned torquestride=0;
      double kappa=0;
      vector<double> force_x, force_y, force_z;
      vector<double> torque_x, torque_y, torque_z;
      // we need to solve c=At*inv(A*At)*L
      arma::mat A;  // dim=6,3*n_atoms
      arma::mat At; // dim=3*n_atoms,6
      arma::vec L;  // dim=6
      arma::vec c;  // 3*n_atoms

      double err_tol=1e-8;
      unsigned null_fails=0;
      bool dumpderivatives;

    public:
      explicit Ghostprobe(const ActionOptions &);
      // active methods:
      void calculate() override;
      void reset();
      void remove_netforcetorque();
      void print_protein();
      void get_init_crd();
      static void registerKeywords(Keywords &keys);
    };

    PLUMED_REGISTER_ACTION(Ghostprobe, "GHOSTPROBE")

    void Ghostprobe::registerKeywords(Keywords &keys)
    {
      Colvar::registerKeywords(keys);
      keys.addFlag("DEBUG", false, "Running in debug mode");
      keys.addFlag("NOCVCALC", false, "skip CV calculation");
      keys.addFlag("NOUPDATE", false, "skip probe update");
      keys.addFlag("NODXFIX", false, "skip derivative correction");
      keys.addFlag("PERFORMANCE", false, "measure execution time");
      keys.addFlag("DUMPDERIVATIVES", false, "print derivatives and corrections");
      keys.addFlag("RESTART_PROBES", false, "Restart probe positions from stored coordinates");
      keys.addFlag("BOTCH_DERIVATIVES", false, "Botch the derivatives to separate them when using HARMONIC potentials ONLY");
      keys.add("atoms", "ATOMS", "Atoms to include in druggability calculations (start at 1)");
      keys.add("atoms", "ATOMS_INIT", "Atoms in which the probes will be initially centered.");
      keys.add("optional","PDB","PDB file of the atoms involved, for hydrophobicity calculations");
      keys.add("optional", "NPROBES", "Number of probes to use");
      keys.add("optional", "PROBESTRIDE", "Print probe coordinates info every PROBESTRIDE steps");
      keys.add("optional", "RMIN", "");
      keys.add("optional", "DELTARMIN", "");
      keys.add("optional", "RMAX", "");
      keys.add("optional", "DELTARMAX", "");
      keys.add("optional", "CMIN", "");
      keys.add("optional", "DELTAC", "");
      keys.add("optional", "PMIN", "");
      keys.add("optional", "DELTAP", "");
      keys.add("optional", "HMIN", "");
      keys.add("optional", "DELTAH", "");
      keys.add("optional", "KPERT", "");
      keys.add("optional", "KXPLOR", "");
      keys.add("optional", "PERTSTRIDE", "Do a full KPERT random perturbation every PERTSTRIDE steps");
      keys.add("optional","REF_LIG","Coordinates od reference ligand atoms to place the probes on.");
      keys.add("optional", "KAPPA", "");
      keys.add("optional", "TORQUESTRIDE", "stride for removal of accummulated forces and torques");
    }

    Ghostprobe::Ghostprobe(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao),
                                                pbc(true),
                                                nocvcalc(false),
                                                noupdate(false),
                                                nodxfix(false),
                                                performance(false),
                                                dumpderivatives(false),
                                                restart_probes(false)
    {
/*
Initialising openMP threads.
This does not seem to be affected by the environment variable $PLUMED_NUM_THREADS
*/    #pragma omp parallel 
      {
      nthreads = omp_get_num_threads();
      ndev = omp_get_num_devices();
      }
      cout << "------------ Available Computing Resources -------------" << endl;
      cout << "Ghostprobe initialised with " << nthreads << " OMP threads " << endl;
      cout << "and " << ndev << " OMP compatible accelerators (not currently used)" << endl;

      addValueWithDerivatives();
      setNotPeriodic();

      bool nopbc = !pbc;
      parseFlag("NOPBC", nopbc);
      pbc = !nopbc;

      parseFlag("NOCVCALC",nocvcalc);
      parseFlag("NOUPDATE", noupdate);
      parseFlag("NODXFIX", nodxfix);
      parseFlag("PERFORMANCE", performance);
      if (performance)
      {
       ofstream wfile;
       wfile.open("performance.txt");
       wfile << "Psi correction" << endl;
       wfile.close();
       
       //ofstream wfile;
       wfile.open("performance_dxfix.txt");
       wfile << "torques A B Bcoot Bt c correction test total" << endl;
       wfile.close();
      }

      parseFlag("BOTCH_DERIVATIVES", botch_derivatives);
      if (botch_derivatives)
      {
       cout << "*****************************************************************************************" << endl;
       cout << "WARNING:Botching derivatives. Use ONLY with RESTRAINT and KAPPA as the biasing method!!!!" << endl;
       cout << "*****************************************************************************************" << endl;
      }

      parseFlag("DUMPDERIVATIVES",dumpderivatives);
      if (dumpderivatives)
      {
        ofstream wfile;
        wfile.open("derivatives.csv");
        wfile << "Step Atom dx dy dz cx cy cz" << endl;
        wfile.close();

        wfile.open("forces_torques.csv");
        wfile << "Step Atom fx fy fz tx ty tz fcx fcy fcz" << endl;
        wfile.close();
      }

      parseFlag("RESTART_PROBES",restart_probes);

      parseAtomList("ATOMS", atoms);
      n_atoms = atoms.size();
      

      parse("NPROBES", nprobes);
      if (!nprobes)
      {
        nprobes = 16;
      }
      
      /*
      The following bit checks if ATOMS_INIT has been specified.
      If so, all probes will be initialised in the centre of ATOMS_INIT.
      Otherwise, each probe will be initialised on a protein atom chosen at random.
      */
      parseAtomList("ATOMS_INIT", atoms_init);
      n_init = atoms_init.size();

      parse("REF_LIG",ref_lig);
      
      if (!ref_lig.empty())
      {
        if( !lig_pdb.read(ref_lig,usingNaturalUnits(),0.1/getUnits().getLength()) )
         {
          cout << "missing input file " << ref_lig << endl;
          exit(1);
         }
        unsigned lig_natoms=lig_pdb.getAtomNumbers().size();
        cout << ref_lig << " has " << lig_natoms << " atoms. Changing NPROBES from " << nprobes << " to " << lig_natoms << endl;
        nprobes=lig_natoms;
      }
      else if (n_init==0 and !restart_probes)
      {
        cout << "geting random atoms_init" << endl;
        for (unsigned i=0; i<nprobes; i++)
        {
          unsigned j_rand=aidefunctions::get_random_integer(0, n_atoms-1);
          init_j.push_back(j_rand);
          cout << j_rand << endl;
        }
        cout << "got random atoms_init" << endl;
      }
      else
      {
        for (unsigned j=0;j<n_init; j++)
        {
          int j_idx=aidefunctions::findIndex(atoms,atoms_init[j]);
          if (j_idx==-1)
          {
            cout << "Atom " << atoms_init[j].index() << " not found in ATOMS. Adding it." << endl;
            atoms.push_back(atoms_init[j]);
            init_j.push_back(atoms.size()-1);
          }
          else
          {
            init_j.push_back(j_idx);
          }
        }
      }

      cout << "Requesting " << atoms.size() << " atoms" << endl;
      requestAtoms(atoms);
      cout << "--------- Initialising Ghostprobe Collective Variable -----------" << endl;

      //Get PDB file for hydrophobicity
      string pdb_protein;
      parse("PDB",pdb_protein);
      if (pdb_protein.empty())
      {
        cout << "PDB not specified. Hydrophobicity coefficients will all be set to 1." << endl;
        h_coeff=vector<double> (n_atoms, 1.0);
      }
      else
      {
        /* This is using push_back(), so the coefficients are in the right order
         even if the atom numbers in the pdb file are not correct (trjconv restarts them and puts them in order)*/
        atoms_pdb.read(pdb_protein, usingNaturalUnits(), 0.1/getUnits().getLength());
        const std::vector<AtomNumber>& atom_numbers = atoms_pdb.getAtomNumbers();
        for (const auto& atom_num : atom_numbers)
        {
          std::string atom_name = atoms_pdb.getAtomName(atom_num);
          if (!atom_name.empty()) {
            char first = atom_name[0];
            if (first == 'C' || first == 'S') {
              h_coeff.push_back(1.0); // Hydrophobicity coefficient for C and S atoms
            } else {
              h_coeff.push_back(0.0); // Default hydrophobicity coefficient for other atoms
            }
          }
        }
      }

      cout << "Using " << nprobes << " spherical probe(s) with the following parameters:" << endl;

      parse("RMIN", Rmin);
      if (!Rmin)
        Rmin = 0.;
      cout << "Rmin = " << Rmin << " nm" << endl;

      parse("DELTARMIN", deltaRmin);
      if (!deltaRmin)
        deltaRmin = 0.15;
      cout << "deltaRmin = " << deltaRmin << " nm" << endl;

      parse("RMAX", Rmax);
      if (!Rmax)
        Rmax = 0.6;
      cout << "Rmax = " << Rmax << " nm" << endl;

      parse("DELTARMAX", deltaRmax);
      if (!deltaRmax)
        deltaRmax = 0.15;
      cout << "deltaRmax = " << deltaRmax << " nm" << endl;

      parse("CMIN", Cmin);
      if (!Cmin)
        Cmin = 0.; 
      cout << "CMIN = " << Cmin << endl;

      parse("DELTAC", deltaC); 
      if (!deltaC)
        deltaC = 5; 
      cout << "DELTAC = " << deltaC << endl;

      parse("PMIN", Pmin);
      if (!Pmin)
        Pmin = 0; 
      cout << "PMIN = " << Pmin << endl;

      parse("DELTAP", deltaP);
      if (!deltaP)
        deltaP = 17; 
      cout << "DELTAP = " << deltaP << endl;

      parse("HMIN", Hmin);
      if (!Hmin)
        Hmin = 0; 
      cout << "HMIN = " << Hmin << endl;

      parse("DELTAH", deltaH);
      if (!deltaH)
        deltaH = 17; 
      cout << "DELTAH = " << deltaH << endl;
       
      /*
      cout << "Hydrophobicity coefficients" << endl;
      const std::vector<AtomNumber>& atom_numbers = atoms_pdb.getAtomNumbers();
      for (unsigned j=0; j<n_atoms; j++)
      {
        cout << atoms_pdb.getAtomName(atom_numbers[j]) << ": " << h_coeff[j] << endl;
      }
      */
      

      parse("KPERT",kpert);
      if (!kpert)
      {
        cout << "****************************************************************************" << endl;
        cout << "KPERT HAS either not been set, or manually set to zero." << endl;
        cout << "WARNING: PROBE WILL NOT BE PERTURBED AND POCKET SEARCH WILL NOT BE PERFORMED" << endl;
        cout << "****************************************************************************" << endl;
      }
      else
      {
      cout << "Perturbations of " << kpert << " nm will be applied to all probes." << endl;
      cout << "Perturbations will go in the direction opposite to the derivatives of the activity" << endl;
      cout << "with respect to the probe, when possible. (Fx=-dV/dx)" << endl;
      }
      
      parse("PERTSTRIDE",pertstride);
      parse("KXPLOR",kxplor);
      if (kxplor)
      {
      cout << "Random perturbations of " << kxplor << " nm will be applied to all probes every" << pertstride << " steps" << endl;
      }

      for (unsigned i = 0; i < nprobes; i++)
      {
        probes.push_back(Probe(i, restart_probes,
                               Rmin, deltaRmin, 
                               Rmax, deltaRmax, 
                               Cmin, deltaC, 
                               Pmin, deltaP,
                               Hmin, deltaH, h_coeff,
                               kpert, kxplor,pertstride,
                               n_atoms));
        cout << "Probe " << i << " initialised" << endl;
      }

      // parameters used to control output
      parse("PROBESTRIDE", probestride);
      if (!probestride)
        probestride = 1;
      cout << "Information to post-process probe coordinates will be printed every " << probestride << " steps" << endl
           << endl;

      // Allocate space for atom coordinates

      atoms_x = vector<double>(n_atoms, 0);
      atoms_y = vector<double>(n_atoms, 0);
      atoms_z = vector<double>(n_atoms, 0);

      cout << "---------Initialisng Ghostprobe and its derivatives---------" << endl;
      Psi = 0;
      d_Psi_dx = vector<double>(n_atoms, 0);
      d_Psi_dy = vector<double>(n_atoms, 0);
      d_Psi_dz = vector<double>(n_atoms, 0);

      if (nocvcalc)
         cout << "WARNING: NOCVCALC flag has been included. CV will NOT be calculated." << endl;

      if (!nodxfix)
      {
        cout << "---------Initialisng correction of Ghostprobe derivatives---------" << endl;
        parse("KAPPA",kappa);
        cout << "KAPPA = " << kappa << " kJ/mol. Please make sure this is the same KAPPA you are using in RESTRAINT" << endl;
        if (!kappa)
        {
          throw std::invalid_argument("Net force and toorque correction has been requested, but the harmonic force constant has not been provided. \
            Please provide the same KAPPA you provided in RESTRAINT (and if you are not using RESTRAINT with KAPPA, this won't work)");
        }
        parse("TORQUESTRIDE",torquestride);
        if (!torquestride)
        {
          torquestride=1;
          cout << "WARNING: Removing net forces and torques at every step. This will be slow." << endl;
        }

        force_x=vector<double>(n_atoms,0);
        force_y=vector<double>(n_atoms,0);
        force_z=vector<double>(n_atoms,0);
        torque_x=vector<double>(n_atoms,0);
        torque_y=vector<double>(n_atoms,0);
        torque_z=vector<double>(n_atoms,0);
        // Accessing individual elements tends to be slow, so we set most of the matrix here (see overleaf Thesis 2020 for equations)
        // Most elements are always 0 and 1, the only ones that change are set in remove_netforcetorque()
        A=arma::mat(6,3*n_atoms,arma::fill::zeros);
        for (unsigned j=0; j<n_atoms;j++)
        {
         A.row(0).col(j+ 0*n_atoms) = 1.0; //cx coefficients
         A.row(1).col(j+ 1*n_atoms) = 1.0; //cy coefficients
         A.row(2).col(j+ 2*n_atoms) = 1.0; //cz coefficients
        }
        At=arma::mat(3*n_atoms,6);
        L=arma::vec(6,arma::fill::zeros); 
        c=arma::vec(3*n_atoms,arma::fill::zeros);


      }
      else
      {
        cout << "Ghostprobe derivatives are not going to be corrected" << endl;
        cout << "Use the NODXFIX flag with care, as this means that" << endl;
        cout << "the sum of forces in the system will not be zero" << endl;
      }

      cout << "--------- Initialisation complete -----------" << endl;
      checkRead();
    }

    // reset Ghostprobe and derivatives to 0
    void Ghostprobe::reset()
    {
      Psi = 0;
      fill(d_Psi_dx.begin(), d_Psi_dx.end(), 0);
      fill(d_Psi_dy.begin(), d_Psi_dy.end(), 0);
      fill(d_Psi_dz.begin(), d_Psi_dz.end(), 0);
    }

    void Ghostprobe::remove_netforcetorque()
    {
     /*
     1) Calculate net force on each atom per probe (needs KAPPA and derivative with respect to atom/coordinate per probe)
     2) Acummulate over atoms/coordinates
     3) Print PLUMED forces and compare (needs botched derivatives)
     4) Accummulate forces and torques over time
     5) if not time to fix: return
     6) else: calculate big correction and reset A to all zeroes
     */
    //1)
    
    #pragma omp parallel for
    for (unsigned j=0; j<n_atoms; j++)
    {
     // individual forces. 
     // these are calculated the way PLUMED calculates them (bacause those are the forces/torques we want to remove)
     // this is independent of wether or not we botch, but only works if using a harmonic potential
     double fxij=-kappa*(Psi-1)*d_Psi_dx[j];
     double fyij=-kappa*(Psi-1)*d_Psi_dy[j];
     double fzij=-kappa*(Psi-1)*d_Psi_dz[j];
     force_x[j]+=fxij;
     force_y[j]+=fyij;
     force_z[j]+=fzij;
     torque_x[j]+=atoms_y[j]*fzij-atoms_z[j]*fyij;
     torque_y[j]+=atoms_z[j]*fxij-atoms_x[j]*fzij;
     torque_z[j]+=atoms_x[j]*fyij-atoms_y[j]*fxij;
    }

    // Build matrix A (see equations on Thesis2020 overleaf)
    for (unsigned j=0; j<n_atoms;j++)
    {
     //cx coefficients
     A.row(4).col(j+ 0*n_atoms) = atoms_z[j];
     A.row(5).col(j+ 0*n_atoms) = -atoms_y[j];
     //cy coefficients 
     A.row(3).col(j+ 1*n_atoms) = -atoms_z[j];
     A.row(5).col(j+ 1*n_atoms) = atoms_x[j];
     //cz coefficients
     A.row(3).col(j+ 2*n_atoms) = atoms_y[j];
     A.row(4).col(j+ 2*n_atoms) = -atoms_x[j];
    }
    // Transpose A
    At=arma::trans(A);
    // Build vector L
    L.fill(0);
    for (unsigned j=0; j<n_atoms; j++)
    {
      L[0]-=force_x[j];
      L[1]-=force_y[j];
      L[2]-=force_z[j];
      L[3]-=torque_x[j];
      L[4]-=torque_y[j];
      L[5]-=torque_z[j];
    }
    L/=-kappa*(Psi-1);
    //cout << "Step " << step << ": L, before correction: " << L[0] << " " << L[1] << " " << L[2] << " " << L[3] << " " << L[4] << " " << L[5] << endl;

    //get constants
    arma::vec c = At*arma::pinv(A*At)*L;
    
    // apply correction
    for (unsigned j = 0; j < n_atoms; j++)
    {
     //cout << c[j + 0 * n_atoms] << " "<< c[j + 1 * n_atoms] << " "<< c[j + 2 * n_atoms] << endl;
     d_Psi_dx[j] += c[j + 0 * n_atoms];
     d_Psi_dy[j] += c[j + 1 * n_atoms];
     d_Psi_dz[j] += c[j + 2 * n_atoms];
    }

    // check results
    L.fill(0);
    for (unsigned j=0; j<n_atoms; j++)
    {
     double fxij=-kappa*(Psi-1)*d_Psi_dx[j];
     double fyij=-kappa*(Psi-1)*d_Psi_dy[j];
     double fzij=-kappa*(Psi-1)*d_Psi_dz[j];
     L[0]+=fxij;
     L[1]+=fyij;
     L[2]+=fzij;
     L[3]+=atoms_y[j]*fzij-atoms_z[j]*fyij;
     L[4]+=atoms_z[j]*fxij-atoms_x[j]*fzij;
     L[5]+=atoms_x[j]*fyij-atoms_y[j]*fxij;
    }
    if (L[0]>1e-8 or L[1]>1e-8 or L[2]>1e-8 or
        L[3]>1e-8 or L[4]>1e-8 or L[5]>1e-8)
    {
      cout << "Error: Removal of net forces and torques failed. Simulation will now end." << endl;
      cout << "Sum forces: "  << L[0] << " " << L[1] << " " << L[2] << endl;
      cout << "Sum torques: " << L[3] << " " << L[4] << " " << L[5] << endl;
      exit(0);
    }

    if (dumpderivatives and step%probestride==0)
    {
     ofstream wfile;
     wfile.open("forces_torques.csv",std::ios_base::app);
     for (unsigned j=0; j<n_atoms; j++)
     {
       wfile << step << " " << j << " " 
             << force_x[j] << " " << force_y[j] << " " << force_z[j] << " " 
             << torque_x[j] << " " << torque_y[j] << " " << torque_z[j] << " "
             << -kappa*(Psi-1)*c[j + 0 * n_atoms] << " " << -kappa*(Psi-1)*c[j + 1 * n_atoms] << " "<< -kappa*(Psi-1)*c[j + 2 * n_atoms] << " "
             << endl;
     }
     wfile.close();

     wfile.open("derivatives.csv",std::ios_base::app);
     for (unsigned j=0; j<n_atoms; j++)
     {
       wfile << step << " " << j << " " 
             << d_Psi_dx[j] << " " << d_Psi_dy[j] << " " << d_Psi_dz[j] << " "
             << c[j + 0 * n_atoms] << " " << c[j + 1 * n_atoms] << " "<< c[j + 2 * n_atoms] << " "
             << endl;
     }
     wfile.close();
    }
     
    //cout << "Removed net forces and torques. New L: " << L[0] << " " << L[1] << " " << L[2] << " " << L[3] << " " << L[4] << " " << L[5] << endl;
    fill(force_x.begin(), force_x.end(), 0.0);
    fill(force_y.begin(), force_y.end(), 0.0);
    fill(force_z.begin(), force_z.end(), 0.0);
    fill(torque_x.begin(), torque_x.end(), 0.0);
    fill(torque_y.begin(), torque_y.end(), 0.0);
    fill(torque_z.begin(), torque_z.end(), 0.0);

    return;
    }

    void Ghostprobe::print_protein()
    {
     string filename = "protein.xyz";
     ofstream wfile;
     if (step==0)
     {
      wfile.open(filename.c_str());
     }
     else
     {
      wfile.open(filename.c_str(),std::ios_base::app);
     }
     wfile << n_atoms << endl;
     wfile << "Step  "<< to_string(step) << endl;
     for (unsigned j=0; j<n_atoms;j++)
     {
      wfile << atoms[j].serial() << " " << std::fixed << std::setprecision(5) << atoms_x[j]*10 << " " << atoms_y[j]*10 << " " << atoms_z[j]*10 << endl;  
     } 
     wfile.close();
    }

    void Ghostprobe::get_init_crd()
    {
      double x=0;
      double y=0;
      double z=0;

      if (restart_probes) //restart probes from input coordinates
      {
       vector<vector<double>> protein_xyz=aidefunctions::read_xyz("protein.xyz",0);
       for (unsigned i=0; i<nprobes; i++)
       {
        cout << "Restarting probe " << i << endl;
        string filename = "probe-";
        filename.append(to_string(i));
        filename.append(".xyz");
        vector<vector<double>> probe_xyz=aidefunctions::read_xyz(filename,0);
        x=probe_xyz[0][0];
        y=probe_xyz[0][1];
        z=probe_xyz[0][2];
        probes[i].place_probe(x,y,z);
        probes[i].get_atoms_restart(protein_xyz);
       }
      }
      else if (!ref_lig.empty()) // place probes on the coordinates of an input ligand
      {
        for (unsigned i=0; i<nprobes; i++)
        {
          x=lig_pdb.getPositions()[i][0];
          y=lig_pdb.getPositions()[i][1];
          z=lig_pdb.getPositions()[i][2];
          probes[i].place_probe(x,y,z);
        }  
      }
      else if (atoms_init.size()!=0) // Place ALL probes in the geometric centre of a set of atoms
      {
       for (unsigned j=0; j<init_j.size();j++)
       {
        //cout << j << " " << init_j[j] << " " << getPosition(init_j[j])[0]<< " " << getPosition(init_j[j])[1]<< " " << getPosition(init_j[j])[2] <<  endl;
        x+=getPosition(init_j[j])[0]/init_j.size();
        y+=getPosition(init_j[j])[1]/init_j.size();
        z+=getPosition(init_j[j])[2]/init_j.size();
       }
       for (unsigned i = 0; i < nprobes; i++)
       {
        probes[i].place_probe(x,y,z);
       }
       cout << "All probes are initialised at point " << x << " " << y << " " << z << endl;
      }
      else 
      {
        for (unsigned i = 0; i < nprobes; i++)
        {
          x=getPosition(init_j[i])[0];
          y=getPosition(init_j[i])[1];
          z=getPosition(init_j[i])[2];
          probes[i].place_probe(x,y,z);
          probes[i].perturb_probe(0);
          cout << "Probe " << i << " centered on atom " << atoms[init_j[i]].serial() << endl;
        }
      }
      return;
    }
    // calculator
    void Ghostprobe::calculate()
    {
      if (pbc)
        makeWhole();
      reset();

      step = getStep();
      
      // Get atom positions
      for (unsigned j = 0; j < n_atoms; j++)
      {
        atoms_x[j] = getPosition(j)[0];
        atoms_y[j] = getPosition(j)[1];
        atoms_z[j] = getPosition(j)[2];
      }

      // At step 0, place the probes using get_init_crd()
      if (step==0 or noupdate)
         get_init_crd();
      
      if (performance and step%probestride==0) start_psi = high_resolution_clock::now();
      
      /////////////////////////////////////////////
      // PSI score
      //////////////////////////////////////////////
      #pragma omp parallel for 
      for (unsigned i = 0; i < nprobes; i++)
      {
        // Update probe coordinates
        if (!noupdate)
        {
         probes[i].move_probe(step, atoms_x, atoms_y, atoms_z);
        }

        //Calculate Psi and its derivatives
        if (!nocvcalc)
        {
          probes[i].calculate_activity(atoms_x, atoms_y, atoms_z);
          #pragma omp critical //avoid race condition
          {
           Psi+=probes[i].activity/nprobes;
          }
        }

        /*
        The following needs to go IN THIS ORDER: print_probe_xyz(), perturb_probe(), print_prove_movement()
        because print_probe_movement() records the type of perturbation.
        */
        if (step%probestride==0) probes[i].print_probe_xyz(step);
        if (kpert>0) probes[i].perturb_probe(step);
        if (step%probestride==0) probes[i].print_probe_movement(step,atoms,n_atoms);
      }
      /////////////////////////////////////////////////
      //DERIVATIVES
      ////////////////////////////////////////////////
      #pragma omp parallel for
      for (unsigned i=0; i<nprobes; i++)
      {
       if (botch_derivatives) probes[i].botch_derivatives(Psi);
       #pragma omp critical // avoid race condition
       {
         for (unsigned j=0;j<n_atoms;j++)
         {
          d_Psi_dx[j]+=probes[i].d_activity_dx[j]/nprobes;
          d_Psi_dy[j]+=probes[i].d_activity_dy[j]/nprobes;
          d_Psi_dz[j]+=probes[i].d_activity_dz[j]/nprobes;
         }
       }
      }
      
      if (performance and step%probestride==0)  end_psi = high_resolution_clock::now();
      
      //Correct the Psi derivatives so that they sum 0
      if (!nodxfix and step%torquestride==0)
      {
       remove_netforcetorque();
      }
   
      //Send Psi and derivatives back to Plumed
      setValue(Psi);
      for (unsigned j=0;j<n_atoms;j++)
      {
        setAtomsDerivatives(j,Vector(d_Psi_dx[j],d_Psi_dy[j],d_Psi_dz[j]));
      }

      //print output for post_processing
       if (step % probestride == 0)
       {
         print_protein();
       }
      
      if (performance and step%probestride==0)
      {
       int psi_time = duration_cast<microseconds>(end_psi - start_psi).count();
       int dxfix_time = duration_cast<microseconds>(end_dxfix - start_dxfix).count();
       ofstream wfile;
       wfile.open("performance.txt",std::ios_base::app);
       wfile << psi_time << " " << dxfix_time << endl;
       wfile.close();
      }

      // if (step>=10) exit(0);

    } // close calculate
  }   // close colvar
} // close plmd
