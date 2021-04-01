#include <iostream>
#include "kernel.h"

using namespace std;

kernel::kernel(unsigned &n_atoms, 
               double ccmin, double deltacc, 
               double srmin, double deltasr)
{
 //parameters
 CCmin=ccmin;
 deltaCC=deltacc;
 SRmin=srmin;
 deltaSR=deltasr;
 
 //distance
 r=vector<double>(n_atoms,0);
 dr_dx=vector<double>(n_atoms,0);
 dr_dy=vector<double>(n_atoms,0);
 dr_dz=vector<double>(n_atoms,0);

 //close contact
 CC=0;
 dCC_dx=vector<double>(n_atoms,0);
 dCC_dy=vector<double>(n_atoms,0);
 dCC_dz=vector<double>(n_atoms,0);

 //hydrophobicity
 hydrophobicity=0;
 d_hydrophobicity_dx=vector<double>(n_atoms,0);
 d_hydrophobicity_dy=vector<double>(n_atoms,0);
 d_hydrophobicity_dz=vector<double>(n_atoms,0);

 activity=0;
 d_activity_dx=vector<double>(n_atoms,0);
 d_activity_dy=vector<double>(n_atoms,0);
 d_activity_dz=vector<double>(n_atoms,0);
}

// Replace all the components of the vectors with zeros so that
// the same kernel can be used for another grid point
void kernel::reset()
{
 activity=0;
 CC=0;
 hydrophobicity=0;
 for (unsigned j=0;j<d_activity_dx.size();j++)
  {
   //reset distance
   r[j]=0;
   dr_dx[j]=0;
   dr_dy[j]=0;
   dr_dz[j]=0;
   //reset close contact
   dCC_dx[j]=0;
   dCC_dy[j]=0;
   dCC_dz[j]=0;
   //reset hydrophobicity
   d_hydrophobicity_dx[j]=0;
   d_hydrophobicity_dy[j]=0;
   d_hydrophobicity_dz[j]=0;
   //reset activity
   d_activity_dx[j]=0;
   d_activity_dy[j]=0;
   d_activity_dz[j]=0;
  }
}

void kernel::calculate_distance(vector<double> &gridpoint_crd, 
                                vector<PLMD::Vector> &atom_crd, 
                                vector<unsigned> &bsite_bin)
{
  for (unsigned j=0; j<atom_crd.size();j++)
  {
    if (bsite_bin[j]==0)
    {   // values are already set to zero at kernel::reset()
        continue;
    }
    r[j]=sqrt(pow((atom_crd[j][0]-gridpoint_crd[0]),2)+
              pow((atom_crd[j][1]-gridpoint_crd[1]),2)+
              pow((atom_crd[j][2]-gridpoint_crd[2]),2));
    
    dr_dx[j]=(atom_crd[j][0]-gridpoint_crd[0])/r[j];
    dr_dy[j]=(atom_crd[j][1]-gridpoint_crd[1])/r[j];
    dr_dz[j]=(atom_crd[j][2]-gridpoint_crd[2])/r[j];
  }
}

void kernel::calculate_CC(vector<unsigned> &bsite_bin)
{
  CC=1;
  double m=0;
  double dSoff_j=0;
  for (unsigned j=0; j<bsite_bin.size();j++)
  {
    if (bsite_bin[j]==0)
    {
      continue;
    }
    else
    {
     m=m_v(r[j],CCmin,deltaCC);
     CC+=Soff_m(m,kCCmin);
     //derivatives are still missing the total CC bit!
     dSoff_j=dSoff_dm(m,kCCmin)*dm_dv(deltaCC);
     dCC_dx[j]=dSoff_j*dr_dx[j];
     dCC_dy[j]=dSoff_j*dr_dy[j];
     dCC_dz[j]=dSoff_j*dr_dz[j];
    }
  }
  CC=1/CC;
  
  for (unsigned j=0; j<bsite_bin.size();j++)
  { 
    if (bsite_bin[j]==0)
       continue;
    else
    {
     dCC_dx[j]*=-pow(CC,2);
     dCC_dy[j]*=-pow(CC,2);
     dCC_dz[j]*=-pow(CC,2);
    }
  }
}

void kernel::calculate_hydrophobicity(vector<unsigned> &bsite_bin, vector<double> &charges)
{
  
}

void kernel::calculate_activity(vector<double> &gridpoint_crd, 
                                vector<PLMD::Vector> &atom_crd,
                                vector<unsigned> &bsite_bin,
                                vector <double> &charges)
{
  reset();
  calculate_distance(gridpoint_crd,atom_crd,bsite_bin);
  calculate_CC(bsite_bin);
  calculate_hydrophobicity(bsite_bin,charges);
  activity=CC*hydrophobicity;
  for (unsigned j=0;j<atom_crd.size();j++)
  {
   d_activity_dx[j]=hydrophobicity*dCC_dx[j]+CC*d_hydrophobicity_dx[j];
   d_activity_dy[j]=hydrophobicity*dCC_dy[j]+CC*d_hydrophobicity_dy[j];
   d_activity_dz[j]=hydrophobicity*dCC_dz[j]+CC*d_hydrophobicity_dz[j];
  }
  
  
}