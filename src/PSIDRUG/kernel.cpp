#include <iostream>
#include "kernel.h"

using namespace std;

kernel::kernel(unsigned &n_atoms)
{
 r=vector<double>(n_atoms,0);
 dr_dx=vector<double>(n_atoms,0);
 dr_dy=vector<double>(n_atoms,0);
 dr_dz=vector<double>(n_atoms,0);

 activity=0;
 d_activity_dx=vector<double>(n_atoms,0);
 d_activity_dy=vector<double>(n_atoms,0);
 d_activity_dz=vector<double>(n_atoms,0);
}

void kernel::calculate_distance(vector<double> &gridpoint_crd, 
                                vector<PLMD::Vector> &atom_crd, 
                                vector<unsigned> &bsite_bin)
{
  for (unsigned j=0; j<atom_crd.size();j++)
  {
    if (bsite_bin[j]==0)
    {
        r[j]=0;
        dr_dx[j]=0;
        dr_dy[j]=0;
        dr_dz[j]=0;
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

void kernel::calculate_activity(vector<double> &gridpoint_crd, 
                                vector<PLMD::Vector> &atom_crd,
                                vector<unsigned> &bsite_bin)
{
  calculate_distance(gridpoint_crd,atom_crd,bsite_bin);
  for (unsigned j=0;j<atom_crd.size();j++)
  {
   activity+=r[j];
   d_activity_dx[j]+=dr_dx[j];
   d_activity_dy[j]+=dr_dy[j];
   d_activity_dz[j]+=dr_dz[j];
  }
}