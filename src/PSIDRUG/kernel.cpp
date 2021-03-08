#include <iostream>
#include "kernel.h"

using namespace std;

kernel::kernel(unsigned &n_atoms)
{
 r=vector<double>(n_atoms,0);
 dr_dx=vector<double>(n_atoms,0);
 dr_dy=vector<double>(n_atoms,0);
 dr_dz=vector<double>(n_atoms,0);

 exp_r=vector<double>(n_atoms,0);
 CC=0;
 dCC_dx=vector<double>(n_atoms,0);
 dCC_dy=vector<double>(n_atoms,0);
 dCC_dz=vector<double>(n_atoms,0);

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
 for (unsigned j=0;j<d_activity_dx.size();j++)
  {
   r[j]=0;
   dr_dx[j]=0;
   dr_dy[j]=0;
   dr_dz[j]=0;
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
 double r_min=999999999;

 sum_exp=0;
 for (unsigned j=0; j<bsite_bin.size();j++)
 {
   //bsite_bin[j]=0 means the atom is at least at rsite from the grid point (so r[j] was not calculated)
   //r[j]>HCmax means that exp_r has become flat at y=1
   if (bsite_bin[j]==0 or r[j]>HCmax)
   {
    exp_r[j]=1;
   }
   else
   {
    if (r[j]<r_min) r_min=r[j]; //for debug
    //exp_r[j]=exp(pow(((r[j]-x_offset)-HCmax),2))-1;
    exp_r[j]=exp(theta/r[j]); //for comparison with original
   }
   sum_exp+=exp_r[j];
 }
 CC=1/log(sum_exp); // we add theta at the end

 dCC_dr=0;
 for (unsigned j=0; j<bsite_bin.size();j++)
 {
   if (bsite_bin[j]==0 or r[j]>HCmax)
   {
    dCC_dx[j]=0;
    dCC_dy[j]=0;
    dCC_dz[j]=0;
   }
   else
   {
    dCC_dr=-pow(CC,2)/sum_exp*exp_r[j]*2*(r[j]-x_offset-HCmax);
    dCC_dx[j]=dCC_dr*dr_dx[j];
    dCC_dy[j]=dCC_dr*dr_dy[j];
    dCC_dz[j]=dCC_dr*dr_dz[j];
   }
 }
 CC*=theta;
 cout << "r_min = " << r_min << " mind = " << CC << endl;
}

void kernel::calculate_activity(vector<double> &gridpoint_crd, 
                                vector<PLMD::Vector> &atom_crd,
                                vector<unsigned> &bsite_bin)
{
  reset();
  calculate_distance(gridpoint_crd,atom_crd,bsite_bin);
  calculate_CC(bsite_bin);
  activity=CC;
  for (unsigned j=0;j<atom_crd.size();j++)
  {
   d_activity_dx[j]=dCC_dx[j];
   d_activity_dy[j]=dCC_dy[j];
   d_activity_dz[j]=dCC_dz[j];
  }
  
  
}