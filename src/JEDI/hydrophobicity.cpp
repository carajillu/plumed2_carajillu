#include <iostream>
#include "hydrophobicity.h"
#include <algorithm>
#include <functional>

using namespace std;

hydrophobicity::hydrophobicity()
{

};

void hydrophobicity::compute_hydrophobicity(vector<double> contacts_apolar,
                                            vector<vector<double>> d_apolar_dx,vector<vector<double>> d_apolar_dy,vector<vector<double>> d_apolar_dz,
                                            vector<double> total_contacts,
                                            vector<vector<double>> d_contacts_total_dx,vector<vector<double>> d_contacts_total_dy,vector<vector<double>> d_contacts_total_dz,
                                            vector<double> activity,
                                            vector<vector<double>> d_activity_dx,vector<vector<double>> d_activity_dy,vector<vector<double>> d_activity_dz,
                                            double sum_activity,
                                            vector<double> d_sum_activity_dx, vector<double> d_sum_activity_dy,vector<double> d_sum_activity_dz)
{
  unsigned size_grid=contacts_apolar.size();
  unsigned size_protein=d_apolar_dx[0].size();

  Ha=0;
  d_Ha_dx=vector<double>(size_protein,0);
  d_Ha_dy=vector<double>(size_protein,0);
  d_Ha_dz=vector<double>(size_protein,0);

  vector<double> d_hydroactivity_dx(size_protein,0);
  vector<double> d_hydroactivity_dy(size_protein,0);
  vector<double> d_hydroactivity_dz(size_protein,0);
  

  for (unsigned i=0; i<size_grid;i++)
  {
   double hydrophobicity_grid_i=contacts_apolar[i]/total_contacts[i];
   
   for(unsigned j=0;j<size_protein;j++)
   {
    double d_hydrogrid_dx_ij=(total_contacts[i]*d_apolar_dx[i][j]-d_contacts_total_dx[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
    double d_hydrogrid_dy_ij=(total_contacts[i]*d_apolar_dy[i][j]-d_contacts_total_dy[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
    double d_hydrogrid_dz_ij=(total_contacts[i]*d_apolar_dz[i][j]-d_contacts_total_dz[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
    d_hydroactivity_dx[j]+=activity[i]*d_hydrogrid_dx_ij+hydrophobicity_grid_i*d_activity_dx[i][j];
    d_hydroactivity_dy[j]+=activity[i]*d_hydrogrid_dy_ij+hydrophobicity_grid_i*d_activity_dy[i][j];
    d_hydroactivity_dz[j]+=activity[i]*d_hydrogrid_dz_ij+hydrophobicity_grid_i*d_activity_dz[i][j];
   }
   Ha+=activity[i]*hydrophobicity_grid_i; //note we are missing the denominator!
  }
  
  for (unsigned j=0;j<size_protein;j++)
  {
    d_Ha_dx[j]=(sum_activity*d_hydroactivity_dx[j]-Ha*d_sum_activity_dx[j])/pow(sum_activity,2);
    d_Ha_dy[j]=(sum_activity*d_hydroactivity_dy[j]-Ha*d_sum_activity_dy[j])/pow(sum_activity,2);
    d_Ha_dz[j]=(sum_activity*d_hydroactivity_dz[j]-Ha*d_sum_activity_dz[j])/pow(sum_activity,2);
  }
  Ha/=sum_activity; //Normalise Ha

 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<activity.activity.size();i++)
 {
     cout << "Point " << i << "Activity = " << activity.activity[i] << " h_i = " << hydrophobicity_grid[i] << endl;
 }

 cout << "Ha = " << Ha << endl;

 for (unsigned j=0; j<atomnames.size();j++)
 {
  cout << "Derivatives with respect to atom j " << j << ": " << d_Ha_dx[j] << " "<< d_Ha_dy[j] << " " << d_Ha_dz[j] << " " << endl;
 }
 */
}
