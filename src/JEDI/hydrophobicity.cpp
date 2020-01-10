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
                                            vector<double> contacts_polar,
                                            vector<vector<double>> d_polar_dx,vector<vector<double>> d_polar_dy,vector<vector<double>> d_polar_dz,
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
  vector<double> d_hydroactivity_dz(size_protein,0);
  vector<double> d_hydroactivity_dy(size_protein,0);

  for (unsigned i=0; i<size_grid;i++)
  {
   double total_contacts=contacts_apolar[i]+contacts_polar[i];
   /* Notice here that d_polar values are added to d_apolar. The trick here is that the derivative
      of a polar contact with respect to an apolar atom will always be 0, and viceversa, so we can effectively
      use the same vector for both polar and apolar.
      This is using a local copy of the object so it should not affect the values of the 
      original contacts object
   */
   transform(d_polar_dx[i].begin(),d_polar_dx[i].end(),d_apolar_dx[i].begin(),d_polar_dx[i].begin(),plus<double>());
   transform(d_polar_dy[i].begin(),d_polar_dy[i].end(),d_apolar_dy[i].begin(),d_polar_dy[i].begin(),plus<double>());
   transform(d_polar_dz[i].begin(),d_polar_dz[i].end(),d_apolar_dz[i].begin(),d_polar_dz[i].begin(),plus<double>());
   double hydrophobicity_grid_i=contacts_apolar[i]/total_contacts;
   
   for(unsigned j=0;j<size_protein;j++)
   {
    double d_hydrogrid_dx_ij=(contacts_polar[i]*d_apolar_dx[i][j]-contacts_apolar[i]*d_polar_dx[i][j])/pow(total_contacts,2);
    double d_hydrogrid_dy_ij=(contacts_polar[i]*d_apolar_dy[i][j]-contacts_apolar[i]*d_polar_dy[i][j])/pow(total_contacts,2);
    double d_hydrogrid_dz_ij=(contacts_polar[i]*d_apolar_dz[i][j]-contacts_apolar[i]*d_polar_dz[i][j])/pow(total_contacts,2);
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
