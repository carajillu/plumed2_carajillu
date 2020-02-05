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
                                            vector<vector<double>> d_contacts_total_dx,vector<vector<double>> d_contacts_total_dy,vector<vector<double>> d_contacts_total_dz)
{
  //cout << "Entering hydrophobicity"<< endl;
  unsigned size_grid=contacts_apolar.size();
  unsigned size_protein=d_apolar_dx[0].size();

  Hydrophobicity_grid=vector<double>(size_grid,0);
  vector<double> d_hydrophobicity(size_protein,0);
  d_Hydrophobicity_dx=vector<vector<double>>(size_grid,d_hydrophobicity);
  d_Hydrophobicity_dy=vector<vector<double>>(size_grid,d_hydrophobicity);
  d_Hydrophobicity_dz=vector<vector<double>>(size_grid,d_hydrophobicity);

  for (unsigned i=0; i<size_grid;i++)
  {
   if (total_contacts[i]==0)
   {
     continue; //avoids NaN
   }

   Hydrophobicity_grid[i]=contacts_apolar[i]/total_contacts[i];

   for(unsigned j=0;j<size_protein;j++)
   {
    d_Hydrophobicity_dx[i][j]=(total_contacts[i]*d_apolar_dx[i][j]-d_contacts_total_dx[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
    d_Hydrophobicity_dy[i][j]=(total_contacts[i]*d_apolar_dy[i][j]-d_contacts_total_dy[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
    d_Hydrophobicity_dz[i][j]=(total_contacts[i]*d_apolar_dz[i][j]-d_contacts_total_dz[i][j]*contacts_apolar[i])/pow(total_contacts[i],2);
   }
  }
 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<activity.activity.size();i++)
 {
     cout << "Point " << i << "Activity = " << activity.activity[i] << " h_i = " << Hydrophobicity_grid[i] << endl;
 }
 */
}
