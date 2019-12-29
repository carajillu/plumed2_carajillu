#include <iostream>
#include "hydrophobicity.h"

using namespace std;

Hydrophobicity::Hydrophobicity()
{

};

void Hydrophobicity::compute_contacts(distances &r_matrix, double &r_hydro, double &deltar_hydro)
{
 int size_grid  = r_matrix.r_matrix.size();
 int atomnumber = r_matrix.r_matrix[0].size();
 

 vector<double> contacts_i(atomnumber,0);
 for (unsigned i=0; i<size_grid;i++)
  {
   contacts.push_back(contacts_i);
   d_contacts_dx.push_back(contacts_i);
   d_contacts_dy.push_back(contacts_i);
   d_contacts_dz.push_back(contacts_i);
  }

 #pragma omp parallel for
 for (unsigned i=0;i<size_grid;i++)
 {
    for (unsigned j=0; j<atomnumber; j++)
    {
      // The S_off class has a vector<double> with length = atomnumber for the derivatives in each direction,
      vector<double> drij_dx(atomnumber,0);
      drij_dx[j]=r_matrix.dr_matrix_dx[i][j];
      vector<double> drij_dy(atomnumber,0);
      drij_dy[j]=r_matrix.dr_matrix_dy[i][j];
      vector<double> drij_dz(atomnumber,0);
      drij_dz[j]=r_matrix.dr_matrix_dz[i][j];
      
      S_off contact_ij(r_matrix.r_matrix[i][j],r_hydro,deltar_hydro,drij_dx,drij_dy,drij_dz);
      contact_ij.compute_S_off();
      contacts[i][j]=contact_ij.S_off_value;
      
      d_contacts_dx[i][j]=contact_ij.d_Soff_dx[j];
      d_contacts_dy[i][j]=contact_ij.d_Soff_dy[j];
      d_contacts_dz[i][j]=contact_ij.d_Soff_dz[j];
      
    }
 }

 //Uncomment the following lines for testing
 /*
 for (unsigned i=0;i<size_grid;i++)
 {
  cout << "grid point " << i << endl;
  for (unsigned j=0; j<atomnumber; j++)
  {
    cout << "rij = " << r_matrix.r_matrix[i][j] << " S_off rij with atom " << j << ": " << contacts[i][j] << endl;
    cout << "Derivatives with respect to atom " << j << " " << d_contacts_dx[i][j] << " " <<
                                                               d_contacts_dy[i][j] << " " <<
                                                               d_contacts_dz[i][j] << " " << endl;
  }
 }
 */
};

void Hydrophobicity::compute_hydrophobicity_grid(vector<string> &atomnames, distances &r_matrix, double &r_hydro, double &deltar_hydro)
{
 int size_grid=r_matrix.r_matrix.size();
 int atomnumber = atomnames.size();
 
 //Ugly way to fill vectors
 vector<double> d_hydro_i(atomnames.size(),0);
 for (unsigned i=0; i<size_grid;i++)
 {
     hydrophobicity_grid.push_back(0);
     d_hydrogrid_dx.push_back(d_hydro_i);
     d_hydrogrid_dy.push_back(d_hydro_i);
     d_hydrogrid_dz.push_back(d_hydro_i);
 }

 //#pragma omp parallel for
 for (unsigned i=0; i<size_grid;i++)
 {
     double apolar_contacts=0;
     vector<double> d_apolarcontacts_dx=d_hydro_i;
     vector<double> d_apolarcontacts_dy=d_hydro_i;
     vector<double> d_apolarcontacts_dz=d_hydro_i;
     double polar_contacts=0;
     vector<double> d_polarcontacts_dx=d_hydro_i;
     vector<double> d_polarcontacts_dy=d_hydro_i;
     vector<double> d_polarcontacts_dz=d_hydro_i;

     for (unsigned j=0; j<atomnames.size();j++)
     {
        if (atomnames[j][0]=='C' or atomnames[j][0]=='S')
         {
          apolar_contacts+=contacts[i][j];
          d_apolarcontacts_dx[j]+=d_contacts_dx[i][j];
          d_apolarcontacts_dy[j]+=d_contacts_dy[i][j];
          d_apolarcontacts_dz[j]+=d_contacts_dz[i][j];
         }
         else if (atomnames[j][0]=='N' or atomnames[j][0]=='O')
         {
          polar_contacts+=contacts[i][j];
          d_polarcontacts_dx[j]+=d_contacts_dx[i][j];
          d_polarcontacts_dy[j]+=d_contacts_dy[i][j];
          d_polarcontacts_dz[j]+=d_contacts_dz[i][j];
          
         }
         
         else
         {
            cout << "Unknown atom detected at position: " << j << ": " << atomnames[j] << ". Exiting." << endl;
            exit(0);
         }
         
     }

     double total_contacts=apolar_contacts+polar_contacts;
     if (total_contacts > 0) 
     {
      hydrophobicity_grid[i]=apolar_contacts/(total_contacts); //otherwise may give NaN
      for (unsigned j=0; j<atomnames.size();j++)
      {
       d_hydrogrid_dx[i][j]=(polar_contacts*d_apolarcontacts_dx[j]-apolar_contacts*d_polarcontacts_dx[j])/pow(total_contacts,2);
       d_hydrogrid_dy[i][j]=(polar_contacts*d_apolarcontacts_dy[j]-apolar_contacts*d_polarcontacts_dy[j])/pow(total_contacts,2);
       d_hydrogrid_dz[i][j]=(polar_contacts*d_apolarcontacts_dz[j]-apolar_contacts*d_polarcontacts_dz[j])/pow(total_contacts,2);
      }
     }
 }

 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<size_grid; i++)
 {
   cout << "Hydrophobicity of grid point " << i << "= " << hydrophobicity_grid[i] << endl;
   for (unsigned j=0; j<atomnames.size(); j++)
   {
        cout << "Derivatives with respect to atom " << j << " with name " << atomnames[j] << ": " << d_hydrogrid_dx[i][j] << " " <<
                                                                   d_hydrogrid_dy[i][j] << " " <<
                                                                   d_hydrogrid_dz[i][j] << " " << endl;
   }
 }
 */

};

void Hydrophobicity::compute_hydrophobicity(vector<string> &atomnames, distances &r_matrix, Activity &activity, double &r_hydro, double &deltar_hydro)
{
 compute_contacts(r_matrix,r_hydro,deltar_hydro);
 compute_hydrophobicity_grid(atomnames, r_matrix, r_hydro, deltar_hydro);
 
 //Ugly way to fill vectors
 vector<double> d_hydro(atomnames.size(),0);
 d_Ha_dx=d_hydro;
 d_Ha_dy=d_hydro;
 d_Ha_dz=d_hydro;

 double activity_sum=0;
 vector<double> d_activity_sum_dx(atomnames.size());
 vector<double> d_activity_sum_dy(atomnames.size());
 vector<double> d_activity_sum_dz(atomnames.size());

 double hydro_activity_sum=0;
 vector<double> d_hydroactivity_sum_dx(atomnames.size());
 vector<double> d_hydroactivity_sum_dy(atomnames.size());
 vector<double> d_hydroactivity_sum_dz(atomnames.size());

 #pragma omp parallel for reduction(+:activity_sum,hydro_activity_sum)
 for (unsigned i=0; i<activity.activity.size();i++)
 {
   activity_sum+=activity.activity[i];
   hydro_activity_sum+=hydrophobicity_grid[i]*activity.activity[i];
   for (unsigned j=0; j<atomnames.size();j++)
   {
    d_activity_sum_dx[j]+=activity.d_activity_dx[i][j];
    d_activity_sum_dy[j]+=activity.d_activity_dy[i][j];
    d_activity_sum_dz[j]+=activity.d_activity_dz[i][j];
    
    d_hydroactivity_sum_dx[j]+=activity.activity[i]*d_hydrogrid_dx[i][j]+hydrophobicity_grid[i]*activity.d_activity_dx[i][j];
    d_hydroactivity_sum_dy[j]+=activity.activity[i]*d_hydrogrid_dy[i][j]+hydrophobicity_grid[i]*activity.d_activity_dy[i][j];
    d_hydroactivity_sum_dz[j]+=activity.activity[i]*d_hydrogrid_dz[i][j]+hydrophobicity_grid[i]*activity.d_activity_dz[i][j];
   }
 }

 if (activity_sum>0) // Otherwise may give NaN
 {
  Ha=hydro_activity_sum/activity_sum; 

  #pragma omp parallel for
  for (unsigned j=0; j<atomnames.size();j++)
  {
   d_Ha_dx[j]=(activity_sum*d_hydroactivity_sum_dx[j]-hydro_activity_sum*d_activity_sum_dx[j])/pow(activity_sum,2);
   d_Ha_dy[j]=(activity_sum*d_hydroactivity_sum_dy[j]-hydro_activity_sum*d_activity_sum_dy[j])/pow(activity_sum,2);
   d_Ha_dz[j]=(activity_sum*d_hydroactivity_sum_dz[j]-hydro_activity_sum*d_activity_sum_dz[j])/pow(activity_sum,2);
  }
 }
 else  Ha=0.;
  
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
