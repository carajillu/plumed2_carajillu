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

 S_off contact_ij(r_matrix.r_matrix[0][0],r_hydro,deltar_hydro,r_matrix.dr_matrix_dx[0],r_matrix.dr_matrix_dy[0],r_matrix.dr_matrix_dz[0]);
 vector<S_off> contacts_i(atomnumber,contact_ij);
 //Ugly way to fill vectors
 for (unsigned i=0;i<size_grid;i++)
 {
     contacts.push_back(contacts_i);
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
      contacts[i][j]=contact_ij;
    }
 }

 //Uncomment the following lines for testing
 for (unsigned i=0;i<size_grid;i++)
 {
  cout << "grid point " << i << endl;
  for (unsigned j=0; j<atomnumber; j++)
  {
    cout << "rij = " << r_matrix.r_matrix[i][j] << " S_off rij with atom " << j << ": " << contacts[i][j].S_off_value << endl;
    for (unsigned l=0; l<atomnumber; l++)
    {
        cout << "Derivatives with respect to atom " << l << " " << contacts[i][j].d_Soff_dx[l] << " " <<
                                                                   contacts[i][j].d_Soff_dy[l] << " " <<
                                                                   contacts[i][j].d_Soff_dz[l] << " " << endl;
    }
  }
  exit(0);
 }
};

void Hydrophobicity::compute_hydrophobicity_score(vector<string> &atomnames, distances &r_matrix, Activity &activity, double &r_hydro, double &deltar_hydro)
{
 compute_contacts(r_matrix,r_hydro,deltar_hydro);
};