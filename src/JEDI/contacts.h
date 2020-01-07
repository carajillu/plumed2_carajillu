#include <vector>
#include "tools/PDB.h"

using namespace std;

#ifndef jedi_contacts_h_
#define jedi_contacts_h_

class contacts_Soff
{
public:
 contacts_Soff();
 vector<vector<double>> contacts_matrix;
 vector<vector<double>> d_contacts_dx;
 vector<vector<double>> d_contacts_dy;
 vector<vector<double>> d_contacts_dz;
 void compute_contacts_S_off(vector<vector<double>> r_matrix, 
                             vector<vector<double>> dr_matrix_dx, 
                             vector<vector<double>> dr_matrix_dy, 
                             vector<vector<double>> dr_matrix_dz,
                             double CC2min, double deltaCC2);

};

class contacts_sum
{
 public:
 contacts_sum();
 vector<double> contacts_apolar;
 vector<vector<double>> d_contacts_apolar_dx;
 vector<vector<double>> d_contacts_apolar_dy;
 vector<vector<double>> d_contacts_apolar_dz;

 vector<double> contacts_polar;
 vector<vector<double>> d_contacts_polar_dx;
 vector<vector<double>> d_contacts_polar_dy;
 vector<vector<double>> d_contacts_polar_dz;

 void compute_contacts_sum(vector<vector<double>> contacts_matrix,
                           vector<vector<double>> d_contacts_dx,
                           vector<vector<double>> d_contacts_dy,
                           vector<vector<double>> d_contacts_dz,
                           vector<string> atomnames);
};
#endif