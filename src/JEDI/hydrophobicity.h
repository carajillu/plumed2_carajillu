#include "tools/PDB.h"
#include "distances.h"
#include "kernel.h"
#include "activity.h"

using namespace std;

#ifndef hydrophobicity_h_
#define hydrophobicity_h_

class hydrophobicity
{
 public:
   hydrophobicity();
   vector<double>  Hydrophobicity_grid;
   vector<vector<double>> d_Hydrophobicity_dx;
   vector<vector<double>> d_Hydrophobicity_dy;
   vector<vector<double>> d_Hydrophobicity_dz;
   void compute_hydrophobicity(vector<double> contacts_apolar,
                               vector<vector<double>> d_apolar_dx,vector<vector<double>> d_apolar_dy,vector<vector<double>> d_apolar_dz,
                               vector<double> total_contacts,
                               vector<vector<double>> d_contacts_total_dx,vector<vector<double>> d_contacts_total_dy,vector<vector<double>> d_contacts_total_dz);
   
   

};
#endif