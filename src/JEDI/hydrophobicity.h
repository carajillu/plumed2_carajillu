#include "tools/PDB.h"
#include "distances.h"
#include "kernel.h"
#include "activity.h"

using namespace std;

class Hydrophobicity
{
 private:
   vector<vector<S_off>> contacts;
   void compute_contacts(distances &r_matrix, double &r_hydro, double &deltar_hydro);
   vector<double> hydrophobicity_grid;
   vector<vector<double>> d_hydrogrid_dx;
   vector<vector<double>> d_hydrogrid_dy;
   vector<vector<double>> d_hydrogrid_dz;
   void compute_hydrophobicity_grid(vector<string> &atomnames, distances &r_matrix, double &r_hydro, double &deltar_hydro);
 public:
   Hydrophobicity();
   double hydrophobicity;
   vector<double> d_hydrophobicity_dx;
   vector<double> d_hydrophobicity_dy;
   vector<double> d_hydrophobicity_dz;
   void compute_hydrophobicity(vector<string> &atomnames, distances &r_matrix, Activity &activity, double &r_hydro, double &deltar_hydro);
   
   

};