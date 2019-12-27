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
 public:
   Hydrophobicity();
   vector<S_off> hydrophobicity_score;
   void compute_hydrophobicity_score(vector<string> &atomnames, distances &r_matrix, Activity &activity, double &r_hydro, double &deltar_hydro);

};