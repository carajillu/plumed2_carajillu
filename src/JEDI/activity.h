#include<vector>
#include "tools/PDB.h"

using namespace std;

#ifndef activity_h_
#define activity_h_

class Activity
{
 private:

   vector<double> ligand_correction;
   vector<double> d_ligand_correction;

   vector<double> close_contact;
   vector<double> d_close_contact;

   
   vector<double> depth;
   vector<double> d_depth;


 public:
   Activity();
   vector<double> activity;
   vector<double> d_activity; //d_lig_i*(close_contact*depth)+((d_close_contact*depth)+(close_contact*d_depth))*lig_i
   vector<vector<int>> neighbours; //move to private after debug
   
   void compute_neighbours(vector<PLMD::Vector> &positions,double &GPmin, double &GPmax);
   
};

#endif