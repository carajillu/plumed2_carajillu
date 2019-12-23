#include <iostream>
#include "activity.h"
#include "kernel.h"

void Activity::compute_close_contact(vector<double> &mindist, 
                           vector<vector<double>> &d_mindist_dx,
                           vector<vector<double>> &d_mindist_dy,
                           vector<vector<double>> &d_mindist_dz,
                           double &CC_min, double &deltaCC)
{
 //Ugly way to fill vectors with zeros
 for (unsigned i=0; i< mindist.size(); i++)
 {
     close_contact.push_back(0);
 }

 vector<double> d_closecontact_i(0,d_mindist_dx.size());
 for (unsigned i=0; i<mindist.size();i++)
 {
  d_closecontact_dx.push_back(d_closecontact_i);
  d_closecontact_dy.push_back(d_closecontact_i);
  d_closecontact_dz.push_back(d_closecontact_i);
 }

 #pragma omp parallel for
 for (unsigned i=0; i<mindist.size();i++)
 {
   S_on close_contact_i(mindist[i],CC_min,deltaCC,d_mindist_dx[i],d_mindist_dy[i],d_mindist_dz[i]);
   close_contact_i.compute_S_on();
   close_contact[i]=close_contact_i.S_on_value;
   d_closecontact_dx[i]=close_contact_i.d_Son_dx;
   d_closecontact_dy[i]=close_contact_i.d_Son_dy;
   d_closecontact_dz[i]=close_contact_i.d_Son_dz;
 }

 // Uncomment the following lines for testing
 /*
 for (unsigned i=0; i<close_contact.size();i++)
 {
     cout << "Mindist point " << i << ": " << mindist[i] << endl;
     cout << "close_contact point "<< i << ": " << close_contact[i] << endl;
     for (unsigned j=0; j<d_closecontact_dx[i].size();j++)
     {
        cout << "Iteration " << j << ": "<< "derivatives with respect to atom " << j << ": " << d_closecontact_dx[i][j] << " " <<
                                                                    d_closecontact_dy[i][j] << " " <<
                                                                    d_closecontact_dz[i][j] << endl;
     }
 }
 */
}