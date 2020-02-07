/*
Class that reads the JEDI parameters file
and transfers the parameters to the main
code in the form of a jediparameters object
*/

using namespace std;


#ifndef _jediparm_h_
#define _jediparm_h_
class jediparameters
{
public:
  jediparameters();
  void readParams(string parameters_file);
  double alpha;
  double beta;
  double gamma;
  double theta;
  double CC_mind;
  double deltaCC;
  double Emin;
  double deltaE;
  double BSmin;
  double deltaBS;
  double CC2_min;
  double deltaCC2;
  double GP_min;
  double GP_max;
  double r_hydro;
  double deltar_hydro;
  double V_max;
  double deltaV_max;
  double V_min;
  double deltaV_min;
  double resolution;
  double r_max;
  double r_max_clust;
};
#endif