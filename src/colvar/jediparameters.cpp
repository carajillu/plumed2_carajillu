#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <cmath>
#include <iterator>

#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"

#include "jediparameters.h"

using namespace std;

jediparameters::jediparameters()
{
  alpha = 0.0;
  beta = 0.0;
  gamma = 0.0;
  theta = 0.0;
  CC_mind = 0.0;
  deltaCC = 0.0;
  Emin = 0.0;
  deltaE = 0.0;
  BSmin = 0.0;
  deltaBS = 0.0;
  CC2_min = 0.0;
  deltaCC2 = 0.0;
  GP_min = 0.0;
  GP_max = 0.0;
  r_hydro = 0.0;
  deltar_hydro = 0.0;
  V_max = 0.0;
  deltaV_max = 0.0;
  V_min = 0.0;
  deltaV_min = 0.0;
  resolution = 0.0;
}

bool jediparameters::readParams(string &parameters_file)
{
  FILE* fp=fopen(parameters_file.c_str(),"r");
  if (!fp) return false;
  string line;
  while (PLMD::Tools::getline(fp, line))
    {
      if (line[0] == '#')
	continue;
      //cout << line << endl;
      istringstream iss(line);
      istream_iterator<string> beg(iss), end;
      vector<string> tokens(beg, end); // done!
      if (tokens.size() < 2)
	{
	  cout << "ABORT ! JEDI parameters file appear malformed check the line below !" << endl; 
	  cout << line << endl;
	  exit(-1);
	}
      string key = tokens[0];
      double item = atof(tokens[2].c_str());
      // set jedi parameters datastruct
      // FIXME: Check that all parameters have been set
      if ( key == string("alpha") )
	alpha = item;
      else if ( key == string("beta") )
	beta = item;
      else if ( key == string("gamma") )
	gamma = item;
      else if ( key == string("theta") )
	theta = item;
      else if ( key == string("CC_mind") )
	CC_mind = item;
      else if ( key == string("deltaCC") )
	deltaCC = item;
      else if ( key == string("Emin") )
	Emin = item;
      else if ( key == string("deltaE") )
	deltaE = item;
      else if ( key == string("BSmin") )
	BSmin = item;
      else if ( key == string("deltaBS") )
	deltaBS = item;
      else if ( key == string("CC2_min") )
	CC2_min = item;
      else if ( key == string("deltaCC2") )
	deltaCC2 = item;
      else if ( key == string("GP_min") )
	GP_min = item;
      else if ( key == string("GP_max") )
	GP_max = item;
      else if ( key == string("r_hydro") )
	r_hydro = item;
      else if ( key == string("deltar_hydro") )
	deltar_hydro = item;
      else if ( key == string("V_max") )
	V_max = item;
      else if ( key == string("deltaV_max") )
	deltaV_max = item;
      else if ( key == string("V_min") )
	V_min = item;
      else if ( key == string("deltaV_min") )
	deltaV_min = item;
      else if ( key == string("grid_resolution") )
	resolution = item;
    }
  fclose(fp);

  cout << "*** Values of the JEDI parameter set loaded in memory: ***" << endl;
  cout << "alpha = " << alpha << endl;
  cout << "beta  = " << beta << endl;
  cout << "gamma = " << gamma << endl;
  cout << "theta = " << theta << endl;
  cout << "CC_mind  = " << CC_mind << endl;
  cout << "deltaCC  = " << deltaCC << endl;
  cout << "Emin  = " << Emin << endl;
  cout << "deltaE  = " << deltaE << endl;
  cout << "BSmin  = " << BSmin << endl;
  cout << "deltaBS  = " << deltaBS << endl;
  cout << "CC2_min  = " << CC2_min << endl;
  cout << "deltaCC2  = " << deltaCC2 << endl;
  cout << "GP_min  = " << GP_min << endl;
  cout << "GP_max  = " << GP_max << endl;
  cout << "r_hydro  = " << r_hydro << endl;
  cout << "deltar_hydro  = " << deltar_hydro << endl;
  cout << "V_max  = " << V_max << endl;
  cout << "deltaV_max  = " << deltaV_max << endl;
  cout << "V_min  = " << V_min << endl;
  cout << "deltaV_min  = " << deltaV_min << endl;
  cout << "grid_resolution = " << resolution << endl;
  
  return true;
}