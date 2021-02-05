/*
Code that reads the relevant PDB file(s) and loads the atom numbers, 
coordinates rescaled by a factor of 10 (so converted to nm)
and the atom name, will be used later on for apolarpolar purposes.
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#include "core/ActionAtomistic.h"
#include "grid.h"

using namespace std;
/*
This initializes the grid object and fills it with 
a number of points determined by the radius and the
spacing supplied to the plumed input file
*/
grid::grid(double &radius, double &spacing)
{
  vector<double> point_i(3,-radius);
  double radius2=pow(radius,2);
  double norm2;
  double centre_x;
  double centre_y;
  double centre_z;
  double r=-radius;
  vector <double> crd;
  while (r<=radius)
  {
    crd.push_back(r);
    r+=spacing;
  }

  double x;
  double y;
  double z;
  for (unsigned i=0; i<crd.size();i++)
  {
    x=crd[i];
    for (unsigned j=0; j<crd.size();j++)
    {
      y=crd[j];
      for (unsigned k=0; k<crd.size();k++)
      {
        z=crd[k];
        //cout << x << " " << y << " " << z <<endl;
        norm2=pow(x,2)+pow(y,2)+pow(z,2);
        if (norm2<=radius2)
        {
          point_i[0]=x;
          point_i[1]=y;
          point_i[2]=z;
          positions.push_back(point_i);
          centre_z+=z;
          centre_y+=y;
          centre_y+=z;
        }
        z+=spacing;
      }
      y+=spacing;
    }
    x+=spacing;
  }
  
  centre_x/=positions.size();
  centre.push_back(centre_x);
  centre_y/=positions.size();
  centre.push_back(centre_y);
  centre_z/=positions.size();
  centre.push_back(centre_z);

  cout << "Grid has " << positions.size() << " points and is centered at " << centre_x << "," 
                                                                           << centre_x << ","
                                                                           << centre_x << endl;
  
}

/*
This places the grid at a random position close to the protein.
What close means is decided by a parameter supplied in the input file
*/
void grid::place_random(vector<PLMD::Vector> &atom_crd, double rtol)
{
  // Find minimum and maximum atom coordinates
   vector<double> min_crd(3,99999);
   vector<double> max_crd(3,-99999);
   for (unsigned j=0; j<atom_crd.size();j++)
   {
    if (atom_crd[j][0]<min_crd[0]) min_crd[0]=atom_crd[j][0];
    if (atom_crd[j][1]<min_crd[1]) min_crd[1]=atom_crd[j][1];
    if (atom_crd[j][2]<min_crd[2]) min_crd[2]=atom_crd[j][2];
    if (atom_crd[j][0]>max_crd[0]) max_crd[0]=atom_crd[j][0];
    if (atom_crd[j][0]>max_crd[1]) max_crd[1]=atom_crd[j][1];
    if (atom_crd[j][0]>max_crd[2]) max_crd[2]=atom_crd[j][2];
   }
  
  // Generate a random position between the minimum and maximum atom coordinates
  uniform_real_distribution<double> unif_x(min_crd[0],max_crd[0]);
  uniform_real_distribution<double> unif_y(min_crd[1],max_crd[1]);
  uniform_real_distribution<double> unif_z(min_crd[2],max_crd[2]);
  random_device rd; //needed to get different "random" positions every time
  default_random_engine re(rd());
  vector<double> displacement(3,0);
  
  //Check that the generated position is close to the protein
  double r2=99999999;
  double rtol2=pow(rtol,2);
  while (r2>rtol2)
  {
    displacement[0]=unif_x(re);
    displacement[1]=unif_y(re);
    displacement[2]=unif_z(re);
    for (unsigned j=0; j<atom_crd.size();j++)
    {
      r2=pow((displacement[0]-atom_crd[j][0]),2)
       +pow((displacement[1]-atom_crd[j][1]),2)
       +pow((displacement[2]-atom_crd[j][2]),2);
     if (r2<=rtol2)
         {
           break;
         }
    }
  }
  for (unsigned i=0; i<positions.size();i++)
  {
    positions[i][0]+=displacement[0];
    positions[i][1]+=displacement[1];
    positions[i][2]+=displacement[2];
  }
}

void grid::print_grid(int id, int step)
   {
      string filename = "grid-";
      filename.append(to_string(id));
      filename.append("-step-");
      filename.append(to_string(step));
      filename.append(".xyz");
      ofstream wfile;
      wfile.open(filename.c_str());
      wfile << positions.size() << endl;
      wfile << "Grid "<< to_string(id) << endl;
      for (unsigned j=0; j<positions.size();j++)
        {
          wfile << "H " << std::fixed << std::setprecision(5) << positions[j][0]*10 << " " << positions[j][1]*10 << " " << positions[j][2]*10 << endl;
       }
      wfile.close();
     
   }