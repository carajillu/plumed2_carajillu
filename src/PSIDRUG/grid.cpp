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

#include "core/ActionAtomistic.h"
#include "grid.h"

using namespace std;

grid::grid(double &radius, double &spacing)
{
  vector<double> point_i(3,-radius);
  double radius2=pow(radius,2);
  double norm2;
  double centre_x;
  double centre_y;
  double centre_z;
  double r=-radius;
  vector <double> coords;
  while (r<=radius)
  {
    coords.push_back(r);
    r+=spacing;
  }

  double x;
  double y;
  double z;
  for (unsigned i=0; i<coords.size();i++)
  {
    x=coords[i];
    for (unsigned j=0; j<coords.size();j++)
    {
      y=coords[j];
      for (unsigned k=0; k<coords.size();k++)
      {
        z=coords[k];
        cout << x << " " << y << " " << z <<endl;
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