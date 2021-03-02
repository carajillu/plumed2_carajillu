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
#include <iterator>

#include "core/ActionAtomistic.h"
#include "grid.h"

using namespace std;
/*
Initialise the grid object
*/
grid::grid(unsigned n_atoms)
{
 // Fill bsite_bin with zeros
 bsite_bin=vector<unsigned>(n_atoms,1);
 PsiGrid=0;
 d_Psigrid_dx=vector<double>(n_atoms,0);
 d_Psigrid_dy=vector<double>(n_atoms,0);
 d_Psigrid_dz=vector<double>(n_atoms,0);
}


/*
Reads a grid from an xyz file and skip the coordinate
generation step. Only used for debugging at the moment
but this could change in the future (maybe post-processing?)
*/

void grid::grid_read(string grid_file)
{
 cout << "Reading Grid" << endl;
 centre=vector<double>(3,0);

 vector<double> point_crd(3,0);
 FILE* fp=fopen(grid_file.c_str(),"r");
 string line;
 int linum=0;
 while (PLMD::Tools::getline(fp, line))
 {
  if (linum<2) //skip first 2 lines of xyz file
  {
  linum++;
  continue;
  }
  istringstream iss(line);
  std::istream_iterator<string> beg(iss), end;
  vector<string> tokens(beg, end);
  if (tokens.size()<4) continue;
  point_crd[0]=atof(tokens[1].c_str())/10;
  point_crd[1]=atof(tokens[2].c_str())/10;
  point_crd[2]=atof(tokens[3].c_str())/10;
  positions.push_back(point_crd);

  centre[0]+=point_crd[0];
  centre[1]+=point_crd[1];
  centre[2]+=point_crd[2];
 }

 centre[0]/=positions.size();
 centre[1]/=positions.size();
 centre[2]/=positions.size();

 size_grid=positions.size();
}

/*
Sets up the grid object with coordinates, center and 
initialises bsite_bin. All is determined by the radius 
and the spacing supplied to the plumed input file
*/
void grid::grid_setup(double &radius, double &spacing, unsigned &n_atoms)
{
  vector<double> point_i(3,-radius);
  double radius2=pow(radius,2);
  double norm2;
  double centre_x;
  double centre_y;
  double centre_z;
  double r=-radius;

  // Generate the vector the grid will be built from (1D)
  vector <double> crd;
  while (r<=radius)
  {
    crd.push_back(r);
    r+=spacing;
  }

  // Build the grid
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
        norm2=pow(x,2)+pow(y,2)+pow(z,2);
        if (norm2<=(radius2+spacing/2))
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
  
  size_grid=positions.size();

  std::cout << "Grid has " << positions.size() << " points and is centered at " << centre_x << "," 
                                                                           << centre_x << ","
                                                                           << centre_x << endl;
}
/*
This places the grid at a random position close to the protein.
What close means is decided by a parameter supplied in the input file
*/
void grid::place_random(vector<PLMD::Vector> &atom_crd, double &rtol)
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
  centre=displacement;
}

// Assign values to vector<unsigned> bsite_bin
//bsite_bin[j]=1 if atom j is closer than or at RSITE nm from the center of the grid
//biste_bool[j]=0 if atom j is further away than RSITE nm from the center of the grid
void grid::assign_bsite_bin(vector<PLMD::Vector> &atom_crd, double &rsite)
{
 double rsite2=pow(rsite,2);
 double r2;
 for (unsigned j=0; j<atom_crd.size();j++)
  {
   r2=pow((centre[0]-atom_crd[j][0]),2)+pow((centre[1]-atom_crd[j][1]),2)+pow((centre[2]-atom_crd[j][2]),2);
   if (r2<=rsite2)
   {
      bsite_bin[j]=1;
   }
   else
      bsite_bin[j]=0;
  }
}

// Center the grid at the center of geometry (not center of mass) of the bining site
void grid::center_grid(vector<PLMD::Vector> &atom_crd)
{
 double x=0;
 double y=0;
 double z=0;
 int total_bsite=0;
 for (unsigned j=0; j<atom_crd.size();j++)
 {
  total_bsite+=bsite_bin[j]; // That gives us the total of elements in bsite_bin that are equal to 1
  x+=bsite_bin[j]*atom_crd[j][0];
  y+=bsite_bin[j]*atom_crd[j][1];
  z+=bsite_bin[j]*atom_crd[j][2];
 }
 x/=total_bsite;
 y/=total_bsite;
 z/=total_bsite;
 
 for (unsigned i=0; i<positions.size();i++)
 {
   positions[i][0]+=(x-centre[0]);
   positions[i][1]+=(y-centre[1]);
   positions[i][2]+=(z-centre[2]);
 }
 centre[0]=x;
 centre[1]=y;
 centre[2]=z;
}
// Output the i-th grid at the current time step in an xyz file
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

void grid::init_psigrid(unsigned &n_atoms)
{
  PsiGrid=0;
  d_Psigrid_dx=vector<double>(n_atoms,0);
  d_Psigrid_dy=vector<double>(n_atoms,0);
  d_Psigrid_dz=vector<double>(n_atoms,0);
}

void grid::add_activity(double &activity,
                        vector<double> &d_activity_dx,
                        vector<double> &d_activity_dy,
                        vector<double> &d_activity_dz)
{
  PsiGrid+=activity;
  transform (d_Psigrid_dx.begin(), d_Psigrid_dx.end(), 
             d_activity_dx.begin(), d_Psigrid_dx.begin(), 
             std::plus<double>());
  transform (d_Psigrid_dy.begin(), d_Psigrid_dy.end(), 
             d_activity_dy.begin(), d_Psigrid_dy.begin(), 
             std::plus<double>());
  transform (d_Psigrid_dz.begin(), d_Psigrid_dz.end(), 
             d_activity_dz.begin(), d_Psigrid_dz.begin(), 
             std::plus<double>());
}