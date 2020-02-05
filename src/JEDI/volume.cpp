#include <iostream>
#include "volume.h"
#include "activity.h"
using namespace std;

Volume::Volume()
{
}

void Volume::compute_hydroactivity(vector<double> activity_grid,
                             vector<vector<double>> d_activity_dx,vector<vector<double>> d_activity_dy,vector<vector<double>> d_activity_dz,
                             vector<double> hydrophobicity_grid,
                             vector<vector<double>> d_hydrophobicity_dx,vector<vector<double>> d_hydrophobicity_dy,vector<vector<double>> d_hydrophobicity_dz)
{
 unsigned size_grid=activity_grid.size();
 unsigned size_protein=d_activity_dx[0].size();
 hydroactivity_grid=vector<double>(size_grid,0);
 vector<double> d_hydroactivity(size_protein,0);
 d_hydroactivity_dx=vector<vector<double>>(size_grid,d_hydroactivity);
 d_hydroactivity_dy=vector<vector<double>>(size_grid,d_hydroactivity);
 d_hydroactivity_dz=vector<vector<double>>(size_grid,d_hydroactivity);

 for (unsigned i=0; i<size_grid;i++)
 {
   hydroactivity_grid[i]=activity_grid[i]*hydrophobicity_grid[i];
   for (unsigned j=0; j<size_protein;j++)
   {
     d_hydroactivity_dx[i][j]= hydrophobicity_grid[i]*d_activity_dx[i][j]+activity_grid[i]*d_hydrophobicity_dx[i][j];
     d_hydroactivity_dy[i][j]= hydrophobicity_grid[i]*d_activity_dy[i][j]+activity_grid[i]*d_hydrophobicity_dy[i][j];
     d_hydroactivity_dz[i][j]= hydrophobicity_grid[i]*d_activity_dz[i][j]+activity_grid[i]*d_hydrophobicity_dz[i][j];
   }
 }
}

void Volume::compute_volume(vector<double> hydroactivity_grid,
                            vector<vector<double>> d_hydroactivity_dx,vector<vector<double>> d_hydroactivity_dy,vector<vector<double>> d_hydroactivity_dz,
                            double resolution)
{
unsigned size_grid=hydroactivity_grid.size();
unsigned size_protein=d_hydroactivity_dx[0].size();
volume=0;
d_volume_dx=vector<double>(size_protein,0);
d_volume_dz=vector<double>(size_protein,0);
d_volume_dy=vector<double>(size_protein,0);
for (unsigned i=0; i<size_grid; i++)
 {
   volume+=hydroactivity_grid[i]*pow(resolution,3);
  for (unsigned j=0; j<size_protein;j++)
  {
   d_volume_dx[j]+=d_hydroactivity_dx[i][j]*pow(resolution,3);
   d_volume_dz[j]+=d_hydroactivity_dy[i][j]*pow(resolution,3);
   d_volume_dy[j]+=d_hydroactivity_dz[i][j]*pow(resolution,3);
  }
 }
}