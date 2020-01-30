#include "clustering_laio.h"
#include <algorithm>
#include <iostream> 
#include <iomanip> // std::setprecision
#include <fstream> // std::ofstream     
using namespace std;

clustering::clustering()
{

}

bool clustering::sort_grid(const laio &a, const laio &b)
{
  if (a.density!=b.density)
    return a.density>b.density;
  else if (a.activity!=b.activity)
    return a.activity>b.activity;
  else 
    return a.idx > b.idx;
} 

void clustering::cluster_grid(vector<double> activity, vector<vector<double>> r_matrix,vector<vector<unsigned>> neighbours, double GP_max, double grid_resolution)
{
  int clust_id=0;
  #pragma omp parallel shared(grid_stats, clust_id)
  {
    vector<laio> grid_stats_private;
    // Build the vector with all the data dor each grid point and calculate all densities
    #pragma omp for
    for (unsigned i=0; i<activity.size();i++)
    {
     if (activity[i]==0) continue;
     laio point_stats;
     point_stats.idx=i;
     point_stats.activity=activity[i];
     point_stats.density=0;
     point_stats.nnhd=-1;
     point_stats.delta=9999999;
     point_stats.cluster=-1;
     for(unsigned k=0; k<neighbours[i].size();k++)
     {
       if (activity[k]>0) point_stats.density++;
     }
     grid_stats_private.push_back(point_stats);
    }
    #pragma omp critical
    {
     grid_stats.insert(grid_stats.end(),grid_stats_private.begin(),grid_stats_private.end()); 
    }
    #pragma omp barrier

    //Sort the vector<laio> using the criteria in clustering::sort_grid
    #pragma omp single
    {
      sort(grid_stats.begin(),grid_stats.end(),sort_grid);
    }
    #pragma omp barrier

    //Calculate delta for each grid point
    #pragma omp for
    for (unsigned i=1; i<grid_stats.size();i++)
    {
      unsigned idx_i=grid_stats[i].idx;
      
      //Calculate delta
      for(unsigned k=0; k<i; k++)
      {
       unsigned idx_k=grid_stats[k].idx;
       if (r_matrix[idx_i][idx_k]==grid_resolution) // that speeds up the search
       {
        grid_stats[i].delta=grid_resolution;
        grid_stats[i].nnhd=k; //Note that here we get the position in the vector<laio>, not the actual index of the grid point
        break;
       }
       else if (r_matrix[idx_i][idx_k]<grid_stats[i].delta)
       {
        grid_stats[i].delta=r_matrix[idx_i][idx_k];
        grid_stats[i].nnhd=k;
       }
      }
    }
  }

  //Assign points to clusters (note that's out of the parallel region!)
  
  //The point with the highest density always is a cluster center
    
  grid_stats[0].cluster=0;
  vector<unsigned> cluster(1,grid_stats[0].idx);
  clusters.push_back(cluster);

  for (unsigned i=1; i<grid_stats.size();i++)
  {
    unsigned idx_i=grid_stats[i].idx;
    //Assign cluster centers
    if (grid_stats[i].delta > GP_max)
      {
        clust_id++;  
        grid_stats[i].cluster=clust_id;
        vector<unsigned> cluster(1,idx_i);
        clusters.push_back(cluster);
      }   
    else
      {
        int nnhd=grid_stats[i].nnhd;
        grid_stats[i].cluster=grid_stats[nnhd].cluster;
        clusters[grid_stats[i].cluster].push_back(idx_i);
      }
      
  }
/*
  for (unsigned i=0; i<grid_stats.size();i++)
  {
    cout << grid_stats[i].idx << " " << grid_stats[i].activity << " " << grid_stats[i].density << " " << grid_stats[i].delta << " " << grid_stats[i].nnhd << " " << grid_stats[i].cluster << endl;
  }
  
  int count=0;
  for (unsigned i=0; i<clusters.size();i++)
  {
    cout << "Cluster " << i << " Has " << clusters[i].size() << " elements" << endl;
    count+=clusters[i].size();
  }
  cout << "Total points clustered: " << count << endl;
  cout << "Total number of grid points " << activity.size() << endl;
  exit(0);
  */
}

void clustering::print_clusters(vector<PLMD::Vector> grid)
{
 for (unsigned k=0; k<clusters.size();k++)
     {
        string filename = "cluster-";
        stringstream num;
        num << k;
        string number = num.str();
        //stringstream out;
        //out << step;
        //string outer=out.str();
        filename.append(number);
        filename.append("-step-");
        //filename.append(outer); 
        filename.append(".xyz");
        ofstream wfile;
        wfile.open(filename.c_str());
        wfile << clusters[k].size() << endl;
        //wfile << jedi_clusters[k] << endl;
        wfile << "Cluster " << k << endl;
        for (unsigned i=0; i<clusters[k].size();i++)
          {
            wfile << "H " << std::fixed << std::setprecision(5) << grid[clusters[k][i]][0]*10 << " " << grid[clusters[k][i]][1]*10 << " " << grid[clusters[k][i]][2]*10 << endl;
          }
        wfile.close();
     }
}