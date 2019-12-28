#include <vector>
#include "tools/PDB.h"

using namespace std;

#ifndef jedi_kernel_h_
#define jedi_kernel_h_

class S_on
{
 private:
    double m_value;
    vector<double> dm_dx;
    vector<double> dm_dy;
    vector<double> dm_dz;
    void compute_m(double &v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz);
 public:
    S_on(double v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz);
    double k;
    double S_on_value;
    vector<double> d_Son_dx;
    vector<double> d_Son_dy;
    vector<double> d_Son_dz;
    void compute_S_on();
};

class S_on_grid
{
 public:
  S_on_grid();
  vector<double> S_on_grid_vector;
  vector<vector<double>> d_Son_grid_dx;
  vector<vector<double>> d_Son_grid_dy;
  vector<vector<double>> d_Son_grid_dz;
  void compute_S_on_grid(vector<double> &mindist, 
                          vector<vector<double>> &d_mindist_dx,
                          vector<vector<double>> &d_mindist_dy,
                          vector<vector<double>> &d_mindist_dz,
                          double &v0, double &delta_v);

};

class S_off
{
 private:
    double m_value;
    vector<double> dm_dx;
    vector<double> dm_dy;
    vector<double> dm_dz;
    void compute_m(double &v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz);
 public:
    S_off(double v, double &v0, double &delta_v, vector<double> &dv_dx, vector<double> &dv_dy, vector<double> &dv_dz);
    double k;
    double S_off_value;
    vector<double> d_Soff_dx;
    vector<double> d_Soff_dy;
    vector<double> d_Soff_dz;
    void compute_S_off();
};

class S_off_grid
{
 public:
  S_off_grid();
  vector<double> S_off_grid_vector;
  vector<vector<double>> d_Soff_grid_dx;
  vector<vector<double>> d_Soff_grid_dy;
  vector<vector<double>> d_Soff_grid_dz;
  void compute_S_off_grid(vector<double> &mindist, 
                         vector<vector<double>> &d_mindist_dx,
                         vector<vector<double>> &d_mindist_dy,
                         vector<vector<double>> &d_mindist_dz,
                         double &v0, double &delta_v);

};

#endif