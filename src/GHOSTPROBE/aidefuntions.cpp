#include "aidefunctions.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include<iterator>
#include <fstream>


using namespace std;

int aidefunctions::findIndex(const vector<PLMD::AtomNumber> arr, PLMD::AtomNumber item) {
    auto ret = find(arr.begin(),arr.end(),item);
    if (ret != arr.end())
        return ret - arr.begin();
    return -1;
}

unsigned aidefunctions::get_random_integer(unsigned start, unsigned end)
{
    random_device rd;                                   // only used once to initialise (seed) engine
    mt19937 rng(rd());                                  // random-number engine used (Mersenne-Twister in this case)
    uniform_int_distribution<unsigned> uni(start,end); // guaranteed unbiased
    auto random_integer = uni(rng);
    return random_integer;
}

vector<vector<double>> aidefunctions::read_xyz(string filename, int frame_id)
{
    vector<vector<double>> xyz;
    vector<double> xyz_j(3);
    int current_frame=-1;
    int linum=-2; // so that 0 is the first line with coordinates
    ifstream fp(filename);
    string line;
    while (getline(fp, line))
    {
     istringstream iss(line);
     istream_iterator<string> beg(iss), end;
     vector<string> tokens(beg, end); // done!
     if (tokens.size() == 1)
     {
      current_frame++;
      linum=-2;
      continue;
     }
     else if (current_frame!=frame_id)
     {
      continue;
     }
     else
     {
      linum++;
      if (linum<0)
          continue;
      xyz_j[0]=atof(tokens[1].c_str())/10;
      xyz_j[1]=atof(tokens[2].c_str())/10;
      xyz_j[2]=atof(tokens[3].c_str())/10;
      xyz.push_back(xyz_j);
      //cout << xyz_j[0] << " " << xyz_j[1] << " " << xyz_j[2] << endl;
     }
    }
    fp.close();
    return xyz;
}

size_t aidefunctions::sample_inverse_weighted_index(const vector<double>& r)
{
    const size_t n = r.size();
    if (n == 0) {
        throw invalid_argument("Input vector r must not be empty.");
    }

    vector<double> weights(n, 0.0);
    size_t positive_count = 0;

    // Parallel: compute weights, discard r[j] == 0
    #pragma omp parallel for reduction(+:positive_count) default(none) shared(r, weights, n)
    for (size_t j = 0; j < n; j++) {
        const double x = r[j];
        if (x > 0.0) {
            weights[j] = 1.0 / x;
            positive_count += 1;
        } else {
            weights[j] = 0.0; // explicitly excluded
        }
    }

    if (positive_count == 0) {
        throw runtime_error("No positive elements in r; cannot select an index.");
    }

    // Serial: random sampling
    static thread_local mt19937 gen{random_device{}()};
    discrete_distribution<size_t> dist(weights.begin(), weights.end());

    return dist(gen);
}

size_t aidefunctions::sample_random_index(const size_t n_atoms)
{
    vector<double> weights(n_atoms,1.0);
    // Serial: random sampling
    static thread_local mt19937 gen{random_device{}()};
    discrete_distribution<size_t> dist(weights.begin(), weights.end());

    return dist(gen);
}