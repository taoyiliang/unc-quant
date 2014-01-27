#include <cmath>
#include <vector>
#include <string>
#include <sstream>
////
// get the neiboring cells in terms of global id, also returns i and j index corresponds to global_id.

namespace FVHelper
{
  void get_neighbors(const int global_id, int &i, int &j, int nx, int ny, int &left, int &right, int &lower, int &upper);
  void get_neighbors(const int global_id, int &ieq, int &i, int &j, int nx, int ny, int &left, int &right, int &lower, int &upper);
  void wynn_epsilon(std::vector<double> &sk, int kl, double &sa, int np);

  void pause();
  
}


