#include <Grid/Grid.h>

using namespace Grid;
using namespace std;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);


  double b_plus_c_inner=1.;
  int Ls_inner=12;
  // double b_plus_c_outer=4.;
  double b_plus_c_outer=1.;
  int Ls_outer=16;
  // int Ls_outer=24;
  
  double lambda_max = 1.42;


  // void computeZmobiusGamma(std::vector<ComplexD> &gamma_out, 
  //      const RealD mobius_param_out, const int Ls_out, 
  //      const RealD mobius_param_in, const int Ls_in,
  //      const RealD lambda_bound);

  std::vector<Grid::ComplexD> gamma_inner;
  Grid::Approx::computeZmobiusGamma(gamma_inner, b_plus_c_inner, Ls_inner, b_plus_c_outer, Ls_outer, lambda_max);
  std::cout << gamma_inner << std::endl;

  // std::cout << "====================================" << std::endl;
  //
  // Grid::Approx::zolotarev_data *zdata = Grid::Approx::higham(1.0,Ls_inner);
  // gamma_inner.resize(Ls_inner);
  // for(int s=0;s<Ls_inner;s++) gamma_inner[s] = zdata->gamma[s];
  // std::cout << gamma_inner << std::endl;
}
