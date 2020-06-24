#include <Grid/Grid.h>
#include "io.h"

using namespace std;
using namespace Grid;

// typedef typename DomainWallFermionR::FermionField FermionField;
// // typedef LatticeFermion FermionField;

struct Params {
  double mass;
  double M5;
  std::vector<std::complex<double>> omega;

  Params(const std::string &ensemble) {
    if(ensemble == "16I") {
      mass = 0.01;
      M5 = 1.8;

      omega.resize(12);
      omega[0] = std::complex<double>(1.0903256131299373, 0);
      omega[1] = std::complex<double>(0.9570283702230611, 0);
      omega[2] = std::complex<double>(0.7048886040934104, 0);
      omega[3] = std::complex<double>(0.48979921782791747, 0);
      omega[4] = std::complex<double>(0.328608311201356, 0);
      omega[5] = std::complex<double>(0.21664245377015995, 0);
      omega[6] = std::complex<double>(0.14121112711957107, 0);
      omega[7] = std::complex<double>(0.0907785101745156, 0);
      omega[8] = std::complex<double>(0.05608303440064219, -0.007537158177840385);
      omega[9] = std::complex<double>(0.05608303440064219, 0.007537158177840385);
      omega[10] = std::complex<double>(0.0365221637144842, -0.03343945161367745);
      omega[11] = std::complex<double>(0.0365221637144842, 0.03343945161367745);
    }
    else assert(0);
  }

}; 



int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate gcoor({16,16,16,32});
  const int Ls = 24;
  std::cout << "Ls: " << Ls << std::endl;

  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
  printf("UGrid=%p UrbGrid=%p FGrid=%p FrbGrid=%p\n", UGrid, UrbGrid, FGrid, FrbGrid);

  std::string config = "/sdcc/u/dguo/config/16I/ckpoint_lat.2000";
  LatticeGaugeField Umu(UGrid);
  readGF(Umu, config);

  typename MobiusFermionR::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);
  boundary_phases[Nd-1] = -1.;
  params.boundary_phases = boundary_phases;

  Params ensem("16I");
  double mob_b = 2.5;
  MobiusFermionR  Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, ensem.mass, ensem.M5, mob_b, mob_b-1., params);
  // ZMobiusFermionR Ddwf(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, ensem.mass, ensem.M5, ensem.omega, 1., 0., params);
  SchurDiagTwoOperator<MobiusFermionR, LatticeFermion> HermOp(Ddwf);
 
  double alpha = 0.002;
  // double beta = 35;
  double beta = 5.5;
  int Npoly = 101;
  Chebyshev<LatticeFermion> Cheby(alpha, beta, Npoly);

  FunctionHermOp<LatticeFermion> OpCheby(Cheby,HermOp);
  PlainHermOp<LatticeFermion> Op(HermOp);

  const int Nstop = 30;
  const int Nk = 40;
  const int Np = 40;
  const int Nm = Nk+Np;
  const int MaxIt= 10000;
  // RealD resid = 1.0e-8;
  RealD resid = 9e-6;
  ImplicitlyRestartedLanczos<LatticeFermion> IRL(OpCheby,Op,Nstop,Nk,Nm,resid,MaxIt);
  
  std::vector<RealD>          eval(Nm);
  std::vector<LatticeFermion> evec(Nm,FrbGrid);

  // Construct source vector
  LatticeFermion    src(FrbGrid);
  {
    src=1.0;
    src.Checkerboard() = Odd;

    // normalize
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0/nn);
  }

  int Nconv;
  IRL.calc(eval,evec,src,Nconv);


  Grid_finalize();
}
