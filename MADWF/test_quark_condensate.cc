#include "../a2a.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate gcoor({16, 16, 16, 32});
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  std::string prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a";
  int nl = 71;
  int nh = 384;
  int traj = 2000;

  double vol = 1.;
  for(int i=0; i<4; ++i) vol *= gcoor[i];

  std::vector<LatticeFermionD> vl(nl, UGrid), wl(nl, UGrid), vh(nh, UGrid), wh(nh, UGrid);
  A2AVectorsIo::read(vl, prefix + "/MADWF_A2AVector_vl", true, traj);
  A2AVectorsIo::read(wl, prefix + "/MADWF_A2AVector_wl", true, traj);
  A2AVectorsIo::read(vh, prefix + "/MADWF_A2AVector_vh", true, traj);
  A2AVectorsIo::read(wh, prefix + "/MADWF_A2AVector_wh", true, traj);

  ComplexD rst = 0.;
  for(int i=0; i<nl; ++i) {
    rst += innerProduct(wl[i], vl[i]);   // wl[i]^dagger vl[i]
    std::cout << "rst per site: " << rst / vol << std::endl;
  }
  std::cout << "rst from low modes: " << rst << std::endl;
  std::cout << "rst per site: " << rst / vol << std::endl;

  for(int i=0; i<nh; ++i) {
    rst += innerProduct(wh[i], vh[i]);   // wh[i]^dagger vh[i]
    std::cout << "rst per site: " << rst / vol << std::endl;
  }
  std::cout << "rst after adding high modes: " << rst << std::endl;
  std::cout << "rst per site: " << rst / vol << std::endl;


  Grid_finalize();

}
