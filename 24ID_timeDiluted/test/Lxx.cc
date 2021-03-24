#include "../a2a.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc, &argv);

  Coordinate gcoor({24, 24, 24, 64});
  int T = gcoor[3];
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  std::string prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID";
  int nl = 2000;
  int nh = 768;
  int traj = 2300;
  // /hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a_ml10//24ID/MADWF_A2AVector_wh_no_timeDilution.2300.bin

  std::vector<LatticeFermionD> v(nl, UGrid), w(nl, UGrid), vh(nh, UGrid), wh(nh, UGrid);
  A2AVectorsIo::read(v, prefix + "/MADWF_A2AVector_vl", false, traj);
  A2AVectorsIo::read(w, prefix + "/MADWF_A2AVector_wl", false, traj);

  A2AVectorsIo::read(vh, prefix + "/MADWF_A2AVector_vh_no_timeDilution", false, traj);
  v.insert(v.end(), vh.begin(), vh.end());          vh.clear();
  print_memory();

  A2AVectorsIo::read(wh, prefix + "/MADWF_A2AVector_wh_no_timeDilution", false, traj);
  w.insert(w.end(), wh.begin(), wh.end());          wh.clear();
  print_memory();

  double vol = 1.;
  for(int i=0; i<4; ++i) vol *= gcoor[i];

  LatticePropagatorD Lxx(UGrid);
  Lxx = Zero();

  for(int i=0; i<nl + nh; ++i) {
    if(i%100==0) std::cout << GridLogMessage << "i = " << i << std::endl;
    Lxx = Lxx + outerProduct(v[i], w[i]);  // v[i](x) w[i]^\dagger(x)


    if(i==nl-1) {
      LatticeComplexD tr_Lxx = trace(Lxx);
      std::cout << "After low modes: quark condensate: " << sum(tr_Lxx) / vol / 12. << std::endl;
    }
  }
  writeScidac(Lxx, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/Lxx/Lxx." + to_string(traj));


  LatticeComplexD tr_Lxx = trace(Lxx);

  using LatticeComplexSite = typename LatticeComplexD::vector_object::scalar_object;
  std::vector<LatticeComplexSite> tr_Lxx_slice_sum;
  sliceSum(tr_Lxx, tr_Lxx_slice_sum, Tdir);
  std::cout << "Time slice sum of tr(L(x,x)): " << tr_Lxx_slice_sum << std::endl;

  std::cout << "quark condensate: " << sum(tr_Lxx) / vol / 12. << std::endl;

  Grid_finalize();

  // ComplexD rst = 0.;
  // for(int i=0; i<nl; ++i) {
  //   rst += innerProduct(wl[i], vl[i]);   // wl[i]^dagger vl[i]
  //   std::cout << "rst per site: " << rst / vol << std::endl;
  // }
  // std::cout << "rst from low modes: " << rst << std::endl;
  // std::cout << "rst per site: " << rst / vol << std::endl;
  //
  // for(int i=0; i<nh; ++i) {
  //   rst += innerProduct(wh[i], vh[i]);   // wh[i]^dagger vh[i]
  //   std::cout << "rst per site: " << rst / vol << std::endl;
  // }
  // std::cout << "rst after adding high modes: " << rst << std::endl;
  // std::cout << "rst per site: " << rst / vol << std::endl;



}
