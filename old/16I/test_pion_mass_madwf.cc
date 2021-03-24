#include "../a2a.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate gcoor({16, 16, 16, 32});
  int T = gcoor[3];
  GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(gcoor, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi());

  std::string prefix = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/16I";
  int nl = 70;
  int nh = 384;
  int traj = 2000;

  double vol = 1.;
  for(int i=0; i<4; ++i) vol *= gcoor[i];

  std::vector<LatticeFermionD> v(nl, UGrid), w(nl, UGrid), vh(nh, UGrid), wh(nh, UGrid);
  A2AVectorsIo::read(v, prefix + "/MADWF_A2AVector_vl", false, traj);
  A2AVectorsIo::read(w, prefix + "/MADWF_A2AVector_wl", false, traj);

  A2AVectorsIo::read(vh, prefix + "/MADWF_A2AVector_vh", false, traj);
  v.insert(v.end(), vh.begin(), vh.end());          vh.clear();
  print_memory();

  A2AVectorsIo::read(wh, prefix + "/MADWF_A2AVector_wh", false, traj);
  w.insert(w.end(), wh.begin(), wh.end());          wh.clear();
  print_memory();

  Gamma g5(Gamma::Algebra::Gamma5);

  std::cout << GridLogMessage << "Calculating meson field" << std::endl;
  vector<vector<vector<ComplexD>>> mf(nl+nh); // mesonfield(i, j)[t]

  for(auto &x: mf) x.resize(nl+nh, vector<ComplexD>(T));

  for(int i=0; i<nl+nh; ++i) 
    for(int j=0; j<nl+nh; ++j) 
      sliceInnerProductVector(mf[i][j], LatticeFermionD(w[i]), LatticeFermionD(g5 * v[j]), Tdir); // For intel compiler, need to expplicitly convert w[i] to LatticeFermionD; do not know why


  std::cout << GridLogMessage << "Calculating pion correlator" << std::endl;

  vector<vector<ComplexD>> rst(T, vector<ComplexD>(T)); // rst(t, tp)
  std::cout << rst.size() << std::endl;
  std::cout << rst[0].size() << std::endl;
  std::cout << mf.size() << std::endl;
  std::cout << mf[0].size() << std::endl;
  std::cout << mf[0][0].size() << std::endl;

  for(int t=0; t<T; ++t) {
    for(int tp=0; tp<T; ++tp) {
    std::cout << GridLogMessage << t << " " << tp << std::endl;
      rst[t][tp] = 0.;
      for(int i=0; i<nl+nh; ++i) {
        for(int j=0; j<nl+nh; ++j) {
          rst[t][tp] += mf[i][j][t] *  mf[j][i][tp];
        }
      }
    }
  }
  std::cout << rst << std::endl;

  std::cout << "Finished! " << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);
  Grid_finalize();

}
