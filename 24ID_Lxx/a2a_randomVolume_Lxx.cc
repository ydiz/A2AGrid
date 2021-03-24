#include "../a2a_MADWF.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{

  // Grid_init(&argc,&argv);
  cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function 

  int traj;
  if( GridCmdOptionExists(argv, argv+argc, "--traj") ) {
    string arg = GridCmdOptionPayload(argv, argv+argc, "--traj");
    GridCmdOptionInt(arg, traj);
  }
  else {
    std::cout << "traj not specified; exiting" << std::endl;
    assert(0);
  }

  JSONReader reader("a2a_madwf_notTimeDiluted.json");

  A2AParams a2a_arg(traj, reader);
  LanczosParams lanc_arg(reader);
  CGParams_MADWF cg_arg(reader);

  MobiusFermion_arg mob_arg("24ID");
  ZMobiusFermion_arg zmob_arg("24ID");

  A2A_MADWF a2a(a2a_arg);
  a2a.setFermion(mob_arg, zmob_arg);

  // a2a.load_compressed_evecs("None");
  // a2a.load_compressed_evecs("/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/compressed_eigenvector/2300/lanczos.output");
  a2a.load_compressed_evecs(a2a_arg.evec_path);
  test_evals_evecs(*a2a.zHermOp_f, a2a.evals, a2a.evecs_f);

  a2a.computeVWlow();



  LatticePropagatorD Lxx(a2a.UGrid);  Lxx = Zero();
  for(int i=0; i<a2a.nl; ++i) Lxx += outerProduct(a2a.vl[i], a2a.wl[i]);
  a2a.vl.clear();
  a2a.wl.clear();
  // print_memory();


  a2a.computeVWhigh(cg_arg);

  for(int i=0; i<a2a.nh; ++i) Lxx += outerProduct(a2a.vh[i], a2a.wh[i]);

///////////////////////////////////////////////////////////////




  // for(int i=0; i<nl + nh; ++i) {
  //   if(i%100==0) std::cout << GridLogMessage << "i = " << i << std::endl;
  //   Lxx = Lxx + outerProduct(v[i], w[i]);  // v[i](x) w[i]^\dagger(x)
  //
  //
  //   // if(i==nl-1) {
  //   //   LatticeComplexD tr_Lxx = trace(Lxx);
  //   //   std::cout << "After low modes: quark condensate: " << sum(tr_Lxx) / vol / 12. << std::endl;
  //   // }
  // }
  writeScidac(Lxx, "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/Lxx/Lxx." + to_string(traj));



  // using LatticeComplexSite = typename LatticeComplexD::vector_object::scalar_object;
  // std::vector<LatticeComplexSite> tr_Lxx_slice_sum;
  // sliceSum(tr_Lxx, tr_Lxx_slice_sum, Tdir);
  // std::cout << "Time slice sum of tr(L(x,x)): " << tr_Lxx_slice_sum << std::endl;

  double vol = 1.;
  for(int i=0; i<4; ++i) vol *= a2a_arg.fdims[i];
  LatticeComplexD tr_Lxx = trace(Lxx);
  std::cout << "quark condensate: " << sum(tr_Lxx) / vol / 12.  << " ( should be 0.002134(3) )"<< std::endl;
  // TBD: print true value


  Grid_finalize();
}
