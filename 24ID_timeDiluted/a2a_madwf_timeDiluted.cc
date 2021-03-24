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


  JSONReader reader("a2a_madwf_timeDiluted.json");

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

  a2a.computeVWlow();    // FIXME: uncomment this line
  a2a.computeVWhigh(cg_arg);

  Grid_finalize();
}
