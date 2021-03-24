#include "../a2a_MADWF.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{

  // Grid_init(&argc,&argv);
  cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function 

  JSONReader reader("a2a_madwf.json");

  A2AParams a2a_arg(reader);
  LanczosParams lanc_arg(reader);
  CGParams_MADWF cg_arg(reader);

  MobiusFermion_arg mob_arg("24ID");
  ZMobiusFermion_arg zmob_arg("24ID");

  A2A_MADWF a2a(a2a_arg);
  a2a.setFermion(mob_arg, zmob_arg);

  // a2a.load_compressed_evecs("None");
  // a2a.load_compressed_evecs("/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/compressed_eigenvector/2300/lanczos.output");
  // test_evals_evecs(*a2a.zHermOp_f, a2a.evals, a2a.evecs_f);

  // a2a.computeVWlow();
  string save_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/only_high_modes";
  a2a.computeVWhigh(cg_arg, save_path);

  Grid_finalize();
}
