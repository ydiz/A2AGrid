#include "../a2a_MADWF.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{

  // Grid_init(&argc,&argv);
  cps::Start(&argc, &argv); // Grid_init(&argc,&argv) is called inside this function 

  JSONReader reader("a2a_madwf.json");

  A2AParams a2a_arg(reader);
  // MADWFParams madwf_arg(reader);
  LanczosParams lanc_arg(reader);
  CGParams_MADWF cg_arg(reader);

  MobiusFermion_arg mob_arg("16I");
  ZMobiusFermion_arg zmob_arg("16I");

  // A2A_MADWF a2a(a2a_arg, madwf_arg, mob_arg, zmob_arg);
  A2A_MADWF a2a(a2a_arg);
  a2a.setFermion(mob_arg, zmob_arg);

  // a2a.load_evecs();
  a2a.load_compressed_evecs("/sdcc/u/dguo/evec/16I/lanczos.output");
  test_evals_evecs(*a2a.zHermOp_f, a2a.evals, a2a.evecs_f);

  a2a.computeVWlow();
  a2a.computeVWhigh(cg_arg);

  Grid_finalize();
}
