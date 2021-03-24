#include "../a2a.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  JSONReader reader("a2a.json");
  A2AParams a2a_arg(reader);
  CGParams cg_arg(reader);

  MobiusFermion_arg mob_arg("24ID");

  A2A a2a(a2a_arg);
  a2a.setFermion(mob_arg);

  // a2a.run_lanczos(lanc_arg);

  // a2a.load_evecs();
  // test_evals_evecs(*a2a.mobHermOp_f, a2a.evals, a2a.evecs_f);

  // a2a.computeVWlow();
  string save_path = "/hpcgpfs01/work/lqcd/qcdqedta/ydzhao/a2a/24ID/only_high_modes";
  a2a.computeVWhigh(cg_arg, save_path);


  std::cout << "Finished" << std::endl;
  Grid_finalize();
}
