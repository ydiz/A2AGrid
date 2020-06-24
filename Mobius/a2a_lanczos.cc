#include "../a2a.h"

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  JSONReader reader("a2a.json");
  LanczosParams lanc_arg(reader);
  A2AParams a2a_arg(reader);
  CGParams cg_arg(reader);

  // ZMobiusFermion_arg fermion_arg("16I");
  MobiusFermion_arg mob_arg("16I");

  A2A a2a(a2a_arg, mob_arg);
  // a2a.setZMobiusFermion(fermion_arg);
  // a2a.setMobiusFermion(fermion_arg);

  a2a.run_lanczos(lanc_arg);
}
