#include <Grid/Grid.h>
#include "a2a_arg.h"



namespace Grid {

struct CGParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CGParams,
    double, resid,    // 1.0e-8
    unsigned int, max_iters  // 10000
  );

  template <class ReaderClass>
  CGParams(Reader<ReaderClass>& reader) {
    read(reader, "CG", *this);
  }
};

struct A2AParams : Serializable {
public:
  std::vector<int> fdims; // lattice size

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AParams,
    int, traj,
    int, Ls,
    std::string, lat,
    std::string, config,
    int, nhits,
    std::string, prefix
  );

  template <class ReaderClass>
  A2AParams(Reader<ReaderClass>& reader) {
    read(reader, "A2A", *this);
    GridCmdOptionIntVector(lat, fdims);
    std::cout << "fdims: " << fdims << std::endl;
  }

};

struct LanczosParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LanczosParams,
    ChebyParams, Cheby,/*Chebyshev*/
    int, Nstop,    /*Vecs in Lanczos must converge Nstop < Nk < Nm*/
    int, Nk,       /*Vecs in Lanczos seek converge*/
    int, Nm,       /*Total vecs in Lanczos include restart*/
    RealD, resid,  /*residual*/
    int, MaxIt 
  );

  template <class ReaderClass>
  LanczosParams(Reader<ReaderClass>& reader) {
    read(reader, "Lanczos", *this);
  }
};

}

using namespace Grid;
using namespace std;

int main(int argc, char **argv)
{
  Grid_init(&argc, &argv);

  JSONReader reader("a2a.json");
  LanczosParams lanc_arg(reader);
  A2AParams a2a_arg(reader);
  CGParams cg_arg(reader);

  std::cout << a2a_arg << std::endl;
  std::cout << lanc_arg << std::endl;
  std::cout << cg_arg << std::endl;
}
