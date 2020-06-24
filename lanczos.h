#pragma once

#include <Grid/Grid.h>

namespace Grid {

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
    std::cout << *this << std::endl;
  }
};

template<typename DiracMatrix>     // Examples of DiracMatrix: MobiusFermionD, ZMobiusFermionD
void lanczos(DiracMatrix &Ddwf, const LanczosParams &lanc_arg, std::vector<double> &evals, std::vector<LatticeFermionD> &evecs) {

  GridBase *grid = Ddwf.FermionRedBlackGrid();

  SchurDiagTwoOperator<DiracMatrix, LatticeFermionD>  HermOp(Ddwf);
  Chebyshev<LatticeFermionD> Cheby(lanc_arg.Cheby);
  FunctionHermOp<LatticeFermionD> OpCheby(Cheby, HermOp);
  PlainHermOp<LatticeFermionD> Op(HermOp);

  ImplicitlyRestartedLanczos<LatticeFermionD> IRL(OpCheby, Op, lanc_arg.Nstop, lanc_arg.Nk, lanc_arg.Nm, lanc_arg.resid, lanc_arg.MaxIt);

  LatticeFermionD    src(grid);
  {
    src = 1.0;
    src.Checkerboard() = Odd; // Must be Odd
    // normalize
    RealD nn = norm2(src);
    nn = Grid::sqrt(nn);
    src = src * (1.0 / nn);
  }

  evals.resize(lanc_arg.Nm); // there is an assert `Nm <= evec.size() && Nm <= eval.size()'
  evecs.resize(lanc_arg.Nm, grid);
  for(int i=0; i<evecs.size(); i++) evecs[i].Checkerboard() = Grid::Odd;

  int Nconv;
  IRL.calc(evals, evecs, src, Nconv);
}




}
