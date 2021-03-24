#pragma once

#include "a2a.h"

namespace Grid {

// // FIXME: figure this out
// template<typename FermionFieldType>
// struct CGincreaseTol : public Grid::MADWFinnerIterCallbackBase{
//  Grid::ConjugateGradient<FermionFieldType> &cg_inner;  
//  Grid::RealD outer_resid;
//
//  CGincreaseTol(Grid::ConjugateGradient<FermionFieldType> &cg_inner,
//           Grid::RealD outer_resid): cg_inner(cg_inner), outer_resid(outer_resid){}
//  
//  void operator()(const Grid::RealD current_resid){
//    std::cout << "CGincreaseTol with current residual " << current_resid << " changing inner tolerance " << cg_inner.Tolerance << " -> ";
//    while(cg_inner.Tolerance < current_resid) cg_inner.Tolerance *= 2;    
//    //cg_inner.Tolerance = outer_resid/current_resid;
//    std::cout << cg_inner.Tolerance << std::endl;
//  }
// };


class A2A_MADWF : public A2A {
public:
  
  A2A_MADWF(const A2AParams &a2a_arg) : A2A(a2a_arg) {}
  
  ZMobiusFermionD *Dzmob;
  ZMobiusFermionF *Dzmob_f;
  SchurDiagTwoOperator<ZMobiusFermionD, LatticeFermionD> *zHermOp = NULL;
  SchurDiagTwoOperator<ZMobiusFermionF, LatticeFermionF> *zHermOp_f = NULL;

  GridCartesian *mobFGrid;
  GridRedBlackCartesian *mobFrbGrid;
  GridCartesian *mobFGrid_f;
  GridRedBlackCartesian *mobFrbGrid_f;

  void setFermion(MobiusFermion_arg &mob_arg, ZMobiusFermion_arg &zmob_arg);

  void run_lanczos(const LanczosParams &lanc_arg);
  void load_evecs();
  void computeVWlow();
  void computeVWhigh(const CGParams_MADWF &cg_arg);
};



void A2A_MADWF::setFermion(MobiusFermion_arg &mob_arg, ZMobiusFermion_arg &zmob_arg) {
  // For inner ZMobius
  FGrid = SpaceTimeGrid::makeFiveDimGrid(zmob_arg.Ls_inner, UGrid);
  FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(zmob_arg.Ls_inner, UGrid);
  FGrid_f = SpaceTimeGrid::makeFiveDimGrid(zmob_arg.Ls_inner, UGrid_f);
  FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(zmob_arg.Ls_inner, UGrid_f);

  // For outer Mobius
  mobFGrid = SpaceTimeGrid::makeFiveDimGrid(mob_arg.Ls, UGrid);
  mobFrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(mob_arg.Ls, UGrid);
  mobFGrid_f = SpaceTimeGrid::makeFiveDimGrid(mob_arg.Ls, UGrid_f);
  mobFrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(mob_arg.Ls, UGrid_f);

  typename ZMobiusFermionD::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);
  boundary_phases[3] = -1.;
  params.boundary_phases = boundary_phases;

  Dmob = new MobiusFermionD(Umu, *mobFGrid, *mobFrbGrid, *UGrid, *UrbGrid, mob_arg.mass, mob_arg.M5, mob_arg.b, mob_arg.b-1., params);
  Dmob_f = new MobiusFermionF(Umu_f, *mobFGrid_f, *mobFrbGrid_f, *UGrid_f, *UrbGrid_f, mob_arg.mass, mob_arg.M5, mob_arg.b, mob_arg.b-1., params);
  mobHermOp = new SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD>(*Dmob);  
  mobHermOp_f = new SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF>(*Dmob_f); 

  Dzmob = new ZMobiusFermionD(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, zmob_arg.mass, zmob_arg.M5, zmob_arg.omega, 1., 0., params);
  Dzmob_f = new ZMobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, zmob_arg.mass, zmob_arg.M5, zmob_arg.omega, 1., 0., params);
  zHermOp = new SchurDiagTwoOperator<ZMobiusFermionD, LatticeFermionD>(*Dzmob);  
  zHermOp_f = new SchurDiagTwoOperator<ZMobiusFermionF, LatticeFermionF>(*Dzmob_f); 

}

void A2A_MADWF::run_lanczos(const LanczosParams &lanc_arg) {

  std::vector<LatticeFermionD> evecs;

  lanczos(*Dzmob, lanc_arg, evals, evecs);

  nl = evals.size();
  write_evecs_evals_d2f(prefix + "/evals_evecs_Zmobius", evals, evecs);
}

void A2A_MADWF::load_evecs() {
  read_evecs_evals(prefix + "/evals_evecs_Zmobius", FrbGrid_f, evals, evecs_f);
  for(int i=0; i<evecs_f.size(); ++i) evecs_f[i].Checkerboard() = Odd;

  nl = evals.size();
  std::cout << "Number of eigenvectors: " << nl << std::endl;
  print_memory();
}

void A2A_MADWF::computeVWlow() {

  A2A::computeVWlow(*Dzmob); // call method of parent class, which has been "name hiding"
  A2AVectorsIo::write(prefix + "/MADWF_A2AVector_vl", vl, false, traj);
  A2AVectorsIo::write(prefix + "/MADWF_A2AVector_wl", wl, false, traj);

  vl.clear();
  wl.clear();
  print_memory();
}


void A2A_MADWF::computeVWhigh(const CGParams_MADWF &cg_arg) {

   assert(nl!=-1);
   // Inner CG (only appear in PV): Mixed precision CG; 
   // outer CG: single precision CG

   // CG_outer is used only for PV
   // ConjugateGradient<LatticeFermionD> CG_outer(cg_arg.resid_outer, cg_arg.max_iters); 
   using MixedCG = MixedPrecisionConjugateGradientOp<LatticeFermionD, LatticeFermionF>;
   MixedCG CG_outer(cg_arg.resid_outer, cg_arg.max_iters, 50, mobFrbGrid_f, *mobHermOp_f, *mobHermOp);
   // using PVtype = PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD>;
   using PVtype = PauliVillarsSolverFourierAccel<LatticeFermionD, LatticeGaugeFieldD, MixedCG>;
   PVtype PV_outer(Umu, CG_outer);

   // //Setup update control
   // CGincreaseTol<LatticeFermionF> tol_control_inner(CG_inner, cg_arg.resid);
   // MADWF<MobiusFermionD, ZMobiusFermionF, PVtype, SchurSolverF, GuessF> madwf(*Dmob, *Dzmob_f, PV_outer, SchurSolver_inner, guesser_inner, cg_arg.resid, 100, &tol_control_inner);

   ConjugateGradient<LatticeFermionF> CG_inner(cg_arg.resid_inner, cg_arg.max_iters, 0);
   using SchurSolverF = SchurRedBlackDiagTwoSolve<LatticeFermionF>;
   SchurSolverF SchurSolver_inner(CG_inner);
        
   using GuessF = DeflatedGuesser<LatticeFermionF>;
   GuessF guesser_inner(evecs_f, evals);

   MADWF<MobiusFermionD, ZMobiusFermionF, PVtype, SchurSolverF, GuessF> madwf(*Dmob, *Dzmob_f, PV_outer, SchurSolver_inner, guesser_inner, cg_arg.resid, 100);

  vh.resize(nh, UGrid); 
  wh.resize(nh, UGrid);
  print_memory();

  for(int i=0; i<nh; ++i) // nh is the number of high modes; nh = nhits * Lt * 12
  {
    std::cout << "High modes: " << i << std::endl;
    // wh is diluted random source
    wh[i] = getDilutedSource(i);

    LatticeFermionD eta(Dmob->FermionGrid()); // eta is actually D_- eta 
    Dmob->ImportPhysicalFermionSource(wh[i], eta); // Project to 5d and multiply it by D_-

    // Calculate vh_5d = D^{-1} \eta
    LatticeFermionD vh_5d(Dmob->FermionGrid());
    // solver(*Dmob, eta, vh_5d);
    madwf(wh[i], vh_5d);   // Input must be src4d, not src5d// Change 1: solver -> madwf

    Dmob->ExportPhysicalFermionSolution(vh_5d, vh[i]); 

    // remove low mode contribution: (LDU) eta
    LatticeFermionD lowmode_contrib = compute_lowmode_contrib(*Dzmob, guesser_inner, eta); // Change 2: Dmob -> Dzmob
    vh[i] -= lowmode_contrib;  // remove low mode contribution from vh
    
    vh[i] = (1. / nhits) * vh[i]; // including 1/nhit for the hit average

  }
  A2AVectorsIo::write(prefix + "/MADWF_A2AVector_vh", vh, false, traj);
  A2AVectorsIo::write(prefix + "/MADWF_A2AVector_wh", wh, false, traj);
}


}
