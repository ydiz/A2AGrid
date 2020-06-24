
#pragma once

// #include "mpi.h"
// #include "stdio.h"
// static int comm_counter=0;
//
// int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)
// {
//       int world_rank;
//       MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//       comm_counter++;
//       printf("%s %i %s %i\n", "MPI_Comm_dup ", comm_counter, " from ", world_rank);
//       return PMPI_Comm_dup(comm, newcomm);
// }
//
//
// int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
// {
//       int world_rank;
//       MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//       comm_counter++;
//       printf("%s %i %s %i\n", "MPI_Comm_split ", comm_counter, " from ", world_rank);
//       return PMPI_Comm_split(comm, color, key, newcomm);
// }
//
// int MPI_Comm_free(MPI_Comm *comm)
// {
//       int world_rank;
//       MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//       comm_counter--;
//       printf("%s %i %s %i\n", "PMPI_Comm_free ", comm_counter, " from ", world_rank);
//       return PMPI_Comm_free(comm);
// }


#include <Grid/Grid.h>

#include "io.h"
#include "a2a_arg.h"
#include "a2a_utils.h"
#include "lanczos.h"

#ifdef USE_CPS
#include "read_compressed.h"
#endif

namespace Grid{

//  We add the following method which is necessary for SchurRedBlackDiagTwoSolve: operator()(LinearOperatorBase<Field> &Linop, const Field &src, Field &psi) 
#ifndef USE_CPS
template<class FieldD,class FieldF, typename std::enable_if< getPrecision<FieldD>::value == 2, int>::type = 0,typename std::enable_if< getPrecision<FieldF>::value == 1, int>::type = 0 >
class MixedPrecisionConjugateGradientOp : public MixedPrecisionConjugateGradient<FieldD,FieldF>, public OperatorFunction<FieldD> {
public:
    MixedPrecisionConjugateGradientOp(RealD tol, Integer maxinnerit, Integer maxouterit, 
                                    GridBase* _sp_grid, LinearOperatorBase<FieldF> &_Linop_f, 
                                    LinearOperatorBase<FieldD> &_Linop_d) :
                          MixedPrecisionConjugateGradient<FieldD,FieldF> (tol, maxinnerit, maxouterit, _sp_grid, _Linop_f, _Linop_d) {};

    void operator() (LinearOperatorBase<FieldD> &Linop, const FieldD &in, FieldD &out){
      this->MixedPrecisionConjugateGradient<FieldD,FieldF>::operator()(in,out);
    }
};
#endif




template <class FermOp>
void test_evals_evecs(Grid::SchurDiagTwoOperator<FermOp, LatticeFermionF> &Herm, 
                      std::vector<double> &evals, 
                      std::vector<LatticeFermionF> &evecs) {
  for (int i=0; i<evals.size(); i+=evals.size()/10) {
  // for (int i=0; i<evals.size(); i++) {
      LatticeFermionF tmp(evecs[0].Grid()); 
      Herm.HermOp(evecs[i], tmp);
      double alph = real(innerProduct(evecs[i], tmp));
      tmp = tmp - evals[i] * evecs[i];
      std::cout << GridLogMessage << __func__  << " evecs[" << i << "]: " << " eval: "  << evals[i] <<  " norm2: "  << norm2(evecs[i]) << " <v, Mv>: " << alph << " residual: " << norm2(tmp) << std::endl;
  }
}




class A2A{
public:

  GridCartesian *UGrid;
  GridRedBlackCartesian *UrbGrid;
  GridCartesian *FGrid;
  GridRedBlackCartesian *FrbGrid;
  
  GridCartesian *UGrid_f;
  GridRedBlackCartesian *UrbGrid_f;
  GridCartesian *FGrid_f;
  GridRedBlackCartesian *FrbGrid_f;

  GridParallelRNG URNG;

  int Lt;
  int nl = -1;
  int nhits; // number of hits
  int nh; // number of high modes   nh = nhits * T * 12
  std::string prefix;

  int traj;

  MobiusFermionD *Dmob = NULL; // D_{mobius}, mobius Dirac matrix
  MobiusFermionF *Dmob_f = NULL;
  SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD> *mobHermOp = NULL;
  SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF> *mobHermOp_f = NULL;
  
  LatticeGaugeFieldD Umu; 
  LatticeGaugeFieldF Umu_f; 

  std::vector<double>          evals;
  std::vector<LatticeFermionF>   evecs_f;

  std::vector<LatticeFermionD> vl;
  std::vector<LatticeFermionD> wl;
  std::vector<LatticeFermionD> vh;
  std::vector<LatticeFermionD> wh;

  A2A(const A2AParams &a2a_arg); // MobiusFermion_arg cannot be const, for Mobius fermion constructor
  ~A2A();

  void setFermion(MobiusFermion_arg &arg);

  void run_lanczos(const LanczosParams &lanc_arg);


  void load_evecs();
  #ifdef USE_CPS
  void load_compressed_evecs(const std::string &evec_dir);
  #endif


  void computeVWlow();

  template<typename DiracMatrix>     // Examples of DiracMatrix: MobiusFermionD, ZMobiusFermionD
  void computeVWlow(DiracMatrix &Ddwf);
  
  void computeVWhigh(const CGParams &cg_arg);

  LatticeFermionD getDilutedSource(int idx);

  template<typename DiracMatrix>
  LatticeFermionD compute_lowmode_contrib(DiracMatrix &Ddwf, DeflatedGuesser<LatticeFermionF> &guesser, const LatticeFermionD &eta);

private:
  void highModeIndexToTimeFlavorSpinColorIndex(int idx, int &hit, int &time_idx, int &spin_idx, int &color_idx);
  void extractTimeSpinColor(LatticeFermionD &ret, const LatticeComplexD &lat, int time_idx, int spin_idx, int color_idx);

};







// Umu/URNG/vHigh has no default constructor, thus can not be initialized in constructor body.
A2A::A2A(const A2AParams &a2a_arg) : 
  UGrid(SpaceTimeGrid::makeFourDimGrid(a2a_arg.fdims, GridDefaultSimd(4, vComplexD::Nsimd()), GridDefaultMpi())),
  UrbGrid(SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid)),
  UGrid_f(SpaceTimeGrid::makeFourDimGrid(a2a_arg.fdims, GridDefaultSimd(4, vComplexF::Nsimd()), GridDefaultMpi())),
  UrbGrid_f(SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid_f)),

  URNG(UGrid),

  Umu(UGrid), Umu_f(UGrid_f)
{
  if(a2a_arg.config.empty()) Umu = 1.0;
  // else readGaugeField(Umu, a2a_arg.config);
  else readGaugeField(Umu_f, a2a_arg.config);

  // precisionChange(Umu_f, Umu);
  precisionChange(Umu, Umu_f);

  traj = a2a_arg.traj;
  Lt = a2a_arg.fdims[3];
  nhits = a2a_arg.nhits;
  nh = a2a_arg.nhits * Lt * 12; // * 2; /// ??? Why 2? for Gparity ??
  prefix = a2a_arg.prefix + "/" + a2a_arg.ensemble_name;

  URNG.SeedFixedIntegers({1,2,3,4});
}

// destructor
A2A::~A2A() {
  delete UGrid;
  delete UrbGrid;
  delete FGrid;
  delete FrbGrid;
}

void A2A::setFermion(MobiusFermion_arg &mob_arg) {
  FGrid = SpaceTimeGrid::makeFiveDimGrid(mob_arg.Ls, UGrid);
  FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(mob_arg.Ls, UGrid);
  FGrid_f = SpaceTimeGrid::makeFiveDimGrid(mob_arg.Ls, UGrid_f);
  FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(mob_arg.Ls, UGrid_f);

  // Must set boundary phase in time direction to -1
  typename MobiusFermionD::ImplParams params;
  std::vector<Complex> boundary_phases(4, 1.);
  boundary_phases[3] = -1.;
  params.boundary_phases = boundary_phases;

  Dmob = new MobiusFermionD(Umu, *FGrid, *FrbGrid, *UGrid, *UrbGrid, mob_arg.mass, mob_arg.M5, mob_arg.b, mob_arg.b-1., params);
  Dmob_f = new MobiusFermionF(Umu_f, *FGrid_f, *FrbGrid_f, *UGrid_f, *UrbGrid_f, mob_arg.mass, mob_arg.M5, mob_arg.b, mob_arg.b-1., params);
  mobHermOp = new SchurDiagTwoOperator<MobiusFermionD, LatticeFermionD>(*Dmob);  
  mobHermOp_f = new SchurDiagTwoOperator<MobiusFermionF, LatticeFermionF>(*Dmob_f); 
}


void A2A::run_lanczos(const LanczosParams &lanc_arg) {

  std::vector<LatticeFermionD> evecs;

  lanczos(*Dmob, lanc_arg, evals, evecs);

  nl = evals.size();
  write_evecs_evals_d2f(prefix + "/evals_evecs", evals, evecs);
}




void A2A::load_evecs() {
  read_evecs_evals(prefix + "/evals_evecs", FrbGrid_f, evals, evecs_f);
  for(int i=0; i<evecs_f.size(); ++i) evecs_f[i].Checkerboard() = Odd;

  nl = evals.size();
  std::cout << "Number of eigenvectors: " << nl << std::endl;
  print_memory();
}

#ifdef USE_CPS
void A2A::load_compressed_evecs(const std::string &evec_dir) {

  USING_NAMESPACE_CPS

  typedef A2ApoliciesDoubleAutoAlloc A2Apolicies;
  typedef A2Apolicies::FgridGFclass LatticeType;

  DoArg do_arg;
  if(!do_arg.Decode("do_arg.vml","do_arg")){
    do_arg.Encode("do_arg.templ","do_arg");
    VRB.Result("","main","Can't open do_arg.vml!\n");exit(1);
  }

  int nthreads = omp_get_max_threads();
  std::cout << "Number of threads " << nthreads << std::endl;
  GJP.Initialize(do_arg);  // This sets the number of threads to 1
  GJP.SetNthreads(nthreads);
  std::cout << "Setting number of threads " << nthreads << std::endl;

  //Initialize FGrid
  FgridParams fgp;
  fgp.mobius_scale = 4.0; //b+c //FIXME: read from argument
  LatticeType lattice(fgp);
  lattice.ImportGauge();

  std::cout << GridLogMessage << "Before reading compressed: " << nl << std::endl;
  zyd_read_compressed(evec_dir, FGrid_f, FrbGrid_f, lattice, evals, evecs_f);

  nl = evals.size();
  std::cout << "Number of eigenvectors: " << nl << std::endl;
}
#endif



void A2A::computeVWlow() {
  computeVWlow(*Dmob);

  A2AVectorsIo::write(prefix + "/A2AVector_vl", vl, false, traj);
  A2AVectorsIo::write(prefix + "/A2AVector_wl", wl, false, traj);
  vl.clear();
  wl.clear();
  print_memory();
}

template<typename DiracMatrix>     // Examples of DiracMatrix: MobiusFermionD, ZMobiusFermionD
void A2A::computeVWlow(DiracMatrix &Ddwf) {
  assert(nl>0);
  SchurDiagTwoOperator<DiracMatrix, LatticeFermionD>  HermOp(Ddwf);

  vl.resize(nl, UGrid);
  wl.resize(nl, UGrid);
  for(int i=0; i<nl; ++i) {
    LatticeFermionD evec(FrbGrid);
    precisionChange(evec, evecs_f[i]);      assert(evec.Checkerboard() == Odd);

    // calculate low modes of V; 
    // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i 
    // v_io = (1/eval_i) * MooInv evec_i
    LatticeFermionD vl_5d(FGrid), vl_5d_e(FrbGrid), vl_5d_o(FrbGrid);
    LatticeFermionD tmp(FrbGrid);

    Ddwf.MooeeInv(evec, vl_5d_o);  // The only difference between SchurMooee and SchurTwo: with SchurMooee, this line should be vl_5d_o = evec 
    setCheckerboard(vl_5d, vl_5d_o);        assert(vl_5d_o.Checkerboard() == Odd);

    Ddwf.Meooe(vl_5d_o, tmp);
    Ddwf.MooeeInv(tmp, vl_5d_e);           
    vl_5d_e = - vl_5d_e;
    setCheckerboard(vl_5d, vl_5d_e);        assert(vl_5d_e.Checkerboard() == Even); 

    Ddwf.ExportPhysicalFermionSolution(vl_5d, vl[i]); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)
    vl[i] = (1. / evals[i]) * vl[i];

    // calculate low modes of W
    // w_ie = - MeeInvDag MoeDag Mpc evec_i, where Mpc = D_oo M_oo^{-1}
    // w_io = Mpc evec_i
    // w_i =  D_-^\dagger (w_ie, w_io)
    LatticeFermionD wl_5d(FGrid), wl_5d_o(FrbGrid), wl_5d_e(FrbGrid);

    HermOp.Mpc(evec, wl_5d_o);             
    setCheckerboard(wl_5d, wl_5d_o);        assert(wl_5d_o.Checkerboard() == Odd); 

    Ddwf.MeooeDag(wl_5d_o, tmp);
    Ddwf.MooeeInvDag(tmp, wl_5d_e);
    wl_5d_e = - wl_5d_e;                    
    setCheckerboard(wl_5d, wl_5d_e);        assert(wl_5d_e.Checkerboard() == Even); 
    
    LatticeFermionD tmp_full = wl_5d;
    Ddwf.DminusDag(tmp_full, wl_5d); //Left-multiply by D-^dag. For Cayley fermion D_- = (1-c*DW); for other fermion, D_- = 1

    Ddwf.ExportPhysicalFermionSource(wl_5d, wl[i]); // 5d->4d, w(x) = P_R w(x,0) + P_L w(x, Ls-1)
  }
  std::cout << GridLogMessage << " Finished computing VW_low "<< std::endl;
}





void A2A::computeVWhigh(const CGParams &cg_arg) {

  assert(nl!=-1);
  
  MixedPrecisionConjugateGradientOp<LatticeFermionD, LatticeFermionF> mCG(cg_arg.resid, cg_arg.max_iters, 50, FrbGrid_f, *mobHermOp_f, *mobHermOp);
  DeflatedGuesser<LatticeFermionF> guesser(evecs_f, evals);
  mCG.useGuesser(guesser);   // Must put guesser here, instead of feeding it to solver(xxx);
  SchurRedBlackDiagTwoSolve<LatticeFermionD> solver(mCG);

  vh.resize(nh, UGrid); 
  wh.resize(nh, UGrid);
  for(int i=0; i<nh; ++i) // nh is the number of high modes; nh = nhits * Lt * 12
  {
    std::cout << "High modes: " << i << std::endl;
    // wh is diluted random source
    wh[i] = getDilutedSource(i);

    LatticeFermionD eta(FGrid); // eta is actually D_- eta
    Dmob->ImportPhysicalFermionSource(wh[i], eta); // Project to 5d and multiply it by D_-

    // Calculate vh_5d = D^{-1} \eta
    LatticeFermionD vh_5d(FGrid);
    solver(*Dmob, eta, vh_5d);

    Dmob->ExportPhysicalFermionSolution(vh_5d, vh[i]); // 5d->4d, v(x) = P_L v(x,0) + P_R v(x, Ls-1)

    // remove low mode contribution: (LDU) eta
    LatticeFermionD lowmode_contrib = compute_lowmode_contrib(*Dmob, guesser, eta);
    vh[i] -= lowmode_contrib;  // remove low mode contribution from vh
    
    vh[i] = (1. / nhits) * vh[i]; // including 1/nhit for the hit average

  }
  A2AVectorsIo::write(prefix + "/A2AVector_vh", vh, false, traj);
  A2AVectorsIo::write(prefix + "/A2AVector_wh", wh, false, traj);
}

template<typename DiracMatrix>
LatticeFermionD A2A::compute_lowmode_contrib(DiracMatrix &Ddwf, DeflatedGuesser<LatticeFermionF> &guesser, const LatticeFermionD &eta) {
  // eta must be 5d random source
  // The Ddwf operator must have correct Ls. For MADWF, its Ls should be inner Ls

  ConjugateGradient<LatticeFermionD> dummy_CG(1., 1.);  // This CG is not used; but it is necessary for the constructor of SchurRedBlackDiagTwoSolve
  SchurRedBlackDiagTwoSolve<LatticeFermionD> solver(dummy_CG);

  std::cout << GridLogMessage << "subtracting low modes contribution" << std::endl;
  LatticeFermionD Ueta_e(FrbGrid);         Ueta_e.Checkerboard() = Even; // U * eta
  LatticeFermionD Ueta_o(FrbGrid);         Ueta_o.Checkerboard() = Odd;
  solver.RedBlackSource(Ddwf, eta, Ueta_e, Ueta_o);  

  LatticeFermionF DUeta_o_f(FrbGrid_f);      DUeta_o_f.Checkerboard() = Odd;
  LatticeFermionF Ueta_o_f(FrbGrid_f);       Ueta_o_f.Checkerboard() = Odd;

  // // converting eigenvectors to double is very slow; so I convert source to float. This is also what mixed CG does internally
  precisionChange(Ueta_o_f, Ueta_o);
  guesser(Ueta_o_f, DUeta_o_f);

  LatticeFermionD DUeta_o(FrbGrid);
  precisionChange(DUeta_o, DUeta_o_f);

  LatticeFermionD DUeta_e(FrbGrid);      DUeta_e.Checkerboard() = Even;
  DUeta_e = Zero();

  LatticeFermionD LDUeta(FGrid);
  solver.RedBlackSolution(Ddwf, DUeta_o, DUeta_e, LDUeta);

  LatticeFermionD LDUeta_4d(UGrid);
  Ddwf.ExportPhysicalFermionSolution(LDUeta, LDUeta_4d); 

  std::cout << GridLogMessage << "Finished subtracting low modes contribution" << std::endl;
  return LDUeta_4d;
}



// eta(x) = exp( 2 PI I t), where t ~ uniform(0, 1)
LatticeFermionD A2A::getDilutedSource(int idx) {
  static bool src_created = false;
  static std::vector<LatticeComplexD> rand;
  if(!src_created) {
    rand.assign(nhits, UGrid);
    for(int i=0; i<nhits; ++i) {     // use the same random number for every hit
      LatticeRealD rand_real(UGrid);     // tmp_real ~ uniform(0, 1)
      random(URNG, rand_real);
      rand[i] = exp(2 * M_PI * timesI(toComplex(rand_real)));
    }
    src_created = true;
  }

  LatticeFermionD src(UGrid);
  int hit, t, s, c;
  highModeIndexToTimeFlavorSpinColorIndex(idx, hit, t, s, c);
  extractTimeSpinColor(src, rand[hit], t, s, c);

  return src;
}

void A2A::highModeIndexToTimeFlavorSpinColorIndex(int idx, int &hit, int &time_idx, int &spin_idx, int &color_idx)
{
  int T = UGrid->FullDimensions()[3];
  int ncolor = 3;
  int nspin = 4;
  color_idx = idx % ncolor; idx /= ncolor;
  spin_idx = idx % nspin; idx /= nspin;
  time_idx = idx % Lt; idx /= T;
  hit = idx;
}


// for dilution
// ret(x,t)|_{t=t0, f=f0,s=s0,c=c0)} = lat(x,t)|_{t=t0}, value only depends on time_idx, not spin and color index
void A2A::extractTimeSpinColor(LatticeFermionD &ret, const LatticeComplexD &lat, int time_idx, int spin_idx, int color_idx)
{
  assert(ret.Grid()->Nd()==4 && lat.Grid()->Nd()==4); // must be 4d fermions

  ret = Zero();  // set ret to 0
  thread_for(ss, lat.Grid()->lSites(), {
    Coordinate lcoor, gcoor;
    localIndexToLocalGlobalCoor(lat.Grid(), ss, lcoor, gcoor);
    if(gcoor[3] == time_idx) {
      typename LatticeFermionD::vector_object::scalar_object ret_site;
      typename LatticeComplexD::vector_object::scalar_object lat_site;  
      peekLocalSite(lat_site, lat, lcoor); // zyd: peekLocalSite is efficient; peekSite is not
      ret_site()(spin_idx)(color_idx) = lat_site()()();
      pokeLocalSite(ret_site, ret, lcoor);
    }
  });
}



} // namespace Grid
