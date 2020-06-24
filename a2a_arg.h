#pragma once

#include <Grid/Grid.h>

namespace Grid{

struct MobiusFermion_arg {
  double mass;
  double M5;
  int Ls;
  double b;

  MobiusFermion_arg(const std::string &ensemble) {
    if(ensemble == "16I") {
      mass = 0.01;
      M5 = 1.8;
      // b = 2.5;
      b = 1.0;  // ? Should be equivalent to Shmir kernel with Ls=16
      Ls = 16;
      // c = 1.5; c = b - 1
    }
    else if(ensemble == "24ID") {
      mass = 0.00107;
      M5 = 1.8;
      b = 2.5; // FIXME: is this right?
      Ls = 24;
      // c = 1.5; c = b - 1
    }
    else assert(0);
  } 
};



struct ZMobiusFermion_arg {
  double mass;
  double M5;
  int Ls_inner;
  std::vector<std::complex<double>> omega;

  ZMobiusFermion_arg(const std::string &ensemble) {
    if(ensemble == "16I") {  // ???? 16I is not Mobius ensemble?? maybe Duo made up those numbers??
      mass = 0.01;
      M5 = 1.8;
      Ls_inner = 12;

      double b_plus_c_inner = 1.;
      double b_plus_c_outer = 1.;
      int Ls_outer = 16;
      double lambda_max = 1.42;
      Grid::Approx::computeZmobiusGamma(omega, b_plus_c_inner, Ls_inner, b_plus_c_outer, Ls_outer, lambda_max);
      std::cout << "omega: " << omega << std::endl;
      // // This is Duo's choice for Mobius b=2.5; 
      // omega.resize(12);
      // omega[0] = std::complex<double>(1.0903256131299373, 0);
      // omega[1] = std::complex<double>(0.9570283702230611, 0);
      // omega[2] = std::complex<double>(0.7048886040934104, 0);
      // omega[3] = std::complex<double>(0.48979921782791747, 0);
      // omega[4] = std::complex<double>(0.328608311201356, 0);
      // omega[5] = std::complex<double>(0.21664245377015995, 0);
      // omega[6] = std::complex<double>(0.14121112711957107, 0);
      // omega[7] = std::complex<double>(0.0907785101745156, 0);
      // omega[8] = std::complex<double>(0.05608303440064219, -0.007537158177840385);
      // omega[9] = std::complex<double>(0.05608303440064219, 0.007537158177840385);
      // omega[10] = std::complex<double>(0.0365221637144842, -0.03343945161367745);
      // omega[11] = std::complex<double>(0.0365221637144842, 0.03343945161367745);
    }
    else if(ensemble == "24ID") {
      mass = 0.00107;
      M5 = 1.8;
      Ls_inner = 12;

      omega.resize(12);
      omega[0] = std::complex<double>(1.0903256131299373, 0);
      omega[1] = std::complex<double>(0.9570283702230611, 0);
      omega[2] = std::complex<double>(0.7048886040934104, 0);
      omega[3] = std::complex<double>(0.48979921782791747, 0);
      omega[4] = std::complex<double>(0.328608311201356, 0);
      omega[5] = std::complex<double>(0.21664245377015995, 0);
      omega[6] = std::complex<double>(0.14121112711957107, 0);
      omega[7] = std::complex<double>(0.0907785101745156, 0);
      omega[8] = std::complex<double>(0.05608303440064219, -0.007537158177840385);
      omega[9] = std::complex<double>(0.05608303440064219, 0.007537158177840385);
      omega[10] = std::complex<double>(0.0365221637144842, -0.03343945161367745);
      omega[11] = std::complex<double>(0.0365221637144842, 0.03343945161367745);
    }
    else assert(0);
  } 
};


struct CGParams : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CGParams,
    double, resid,    // 1.0e-8
    unsigned int, max_iters  // 10000
  );

  template <class ReaderClass>
  CGParams(Reader<ReaderClass>& reader) {
    read(reader, "CG", *this);
    std::cout << *this << std::endl;
  }
};

struct CGParams_MADWF : Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(CGParams_MADWF,
    double, resid,    // 1.0e-8
    double, resid_outer,    
    double, resid_inner,   
    unsigned int, max_iters  // 10000
  );

  template <class ReaderClass>
  CGParams_MADWF(Reader<ReaderClass>& reader) {
    read(reader, "CG_MADWF", *this);
    std::cout << *this << std::endl;
  }
};



struct A2AParams : Serializable {
public:
  std::vector<int> fdims; // lattice size

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AParams,
    int, traj,
    std::string, lat,
    std::string, config,
    int, nhits,
    std::string, ensemble_name,
    std::string, prefix
  );

  template <class ReaderClass>
  A2AParams(Reader<ReaderClass>& reader) {
    read(reader, "A2A", *this);
    GridCmdOptionIntVector(lat, fdims);
    std::cout << "fdims: " << fdims << std::endl;
    std::cout << *this << std::endl;
  }

};

// struct MADWFParams : Serializable {
// public:
//   GRID_SERIALIZABLE_CLASS_MEMBERS(MADWFParams,
//     int, Ls_inner
//   );
//
//   template <class ReaderClass>
//   MADWFParams(Reader<ReaderClass>& reader) {
//     read(reader, "MADWF", *this);
//     std::cout << *this << std::endl;
//   }
//
// };







// struct CG_arg {
//   double resid = 1.0e-8;
//   unsigned int max_iters = 10000;
// };
//
//
// struct A2A_arg {
//   int traj;
//   int Ls;
//   std::vector<int> lat; // lattice size
//   std::string config;
//   int nhits;
//   std::string prefix;
// };
//
// // Chris:
// // double lo = lanc_arg.ch_beta * lanc_arg.ch_beta;
// // double hi = lanc_arg.ch_alpha * lanc_arg.ch_alpha;
// // int ord = lanc_arg.ch_ord + 1; //different conventions
// // Grid::Chebyshev<GridFermionField> Cheb(lo,hi,ord);
// struct Lanczos_arg {
//   // For Chebyshev
//   double ch_lo;
//   double ch_hi;
//   double ch_Npoly; 
//
//   // For Lanczos
//   int Nstop; // N_true_get  // evecs.size() == lanc_arg.N_true_get
//   int Nk; // N_get
//   int Nm; // N_use
//   int MaxIt; //MaxIt
//   double resid;  //resid
// };
// struct Lanczos_arg {
//   // For Chebyshev
//   double ch_beta = std::sqrt(0.2);
//   double ch_alpha = std::sqrt(5.0);
//   double ch_ord = 10; //int ord = lanc_arg.ch_ord + 1; //different conventions
//
//   // For Lanczos
//   int Nstop = 30; // N_true_get  // evecs.size() == lanc_arg.N_true_get
//   int Nk = 40; // N_get
//   int Nm = 80; // N_use
//   int maxits = 10000; //MaxIt
//   double stop_rsd = 1.0e-8;  //resid
//
//   // bool precon = true;
// };
// usage:
// int Nstop;   // Number of evecs checked for convergence
// int Nk;      // Number of converged sought
// int Np;      // Np -- Number of spare vecs in kryloc space
// int Nm;      // Nm -- total number of vectors
/////////////////
// const int Nstop = lanc_arg.N_true_get;
// const int Nk = lanc_arg.N_get;
// const int Np = lanc_arg.N_use - lanc_arg.N_get;
// const int Nm = lanc_arg.N_use;
// const int MaxIt= lanc_arg.maxits;
// Grid::RealD resid = lanc_arg.stop_rsd;



// std::ostream &operator<<(std::ostream& out, const Lanczos_arg & lanc_arg)
// {
//   out << "===========Chebyshev===========" << std::endl;
//   out << "ch_beta: " << lanc_arg.ch_beta << std::endl;
//   out << "ch_alpha: " << lanc_arg.ch_alpha << std::endl;
//   out << "ch_ord; " << lanc_arg.ch_ord << std::endl;
//   out << "===========Lanczos===========" << std::endl;
//   out << "N_true_get: " << lanc_arg.N_true_get << std::endl;
//   out << "N_get: " << lanc_arg.N_get << std::endl;
//   out << "N_use: " << lanc_arg.N_use << std::endl;
//   out << "maxits: " << lanc_arg.maxits << std::endl;
//   out << "stop_rsd: " << lanc_arg.stop_rsd << std::endl;
//   out << "precon: " << lanc_arg.precon << std::endl;
//   return out;
// }


}
