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
      // mass = 10;   // FIXME: for testing quark condensate
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
  // std::vector<int> fdims; // lattice size
    std::string config;
    std::string evec_path;
    int traj;

public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AParams,
    // int, traj,
    // std::string, lat,
    std::vector<int>, fdims, // lattice size
    std::string, config_prefix,
    std::string, evec_prefix,
    int, nhits,
    std::string, ensemble_name,
    std::string, prefix,
    bool, doTimeDilution
  );

  template <class ReaderClass>
  A2AParams(int _traj, Reader<ReaderClass>& reader) {
    traj = _traj;

    read(reader, "A2A", *this);

    // GridCmdOptionIntVector(lat, fdims);
    config = config_prefix + "/ckpoint_lat." + std::to_string(traj);
    evec_path = evec_prefix + "/" + std::to_string(traj) + "/lanczos.output";

    std::cout << "trajectory: " << traj  << std::endl;
    std::cout << "config path: " << config << std::endl;
    std::cout << "evec path: " << evec_path << std::endl;

    // std::cout << "fdims: " << fdims << std::endl;
    std::cout << *this << std::endl;
  }

};

}
