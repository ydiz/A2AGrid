#pragma once

#include <Grid/Grid.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdlib>


// #include "/sdcc/u/ydzhao/decays/headers/io.h"

// bool dirExists(const std::string &path){
//   struct stat info;
//   if( stat( path.c_str(), &info ) == 0 ) return true; // dir does exist
//   else return false;
// }


namespace Grid {

template<class T>
void print_grid_field_site(const T &field, const std::vector<int> coor) {
  using namespace Grid;
  std::cout << "[ " << coor[0] << " " << coor[1] << " " << coor[2] << " " << coor[3] << " ]" << std::endl;
  typename T::vector_object::scalar_object site;
  peekSite(site, field, coor);
  std::cout << site << std::endl;
}


bool dirExists(const std::string &path){
  struct stat info;
  if( stat( path.c_str(), &info ) == 0 ) return true; // dir does exist
  else return false;
}


void readGF(LatticeGaugeField &U, const std::string &filename)
{
  FieldMetaData header;
  NerscIO::readConfiguration(U,header,filename);
}

// template<class T>
// void writeScidac(T& field, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
//   if(field.Grid()->IsBoss()) {
//     std::string base_dir = filename.substr(0, filename.rfind('/'));
//     std::cout << "base_dir: " << base_dir << std::endl;
//     assert(dirExists(base_dir));
//     system(("mkdir -p " + base_dir).c_str());      //  sometimes has error
//   }
//   MPI_Barrier(MPI_COMM_WORLD);
//
//   emptyUserRecord record;
//   ScidacWriter WR(field.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
//   WR.open(filename);
//   WR.writeScidacFieldRecord(field, record);
//   WR.close();
// };
//
// template<class T>
// void readScidac(T& field, const std::string &filename){
//   emptyUserRecord record;
//   ScidacReader RD;
//   RD.open(filename);
//   RD.readScidacFieldRecord(field, record);
//   RD.close();
// };




/////////////////////////////////////
// Taken from  Hadrons/Hadrons/A2AVectors.hpp, Hadrons/Hadrons/Global.cpp 
// usage: A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
/////////////////////////////////////

#define MAX_PATH_LENGTH 512u

// recursive mkdir /////////////////////////////////////////////////////////////
int mkdir(const std::string dirName)
{
  if (!dirName.empty() and access(dirName.c_str(), R_OK|W_OK|X_OK))
  {
    mode_t mode755;
    char   tmp[MAX_PATH_LENGTH];
    char   *p = NULL;
    size_t len;

    mode755 = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;

    snprintf(tmp, sizeof(tmp), "%s", dirName.c_str());
    len = strlen(tmp);
    if(tmp[len - 1] == '/')
    {
      tmp[len - 1] = 0;
    }
    for(p = tmp + 1; *p; p++)
    {
      if(*p == '/')
      {
        *p = 0;
        ::mkdir(tmp, mode755);
        *p = '/';
      }
    }

    return ::mkdir(tmp, mode755);
  }
  else
  {
    return 0;
  }
}

std::string dirname(const std::string &s)
{
  constexpr char sep = '/';
  size_t         i   = s.rfind(sep, s.length());

  if (i != std::string::npos)
  {
    return s.substr(0, i);
  }
  else
  {
    return "";
  }
}

void makeFileDir(const std::string filename, GridBase *g)
{
  bool doIt = true;

  if (g)
  {
    doIt = g->IsBoss();
  }
  if (doIt)
  {
    std::string dir    = dirname(filename);
    int         status = mkdir(dir);

    if (status)
    {
      std::cout <<  "cannot create directory '" + dir << std::endl;
    }
  }
}


class A2AVectorsIo
{
  public:
    struct Record: Serializable
  {
    GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
        unsigned int, index);
    Record(void): index(0) {}
  };
  public:
    template <typename Field>
      static void write(const std::string fileStem, std::vector<Field> &vec, 
          const bool multiFile, const int trajectory = -1);
    template <typename Field>
      static void read(std::vector<Field> &vec, const std::string fileStem,
          const bool multiFile, const int trajectory = -1);
  private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
        const bool multiFile)
    {
      std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

      if (multiFile)
      {
        return stem + t;
      }
      else
      {
        return stem + t + ".bin";
      }
    }
};


  template <typename Field>
void A2AVectorsIo::write(const std::string fileStem, std::vector<Field> &vec, 
    const bool multiFile, const int trajectory)
{
  Record       record;
  GridBase     *grid = vec[0].Grid();
  ScidacWriter binWriter(grid->IsBoss());
  std::string  filename = vecFilename(fileStem, trajectory, multiFile);

  if (multiFile)
  {
    std::string fullFilename;

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

      std::cout << GridLogMessage << "Writing vector " << i << std::endl;
      makeFileDir(fullFilename, grid);
      binWriter.open(fullFilename);
      record.index = i;
      binWriter.writeScidacFieldRecord(vec[i], record);
      binWriter.close();
    }
  }
  else
  {
    makeFileDir(filename, grid);
    binWriter.open(filename);
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      std::cout << GridLogMessage << "Writing vector " << i << std::endl;
      record.index = i;
      binWriter.writeScidacFieldRecord(vec[i], record);
    }
    binWriter.close();
  }
}

  template <typename Field>
void A2AVectorsIo::read(std::vector<Field> &vec, const std::string fileStem, 
    const bool multiFile, const int trajectory)
{
  Record       record;
  ScidacReader binReader;
  std::string  filename = vecFilename(fileStem, trajectory, multiFile);

  if (multiFile)
  {
    std::string fullFilename;

    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

      // LOG(Message) << "Reading vector " << i << std::endl;
      std::cout << "Reading vector " << i << std::endl;
      binReader.open(fullFilename);
      binReader.readScidacFieldRecord(vec[i], record);
      binReader.close();
      if (record.index != i)
      {
        std::cout << "vector index mismatch" << std::endl;
        // HADRONS_ERROR(Io, );
      }
    }
  }
  else
  {
    binReader.open(filename);
    for (unsigned int i = 0; i < vec.size(); ++i)
    {
      // LOG(Message) << "Reading vector " << i << std::endl;
      std::cout << "Reading vector " << i << std::endl;
      binReader.readScidacFieldRecord(vec[i], record);
      if (record.index != i)
      {
        std::cout << "vector index mismatch" << std::endl;
        // HADRONS_ERROR(Io, "vector index mismatch");
      }
    }
    binReader.close();
  }
}





template<class T>
void writeScidac(T& field, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const
  std::cout << "writing to" << filename << std::endl;
  makeFileDir(filename, field.Grid());
  // if(field.Grid()->IsBoss()) {
  //   std::string base_dir = filename.substr(0, filename.rfind('/'));
  //   std::cout << "base_dir: " << base_dir << std::endl;
  //   assert(dirExists(base_dir));
  //   // system(("mkdir -p " + base_dir).c_str());      //  sometimes has error
  // }
  // MPI_Barrier(MPI_COMM_WORLD);

  emptyUserRecord record;
  ScidacWriter WR(field.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
  WR.writeScidacFieldRecord(field, record);
  WR.close();
};

template<class T>
void readScidac(T& field, const std::string &filename){
  std::cout << "reading from" << filename << std::endl;
  assert(dirExists(filename));
  emptyUserRecord record;
  ScidacReader RD;
  RD.open(filename);
  RD.readScidacFieldRecord(field, record);
  RD.close();
};







// For Five dimensional rbGrid, convert double -> float
GridRedBlackCartesian *FrbGrid_d2f(GridRedBlackCartesian *FrbGrid_d) {
  std::vector<int> tmp = FrbGrid_d->FullDimensions().toVector();
  int Ls = tmp[0];
  tmp.erase(tmp.begin());
  Coordinate fcoor(tmp);

  tmp = FrbGrid_d->_processors.toVector();   // 5d
  tmp.erase(tmp.begin());
  Coordinate processors(tmp);

  GridCartesian *UGrid_f = SpaceTimeGrid::makeFourDimGrid(fcoor, GridDefaultSimd(Nd,vComplexF::Nsimd()), processors); 
  GridRedBlackCartesian *FrbGrid_f = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid_f);
  return FrbGrid_f;
}


void writeScidac_evec_d2f(LatticeFermionD& evec, const std::string &filename){ // because of writeScidacFieldRecord, field cannot be const

  static GridRedBlackCartesian *FrbGrid_f = FrbGrid_d2f((GridRedBlackCartesian *)evec.Grid());
  LatticeFermionF evec_f(FrbGrid_f);
  precisionChange(evec_f, evec);


  emptyUserRecord record;
  ScidacWriter WR(evec_f.Grid()->IsBoss()); // the parameter is necessary for writer(but not for reader) when using multiple nodes
  WR.open(filename);
  WR.writeScidacFieldRecord(evec_f, record);
  WR.close();
};





void readScidac_evec_f2d(LatticeFermionD& evec, const std::string &filename){

  static GridRedBlackCartesian *FrbGrid_f = FrbGrid_d2f((GridRedBlackCartesian *)evec.Grid());

  LatticeFermionF evec_f(FrbGrid_f);

  assert(dirExists(filename));
  emptyUserRecord record;
  ScidacReader RD;
  RD.open(filename);
  RD.readScidacFieldRecord(evec_f, record);
  RD.close();

  precisionChange(evec, evec_f);
};

void write_evecs_evals(const std::string &base_dir, const std::vector<double> &evals, std::vector<LatticeFermionD> &evecs) {
  system(("mkdir -p " + base_dir).c_str());

    // write evals
  if(evecs[0].Grid()->IsBoss()) {
    std::ofstream f(base_dir + "/evals.txt");
    for (const double &e: evals) f << e << "\n";
    f.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // write evecs
  for(int i=0; i<evecs.size(); ++i) {
    writeScidac(evecs[i], base_dir + "/evec_" + std::to_string(i));
  }
}




void write_evecs_evals_d2f(const std::string &base_dir, const std::vector<double> &evals, std::vector<LatticeFermionD> &evecs) {
  system(("mkdir -p " + base_dir).c_str());

    // write evals
  if(evecs[0].Grid()->IsBoss()) {
    std::ofstream f(base_dir + "/evals.txt");
    for (const double &e: evals) f << e << "\n";
    f.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // write evecs
  for(int i=0; i<evecs.size(); ++i) {
    writeScidac_evec_d2f(evecs[i], base_dir + "/evec_" + std::to_string(i));
  }
}


// read float -> float, or double -> double
template <class FermionType>
void read_evecs_evals(const std::string &base_dir, GridRedBlackCartesian *FrbGrid_fd, std::vector<double> &evals, std::vector<FermionType> &evecs) {
    // read evals
    evals.clear();
    std::ifstream f(base_dir + "/evals.txt");
    double tmp;
    while(f>>tmp) {
      evals.push_back(tmp);
    }

    int N_eval = evals.size();
    evecs.resize(N_eval, FrbGrid_fd);
    // read evecs
    for(int i=0; i<N_eval; ++i) {
      readScidac(evecs[i], base_dir + "/evec_" + std::to_string(i));
      evecs[i].Checkerboard() = Odd;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    f.close();

    std::cout << "norm2: evecs[0]: " << norm2(evecs[0]) << " evecs[-1]: " << norm2(evecs.back())  << std::endl;
}


















}
