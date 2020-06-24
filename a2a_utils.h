#pragma once

#include <Grid/Grid.h>

void print_memory() {
  using namespace std;

  map<string, double> memInfo; // item -> memory in GB
  ifstream f("/proc/meminfo");
  while(f) {
    string item, tmp;
    long long int size;

    f >> item >> size;
    getline(f, tmp); // move to the next line

    if(!item.empty()) {
      item.pop_back(); // remove last character
      memInfo[item] = double(size) / 1024 / 1024;
    }
  }
  // for(auto [s, i]: memInfo) std::cout << s << " " << i << " KB" << std::endl;

  std::cout << "MemTotal: " << memInfo["MemTotal"] << " GB " << std::endl;
  std::cout << "Non cache/buffer MemUsed: " << memInfo["MemTotal"] - memInfo["MemFree"] - memInfo["Buffers"] - memInfo["Cached"] << " GB " << std::endl;
  std::cout << "Buffers/Cached: " << memInfo["Buffers"] + memInfo["Cached"] << " GB " << std::endl;
  std::cout << "MemAvailable: " << memInfo["MemAvailable"] << " GB " << std::endl; // Memavailable is roughly memfree + buffers + cached,

}


namespace Grid {
// A flaw in grid, the field argument of write cannot be const.
void writeGaugeField(LatticeGaugeField &U, const std::string &fileName){
  NerscIO::writeConfiguration(U,fileName,0,0);
}

void readGaugeField(LatticeGaugeField &U, const std::string &fileName){
  FieldMetaData header;
  NerscIO::readConfiguration(U,header,fileName);
}


void localIndexToLocalGlobalCoor(GridBase *grid, int ss, Coordinate &lcoor, Coordinate &gcoor) {
  // ss is local index; parallel_for(int ss=0; ss<ret.Grid()->lSites(); ss++)
  lcoor.resize(4);
  gcoor.resize(4);
  grid->LocalIndexToLocalCoor(ss, lcoor);
  Coordinate processor_coor;
  grid->ProcessorCoorFromRank(grid->ThisRank(), processor_coor);
  grid->ProcessorCoorLocalCoorToGlobalCoor(processor_coor, lcoor, gcoor);
}

}
