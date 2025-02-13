#include <fstream>
#include <iostream>
#include <string>


#include "cloud_optics_fn_438.h"
#include "serdeco.h"

constexpr char ROOT[] = "/home/alex/backup-icon-model/experiments/exclaim_ape_R2B09/cloud_optics/";

int main() {
  printf("Here\n");
  fflush(stdout);
  config_type config;
  {
    std::ifstream data(std::string(ROOT) + "config_t0.1.data");
    serde::deserialize(&config, data);
  }
  thermodynamics_type thermodynamics;
  {
    std::ifstream data(std::string(ROOT) + "thermodynamics_t0.1.data");
    serde::deserialize(&thermodynamics, data);
  }
  cloud_type cloud;
  {
    std::ifstream data(std::string(ROOT) + "cloud_t0.1.data");
    serde::deserialize(&cloud, data);
  }
  
  double * od_lw_cloud=(double *) malloc(46080*sizeof(double));
  
  double* od_sw_cloud=(double *) malloc(40320*sizeof(double));
  
  double* ssa_sw_cloud=(double *) malloc(40320*sizeof(double));

  double ssa_lw_cloud[1];
  double g_lw_cloud[1];
  double *g_sw_cloud=(double *) malloc(40320*sizeof(double));

  auto* h = __dace_init_cloud_optics_fn_438(&cloud,&config,g_lw_cloud,g_sw_cloud,od_lw_cloud,od_sw_cloud,ssa_lw_cloud,ssa_sw_cloud,&thermodynamics,32,1,90,90);
  
  __program_cloud_optics_fn_438(h,&cloud,&config,g_lw_cloud,g_sw_cloud,od_lw_cloud,od_sw_cloud,ssa_lw_cloud,ssa_sw_cloud,&thermodynamics,32,1,90,90);
  
  __dace_exit_cloud_optics_fn_438(h);

  {
    std::ofstream data("od_lw_cloud.got");
    data << "# entries" << std::endl;
    for (int i=0; i<46080; ++i) data << od_lw_cloud[i] << std::endl;
  }

  {
    std::ifstream data(std::string(ROOT) + "od_lw_cloud_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 46080; ++i) {
      std::getline(data, line);
      od_lw_cloud[i] = std::stod(line);
    }
  }
  

  {
    std::ofstream data("od_lw_cloud.want");
    data << "# entries" << std::endl;
    for (int i=0; i<46080; ++i) data << od_lw_cloud[i] << std::endl;
  }

  
  {
    std::ofstream data("od_sw_cloud.got");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << od_sw_cloud[i] << std::endl;
  }



  {
    std::ifstream data(std::string(ROOT) + "od_sw_cloud_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 40320; ++i) {
      std::getline(data, line);
      od_sw_cloud[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("od_sw_cloud.want");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << od_sw_cloud[i] << std::endl;
  }

  {
    std::ofstream data("ssa_sw_cloud.got");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << ssa_sw_cloud[i] << std::endl;

  }

  {
    std::ifstream data(std::string(ROOT) + "ssa_sw_cloud_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 40320; ++i) {
      std::getline(data, line);
      ssa_sw_cloud[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("ssa_sw_cloud.want");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << ssa_sw_cloud[i] << std::endl;
  }


  {
    std::ofstream data("g_sw_cloud.got");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << g_sw_cloud[i] << std::endl;

  }

  
  {
    std::ifstream data(std::string(ROOT) + "g_sw_cloud_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 40320; ++i) {
      std::getline(data, line);
      g_sw_cloud[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("g_sw_cloud.want");
    data << "# entries" << std::endl;
    for (int i=0; i<40320; ++i) data << g_sw_cloud[i] << std::endl;
  }
  free(od_lw_cloud);
  free(od_sw_cloud);
  free(ssa_sw_cloud);
  free(g_sw_cloud);

  return EXIT_SUCCESS;
}
