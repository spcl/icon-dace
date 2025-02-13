#include <fstream>
#include <iostream>
#include <string>


#include "add_aerosol_optics.h"
#include "serdeae.h"

constexpr char ROOT[] = "/home/alex/backup-icon-model/experiments/exclaim_ape_R2B09/add_aerosol_optics/";

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
  gas_type gas;
  {
    std::ifstream data(std::string(ROOT) + "gas_t0.1.data");
    serde::deserialize(&gas, data);
  }
  aerosol_type aerosol;
  {
    std::ifstream data(std::string(ROOT) + "aerosol_t0.1.data");
    serde::deserialize(&aerosol, data);
  }
  
  printf("Here2\n");
  fflush(stdout);
  double * od_lw=(double *) malloc(403200*sizeof(double));
  
  double* od_sw=(double *) malloc(322560*sizeof(double));
  
  double* ssa_sw=(double *) malloc(322560*sizeof(double));

  double ssa_lw[1];
  double g_lw[1];
  double *g_sw=(double *) malloc(322560*sizeof(double));

  {
    std::ifstream data(std::string(ROOT) + "od_lw_t0.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 403200; ++i) {
      std::getline(data, line);
      od_lw[i] = std::stod(line);
    }
  }
  {
    std::ifstream data(std::string(ROOT) + "od_sw_t0.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 322560; ++i) {
      std::getline(data, line);
      od_sw[i] = std::stod(line);
    }
  }
  {
    std::ifstream data(std::string(ROOT) + "ssa_sw_t0.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 322560; ++i) {
      std::getline(data, line);
      ssa_sw[i] = std::stod(line);
    }
  }

  printf("Here3\n");
  fflush(stdout);
  auto* h = __dace_init_add_aerosol_optics(&aerosol,&config,g_lw,g_sw,&gas,od_lw,od_sw,ssa_lw,ssa_sw,&thermodynamics,32,1,90,90);
  
  printf("Here\4");
  fflush(stdout);
  __program_add_aerosol_optics(h,&aerosol,&config,g_lw,g_sw,&gas,od_lw,od_sw,ssa_lw,ssa_sw,&thermodynamics,32,1,90,90);
  
  printf("Here5\n");
  fflush(stdout);
  __dace_exit_add_aerosol_optics(h);

  {
    std::ofstream data("od_lw.got");
    data << "# entries" << std::endl;
    for (int i=0; i<403200; ++i) data << od_lw[i] << std::endl;
  }

  {
    std::ifstream data(std::string(ROOT) + "od_lw_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 403200; ++i) {
      std::getline(data, line);
      od_lw[i] = std::stod(line);
    }
  }
  

  {
    std::ofstream data("od_lw.want");
    data << "# entries" << std::endl;
    for (int i=0; i<403200; ++i) data << od_lw[i] << std::endl;
  }

  
  {
    std::ofstream data("od_sw.got");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << od_sw[i] << std::endl;
  }



  {
    std::ifstream data(std::string(ROOT) + "od_sw_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 322560; ++i) {
      std::getline(data, line);
      od_sw[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("od_sw.want");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << od_sw[i] << std::endl;
  }

  {
    std::ofstream data("ssa_sw.got");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << ssa_sw[i] << std::endl;

  }

  {
    std::ifstream data(std::string(ROOT) + "ssa_sw_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 322560; ++i) {
      std::getline(data, line);
      ssa_sw[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("ssa_sw.want");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << ssa_sw[i] << std::endl;
  }


  {
    std::ofstream data("g_sw.got");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << g_sw[i] << std::endl;

  }

  
  {
    std::ifstream data(std::string(ROOT) + "g_sw_t1.1.data");
    std::string line;
    std::getline(data, line);
    for (int i = 0; i < 322560; ++i) {
      std::getline(data, line);
      g_sw[i] = std::stod(line);
    }
  }

  {
    std::ofstream data("g_sw.want");
    data << "# entries" << std::endl;
    for (int i=0; i<322560; ++i) data << g_sw[i] << std::endl;
  }
  free(od_lw);
  free(od_sw);
  free(ssa_sw);
  free(g_sw);

  return EXIT_SUCCESS;
}
