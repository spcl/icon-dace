#include <fstream>
#include <iostream>
#include <string>

#include "crop_cloud_fraction.h"
#include "serde.h"

constexpr char ROOT[] = "/home/alex/backup-icon-model/experiments/exclaim_ape_R2B09/";

int main() {
  cloud_type cloud;
  {
    std::ifstream data(std::string(ROOT) + "cloud_t0.1.data");
    serde::deserialize(&cloud, data);
  }

 
  auto* h = __dace_init_crop_cloud_fraction(&cloud,9.9999999999999995E-007,     1.0000000000000001E-009, 32, 1, 32,1);
  __program_crop_cloud_fraction(h,&cloud,9.9999999999999995E-007,1.0000000000000001E-009, 32, 1, 32,1);
  __dace_exit_crop_cloud_fraction(h);

  {
    std::ofstream data("cloud.got");
    data << serde::serialize(&cloud)<<std::endl;
  }
  
  cloud_type cloud2;
  {
    std::ifstream data(std::string(ROOT) + "cloud_t1.1.data");
    serde::deserialize(&cloud2, data);
  }

  {
    std::ofstream data("cloud.want");
    data << serde::serialize(&cloud2)<<std::endl;
  }
  return EXIT_SUCCESS;
}
