#!/bin/ksh
specs='CH3OH_VMR_inst,CH3OH,SO2_VMR_inst,SO2,N2O_VMR_inst,N2O,CO_VMR_inst,CO,H2O2_VMR_inst,H2O2,CH4_VMR_inst,CH4,CH3OOH_VMR_inst,CH3OOH,N2O5_VMR_inst,N2O5,HNO3_VMR_inst,HNO3,CH3O2_VMR_inst,CH3O2,NO3_VMR_inst,NO3,O3_VMR_inst,O3,OH_VMR_inst,OH,HO2_VMR_inst,HO2,NO_VMR_inst,NO,NO2_VMR_inst,NO2'
cdo chname,$specs mozart4geos5_r2b4.nc mozart4geos5_r2b4_iconnames.nc
