CRUJRA v1.1 - v2.2 -- Julia E.M.S. Nabel, July 2018-2021 (julia.nabel@mpimet.mpg.de)
CRUJRA v2.3 -- Stefanie Falk, July 2022 (stefanie.falk@lmu.de)

*** crujra (raw) data is provided via the trendy ftp server (trendy-vXX@trendy.ex.ac.uk)
sftp trendy-vXX@trendy.ex.ac.uk
cd /input/CRUJRA20YY 
get *var* (var in dlwrf pre spfh tmax tmin vgrd ugrd)
* no more dswrf => now tswrf and fd (TR wrote in an email on 6/22/21 that it is not sensible to use fd in jsbach 3.2)
cd /input/CRUJRA20YY/Radiation-fields
get *tswrf*

*** gunzip with 
vars=( tmin tmax dlwrf pre spfh vgrd ugrd )
for var in ${vars[@]}; do 
  sbatch gunzip_climate_data_CRUJRA.bash ${var}
done
* rename tswrf file to be able to handle it the same way as the other vars
for year in {1901..20ZZ}; do
  mv tswrf_vXX_${year}.nc crujra.vx.y.5d.tswrf.${year}.365d.noc.nc
done

*** remap with 
remap_CRUJRA.bash 
* possibly separately remap tswrf

*** prepare for jsbach with 
./calculate_jsbach_forcing_CRUJRA.bash

*** created a global CO2 jsbach input file from the GCB20YY global_co2_ann_1700_20ZZ.txt file (download from trendy-vXX@trendy.ex.ac.uk)
create_global_CO2_file.bash

Data for v2.3 can be found here:
/pool/data/JSBACH/jsbalone_forcing/T63/CRUJRA
older versions archived:
packems -o CRUJRA_v2.0_T63 CRUJRA/data/crujra_v2.0 \
        CRUJRA/doc/README_v2.0.txt CRUJRA/scripts/scripts_v2.0 \
        -S /arch/mj0060/jsbalone_forcing -I INDEX_CRUJRA_v2.0_T63.txt

Japanese Reanalysis:
Kobayashi, S., Ota, Y., Harada, Y., Ebita, A., Moriya, M., Onoda,
H., Onogi, K., Kamahori, H., Kobayashi, C., Endo, H., Miyaoka,
K., and Takahashi, K.: The JRA-55 Reanalysis: General Specifications
and Basic Characteristics, J. Meteorol. Soc. Jpn., 93,
5–48, https://doi.org/10.2151/jmsj.2015-001, 2015.

CRU TS v4.03 dataset:
Harris, I., Osborn, T. J., Jones, P., and Lister, D.: Version 4 of the
CRU TS monthly high-resolution gridded multivariate climate
dataset, Sci. Data, 7, 109,  https://doi.org/10.1038/s41597-020-
0453-3, 2021.