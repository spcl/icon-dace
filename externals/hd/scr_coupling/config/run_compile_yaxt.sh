autoreconf -i
 
/configure --prefix=/work/gg0302/g260062/GCOAST_oas5/yaxt --without-regard-for-quality --disable-shared
make -j 4
make install
