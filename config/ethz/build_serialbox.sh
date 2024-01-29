git clone git@github.com:GridTools/serialbox.git
cd serialbox
cmake \
    -DBOOST_ROOT:PATH=${BOOST_ROOT} \
    -DBoost_NO_BOOST_CMAKE=OFF \
    -DSERIALBOX_ENABLE_C=ON \
    -DSERIALBOX_ENABLE_FORTRAN=ON \
    -DSERIALBOX_BUILD_SHARED=ON \
    -DSERIALBOX_ENABLE_EXPERIMENTAL_FILESYSTEM=ON \
    -DCMAKE_INSTALL_PREFIX:PATH=$(pwd) \
    .

make -j 8 && make install

cd -
