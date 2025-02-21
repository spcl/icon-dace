INCS="-I./src/include -I./externals/ecrad/include -I./externals/ecrad/radiation -I./externals/yaxt/include -I./externals/comin/include -I./externals/ppm/include/f77"
DEFS="-DNOMPI -D__NO_ICON_COMIN__"

# For ecrad.
for f in $(find externals/ecrad | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
for f in $(find externals/rte-rrtmgp | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done

# For velocity tendencies.
for f in $(find src | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
for f in $(find externals/comin/src | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
for f in $(find externals/fortran-support/src | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
for f in $(find externals/mtime/src | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
for f in $(find externals/ppm/src | grep -E "\.F90\$|\.f90\$"); do
  gfortran $INCS $DEFS -cpp -E -P "${f}" > "${f}.tmp" && mv "${f}.tmp" "${f}" && sed -i '' '/./,$!d' "${f}"
done
