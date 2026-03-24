set -xe
mkdir -p $SCRATCH/iconbk
[ -d experiments ] && mv experiments $SCRATCH/iconbk/experiments_$(date +%Y%m%d_%H%M%S)
rm -rf bin experiments lib mod externals
