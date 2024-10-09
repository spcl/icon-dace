# create_hdmodel_tarfile.com - Generates tar file of the HD model repository for publication in Zenodo 
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
# This script generates a HD model tar file based on the git commit $1
# e.g. f888528
# Whether a tag can be used instead, should be checked, e.g HD5.0

set VS=5_2_2
# 
cd ..
git archive --prefix=hd_vs5/ -o hdmodel_${VS}.tgz --worktree-attributes $1
mv hdmodel_${VS}.tgz /work/gg0302/g260122/HD/input
echo "hdmodel_${VS}.tgz was generated in /work/gg0302/g260122/HD/input"
#

