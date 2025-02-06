#!/bin/bash

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

tag=$1

echo $tag

# we could also checkout from remote
git archive --format=zip -o yac-$tag.zip --prefix=yac-$tag/ $tag

mv yac-$tag.zip ..
cd ..

# unpack archive and clean it up

unzip yac-$tag.zip

cd yac-$tag

# a) clean up main directory

rm -f yac-distribution.sh

# b) clean up example directory

rm -f *.sh

# create an archive file for distributing

rm yac-$tag.zip

tar -czf yac-$tag.tgz yac-$tag

