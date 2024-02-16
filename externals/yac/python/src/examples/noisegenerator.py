#!/usr/bin/env python3

##
# @file noisegenerator.py
# @brief a model component adding a field `noise`, filled with random values
#
# @copyright Copyright  (C)  2022 DKRZ, MPI-M
#
# @author Nils-Arne Dreier <dreier@dkrz.de>
#

#
# Keywords:
# Maintainer: Nils-Arne Dreier <dreier@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
#
# This file is part of YAC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from yac import *
import numpy as np

class NoiseGenerator:
    def __init__(self, timestep, yac = None):
        self.yac = yac or YAC.default_instance
        self.comp = self.yac.predef_comp("noisegenerator")
        self.timestep = timestep

    def setup(self):
        global grid, noise_field
        x = np.linspace(0,2*np.pi,360)
        y = np.linspace(-0.5*np.pi,0.5*np.pi,180)
        grid = Reg2dGrid(f"noise_grid", x, y)
        points = grid.def_points(Location.CORNER, x, y)

        noise_field = Field.create("noise", self.comp, points, 1,
                                   self.timestep, TimeUnit.ISO_FORMAT)

    def def_couples(self):
        pass

    def step(self):
        print(f"NoiseGenerator: Lets make some noise!!!")
        noise_field.put(np.random.rand(grid.nbr_corners, 1))
        return noise_field.datetime
