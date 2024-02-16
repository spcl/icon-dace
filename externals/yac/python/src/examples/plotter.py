#!/usr/bin/env python3

##
# @file plotter.py
# @brief defining a model component that plots fields from other components
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


import os
from yac import *
import numpy as np
import matplotlib.pyplot as plt

class Plotter:
    def __init__(self, variables, outdir = ".",
                 bounds = [0, 2*np.pi, -0.5*np.pi, 0.5*np.pi],
                 resolution = (360,180),
                 yac = None):
        self.yac = yac or YAC.default_instance
        self.comp = self.yac.predef_comp("plotter")
        self.variables = variables
        self.bounds = bounds
        self.resolution = resolution
        self.outdir = outdir

    def setup(self):
        self.x = np.linspace(self.bounds[0],self.bounds[1],self.resolution[0])
        self.y = np.linspace(self.bounds[2],self.bounds[3],self.resolution[1])
        grid = Reg2dGrid(f"plot_grid", self.x, self.y)
        self.points = grid.def_points(Location.CORNER, self.x, self.y)

    def def_couples(self):
        nnn = InterpolationStack()
        nnn.add_nnn(NNNReductionType.AVG, 1, 1.)

        self.fields = []
        for var in self.variables:
            timestep = self.yac.get_field_timestep(*var)
            collection_size = self.yac.get_field_collection_size(*var)
            self.fields.append(Field.create(var[2], self.comp, self.points, 1,
                                            timestep, TimeUnit.ISO_FORMAT))
            self.yac.def_couple(*var,
                                "plotter", "plot_grid", var[2],
                                timestep, TimeUnit.ISO_FORMAT, 0, nnn)

        for field in self.fields:
            for i in range(field.collection_size):
                os.makedirs(self.outdir +"/"+ f"{field.name}_{i}", exist_ok=True)

    def step(self):
        for field in self.fields:
            time = field.datetime
            print(f"plotting {field.name} at {time}")
            buf, info = field.get()
            for i in range(buf.shape[0]):
                plt.imshow(buf[i,:].reshape((len(self.y),len(self.x)))[::-1,:],
                           extent=self.bounds)
                plt.title(f"{field.name} - {time}")
                plt.colorbar()
                plt.savefig(self.outdir +"/"+ f"{field.name}_{i}/{time}.png")
                plt.clf()
        return field.datetime
