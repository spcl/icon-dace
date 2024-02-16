#!/usr/bin/env python3

##
# @file driver.py
# @brief a driver code combining different compoents
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

import sys
import datetime
import numpy as np
from threading import Thread

from yac import *

## Initialization
def_calendar(Calendar.PROLEPTIC_GREGORIAN)

class Driver:
    def __init__(self, start=None, end=None, verbose=0, yac = None, multithreading = False):
        self.yac = yac or YAC.default_instance
        self.yac.def_datetime(start, end)
        self.verbose = verbose
        self.multithreading = multithreading

    def run(self, *components):
        # finish component definition by adding a own component
        self.yac.def_comps()
        if self.verbose > 0: print("Components: ", self.yac.component_names)

        # let the components setup their grids, fields etc.
        for comp in components:
            comp.setup()

        self.yac.sync_def()

        # some components might add something after synchronization
        for comp in components:
            comp.def_couples()

        self.yac.enddef()

        end_datetime = datetime.datetime.fromisoformat(self.yac.end_datetime)
        for c in components: c.datetime = datetime.datetime.fromisoformat(self.yac.start_datetime)

        if self.multithreading:
            def run_comp(comp):
                while comp.datetime <= end_datetime:
                    comp.datetime = datetime.datetime.fromisoformat(comp.step())
            threads = [Thread(target=run_comp, args=(comp,)) for comp in components]
            for t in threads: t.start()
            for t in threads: t.join()
        else:
            while True:
                min_comp = min(components, key=lambda c: c.datetime)
                if min_comp.datetime > end_datetime:
                    break
                nxt = min_comp.step()
                min_comp.datetime = datetime.datetime.fromisoformat(nxt)
