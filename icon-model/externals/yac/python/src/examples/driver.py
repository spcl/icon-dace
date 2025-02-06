#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

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
