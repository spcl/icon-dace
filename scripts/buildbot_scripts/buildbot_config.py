# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

from util import relate_with_config, relate_without_config
from builder import Builder
from experiment import Experiment

import pandas as pd
import numpy as np
import pickle
import sys
import os

class BuildbotConfig(object):
    def __init__(self):
        self.data = pd.DataFrame()
        self.builder_meta = {}
        self.machine_meta = {}

    @classmethod
    def from_pickle(cls, filename):
        picklefile = open(str(filename), "rb")
        return pickle.load(picklefile)

    def to_pickle(self, filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        picklefile = open(str(filename), "wb")
        pickle.dump(self, picklefile)
    
    def builder_by_config(self, config, relate):
        return [b for b, v in self.builder_meta.items() if relate(config, v.config)]

    def get_indices_from_args(self, machines=None, builders=None, with_config=None, without_config=None):
        midx = machines if machines else slice(None)

        # note: utility dict comparisons return True if lhs is None (treat as empty set)
        config_builders = self.builder_by_config(with_config, relate_with_config)
        without_config_builders = self.builder_by_config(without_config, relate_without_config)

        # add experiment to all builders according to config flags (with and without possibly)
        bidx = set(config_builders).intersection(set(without_config_builders))

        # only add experiments to builders with flag=Active (at the time when experiment is added)
        active_builders = set([b for b, v in self.builder_meta.items() if v.flag and v.flag.lower() == "active"])
        bidx = bidx.intersection(active_builders)
        
        # if a list of builders is given, limit to that list
        if builders:
            bidx = set(builders).intersection(bidx)

        # convert bidx to list for further use
        bidx = list(bidx)

        return midx, bidx

    def add_experiments(self, exps, builders=None, machines=None, with_config=None, without_config=None, run_flags=None):
        for exp in exps:
            self.add_experiment(exp, builders=builders, machines=machines, with_config=with_config, without_config=without_config, run_flags=run_flags)

    def add_experiment(self, exp, builders=None, machines=None, with_config=None, without_config=None, run_flags=None):
        midx, bidx = self.get_indices_from_args(machines, builders, with_config, without_config)


        # if experiment exists, add object to selected builders/machines
        if exp not in self.data.index:
            df = pd.DataFrame(None, columns=self.data.columns, index=[exp])
            self.data = pd.concat([self.data, df], axis=0)
            self.data.sort_index(axis=0, inplace=True)

        try:
            # first get all the builder names that are in the intersection of
            # machines and builders determined by with/without configure flags.
            # Then add a unique experiment to the table. The experiments take
            # the builder object in their constructor.
            builder_names = self.data.loc[(exp, (midx, bidx))].index.get_level_values(1)
            self.data.loc[(exp, (midx, bidx))] = [Experiment(exp, self.builder_meta[bname], run_flags) for bname in builder_names]

            print("added experiment {} to {} configuration".format(exp, ", ".join(bidx)))
        except:
            print("Could not find machine or builder - skipping experiments")


    def remove_experiments(self, exps, builders=None, machines=None, with_config=None, without_config=None):
        for exp in exps:
            self.remove_experiment(exp, builders=builders, machines=machines, with_config=with_config, without_config=without_config)

    def remove_experiment(self, exp, builders=None, machines=None, with_config=None, without_config=None):
        midx, bidx = self.get_indices_from_args(machines, builders, with_config, without_config)

        # if experiment exists, remove object from selected builders/machines
        if exp not in self.data.index:
            print("requested to remove experiment {}, but does not exist.".format(exp))
        else:
            self.data.loc[exp, (midx, bidx)] = np.nan

        # clean up the data by removing all rows that only contain nan
        self.data.dropna(how="all", inplace=True)


    def add_builders(self, builders, machine, script, config=None, flag=None):
        if not isinstance(builders, list):
            builders = [builders]

        if machine not in self.machine_meta.keys():
            print("Could not find machine {} when adding builders {}. Check for typos or add machine {} first.".format(machine, " ".join(builders), machine))
            sys.exit(1)

        mcol = pd.MultiIndex.from_product([[machine], builders], names=("machine", "builder"))
        df = pd.DataFrame(None, columns=mcol)
        self.data = pd.concat([self.data, df], axis=1, sort=False)
        self.data.sort_index(axis=1, inplace=True)

        for builder in builders:
            self.builder_meta[builder] = Builder(builder, machine, script, config, flag)

        print("added builder(s) {} to {} configuration".format(", ".join(builders), machine))


    def add_machine(self, machine, queue=None):
        mdict = {
            "queue": queue
        }
        self.machine_meta[machine] = mdict


    def set_builder_flag(self, builders, flag):
        for b, bobj in self.builder_meta.items():
            if b in builders:
                bobj.flag = flag

    def build_builder(self, builder):
        if builder not in self.builder_meta.keys():
            print("builder is not part of this buildbot configuration")
            sys.exit(1)

        bobj = self.builder_meta[builder]
        status = bobj.configure_make()

        if status != 0:
            print("failed to build builder {}".format(builder))
            sys.exit(1)


    def make_builder_runscripts(self, builder):
        exps = self.get_experiments_by_builder(builder)
        exp_list = []
        for exp in exps.flatten():
            status = exp.make_runscript()
            if status == 0: exp_list.append(exp.run_name)
            else:
                print("failed to create experiment {}".format(exp.name))
                sys.exit(1)

        return " ".join(exp_list)

    def get_experiments_by_builder(self, builder):
        if builder not in self.builder_meta.keys():
            print("builder {} is not part of this buildbot configuration".format(builder))
            sys.exit(1)

        bobj = self.builder_meta[builder]
        machine = bobj.machine

        exps = self.data.loc[:, (machine, builder)].dropna(how="all").values

        return exps

    def list_experiments_by_builder(self, builder):
        exps = self.get_experiments_by_builder(builder)

        return [exp.name for exp in exps.flatten()]

    def add_dependency_manager(self, source_experiment, target_experiment, source_builder=None, target_builder=None, builders=None, machines=None, with_config=None, without_config=None):
        # explicit mode for cross-builder support (not yet implemented)
        if(source_builder and target_builder):
            if(len(source_experiment) == 1 and len(target_experiment) == 1):
                self.add_dependency(source_experiment[0], target_experiment[0], source_builder, target_builder)
            elif(len(source_experiment) > 1 and len(target_experiment) == 1):
                self.add_dependency_many2one(source_experiment, target_experiment[0], source_builder, target_builder)
            elif(len(source_experiment) == 1 and len(target_experiment) > 1):
                self.add_dependency_one2many(source_experiment[0], target_experiment, source_builder, target_builder)
            else:
                print("Cannot process dependency from {} to {}. Only 1:1, 1:N and N:1 relations are allowed.".format(source_experiment, target_experiment))
                sys.exit(1)

        # if both source and target builders are undefined, assume same builder dependencies
        elif(not source_builder and not target_builder):
            midx, bidx = self.get_indices_from_args(machines, builders, with_config, without_config)
            builder_names = self.data.loc[(slice(None), (midx, bidx))].columns.get_level_values(1)
            for bn in builder_names:
                if(len(source_experiment) == 1 and len(target_experiment) == 1):
                    self.add_dependency(source_experiment[0], target_experiment[0], bn, bn)
                elif(len(source_experiment) > 1 and len(target_experiment) == 1):
                    self.add_dependency_many2one(source_experiment, target_experiment[0], bn, bn)
                elif(len(source_experiment) == 1 and len(target_experiment) > 1):
                    self.add_dependency_one2many(source_experiment[0], target_experiment, bn, bn)
                else:
                    print("Cannot process dependency from {} to {}. Only 1:1, 1:N and N:1 relations are allowed.".format(source_experiment, target_experiment))
                    sys.exit(1)

        else:
            print("source_builder and target_builder have to be defined or undefined to switch between explicit and assumed dependency mode")
            sys.exit(1)


    def add_dependency_many2one(self, source_experiments, target_experiment, source_builder, target_builder):
        for se in source_experiments:
            self.add_dependency(se, target_experiment, source_builder, target_builder)


    def add_dependency_one2many(self, source_experiment, target_experiments, source_builder, target_builder):
        for te in target_experiments:
            self.add_dependency(source_experiment, te, source_builder, target_builder)


    def add_dependency(self, source_experiment, target_experiment, source_builder, target_builder):
        source_bobj = self.builder_meta[source_builder]
        target_bobj = self.builder_meta[target_builder]
        
        # get machine from builder. Note that builders must be unique across machines
        source_machine = source_bobj.machine
        target_machine = target_bobj.machine

        # get the experiments involved
        source_eobj = self.data.loc[source_experiment, (source_machine, source_builder)]
        target_eobj = self.data.loc[target_experiment, (target_machine, target_builder)]

        if(not isinstance(source_eobj, Experiment)):
            print("dependency source experiment {} does not exist on builder {}".format(source_experiment, source_builder)) 
        elif(not isinstance(target_eobj, Experiment)):
            print("dependency target experiment {} does not exist on builder {}".format(target_experiment, target_builder)) 
        else:
            # connect experiments by adding to respective parents/children
            source_eobj.add_parent(target_eobj)
            target_eobj.add_child(source_eobj)

    def run_builder(self, builder):
        self.get_experiments_by_builder(builder)


    def to_string(self):
        machines = self.data.columns.levels[0]

        col_idx = self.data.columns.to_frame()

        out = ""
        for midx in machines:
            out += "{} ():\n".format(midx)
            builders = col_idx.loc[(midx, slice(None)), "builder"].sort_index().values
            for bidx in builders:
                bobj = self.builder_meta[bidx]
                out += "  {}\n".format(bobj.to_string())
                if bobj.flag.lower() == "active":
                    experiments = self.data.loc[:, (midx, bidx)].sort_index().values
                    for exp in experiments:
                        if isinstance(exp, Experiment):
                            out += "    {}\n".format(exp.to_string())

        return out

    def list_builders(self):
        machines = self.data.columns.levels[0]

        col_idx = self.data.columns.to_frame()

        builderNames = []
        for midx in machines:
            builders = col_idx.loc[(midx, slice(None)), "builder"].sort_index().values
            for b in builders:
                builderNames.append(b)

        return ' '.join(builderNames)
