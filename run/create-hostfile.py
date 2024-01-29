#! /usr/bin/env python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#______________________________________________________________________________

import subprocess
import argparse

parser = argparse.ArgumentParser(prog='create_hostfile',
                                 description='Create hostfile')

parser.add_argument('--tasks-per-node', type=int, default=1,
                    help='number of tasks per node', required=True)
parser.add_argument('--io-tasks', type=int, default=0,
                    help='number of separate output tasks', required=True)

args = parser.parse_args()

#______________________________________________________________________________
#
def nodenames():
    #return ["m1", "m2", "m3" ]
    return subprocess.run(['scontrol', 'show', 'hostnames'], stdout=subprocess.PIPE).stdout.decode('UTF-8').splitlines()
#______________________________________________________________________________
#

tasks_per_node = args.tasks_per_node
io_tasks       = args.io_tasks

hostlist = nodenames()
nodes=len(hostlist)

io_nodes = (io_tasks + tasks_per_node - 1) // tasks_per_node
compute_nodes = nodes - io_nodes

for hostname in hostlist[0:compute_nodes]:
    for task in range(tasks_per_node):
        print(hostname)        
        # print("Compute ",hostname)
        
last_node = (nodes-1, nodes) [(io_tasks % tasks_per_node) == 0]
    
for hostname in hostlist[compute_nodes:last_node]:
    for task in range(tasks_per_node):
        print(hostname)        
        # print("I/O     ", hostname)

for r in range(io_tasks % tasks_per_node):
    print(hostlist[-1])        
    # print("I/O     ", hostlist[-1])
