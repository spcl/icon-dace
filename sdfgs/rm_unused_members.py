import sys
import dace
import os
from dace import nodes, dtypes, memlet
import copy

velo_base_name = "velocity_tendencies_simplified_f.sdfgz"

print("LOAD BASE")
velo_tendencies = dace.SDFG.from_file(velo_base_name)

current_dir = os.path.dirname(os.path.abspath(__file__))

# Add the parent directory to the Python path
sys.path.append(current_dir)

from dace.sdfg.sdfg import InterstateEdge
from modules.map_over_tasklet import MapOverFreeTasklet
from modules.clean_unused_members import clean_unused_members

clean_unused_members(velo_tendencies)

velo_tendencies.save("sdfgs/velocity_tendencies_simplified_f_pruned.sdfgz")