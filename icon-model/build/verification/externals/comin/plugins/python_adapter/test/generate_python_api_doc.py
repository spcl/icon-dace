# This Python plugin generates a rudimentary markdown documentation of
# the available ComIn Python API.
#
# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

import comin
import inspect
import numpy as np

ignore = ["dataclass", "ep", "name", "register_callback", "simulation_interval", "plugin_info"]

def describe(obj,f,prefix=""):
    for element_name in dir(obj):
        if (not element_name in ignore):
            element = getattr(obj, element_name)
            if (element_name.startswith('_')):
                pass
            elif inspect.isclass(element):
                print(f"class {prefix}{element_name}", file=f)
            elif inspect.ismodule(element):
                print(f"ismodule {prefix}{element_name}", file=f)
                pass
            elif callable(element):
                #if not inspect.isbuiltin(element):
                print(f"- callable `{prefix}{element_name}`: {element.__doc__}", file=f)
                pass
            else:
                try:
                    dims = [":"] * element.ndim
                    print(f"- `{prefix}{element_name}({','.join(dims)})` ({np.asarray(element).dtype})", file=f)
                except:
                    print(f"- `{prefix}{element_name}` ({type(element).__name__})", file=f)

with open("comin_python_api.md", 'w') as f:
    print("# ComIn Python API", file=f)
    print("\n\n## Global ComIn functions, variables and constants", file=f)
    describe(comin,f)
    print("\n\n## Members of data type plugin_info", file=f)
    plugin_info = comin.current_get_plugin_info()
    describe(plugin_info,f,prefix="plugin_info.")
    print("\n\n## Global descriptive data", file=f)
    glob = comin.descrdata_get_global()
    describe(glob,f,prefix="global.")
    print("\n\n## Members of data type simulation_interval", file=f)
    simulation_interval = comin.descrdata_get_simulation_interval()
    describe(simulation_interval,f,prefix="simulation_interval.")
    print("\n\n## Descriptive data for domains", file=f)
    domain = comin.descrdata_get_domain(1)
    describe(domain,f,prefix="domain.")
    print("\n\n## Members of data type domain.cells", file=f)
    describe(domain.cells,f,prefix="domain.cells.")
    print("\n\n## Members of data type domain.edges", file=f)
    describe(domain.edges,f,prefix="domain.edges.")
    print("\n\n## Members of data type domain.verts", file=f)
    describe(domain.verts,f,prefix="domain.verts.")
