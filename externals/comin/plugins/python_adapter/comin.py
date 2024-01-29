# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

from _comin import *
import _comin
from dataclasses import dataclass
import numpy as _np


def register_callback(ep):
    def __callback(fun):
        _comin._callback_register(ep, fun)
        return fun
    return __callback


COMIN_ZAXIS_NONE  = 0
COMIN_ZAXIS_2D    = 1
COMIN_ZAXIS_3D    = 2
COMIN_ZAXIS_UNDEF = 3


class _variable:
    def __init__(self, handle):
        self._handle = handle

    def __array__(self):
        return _np.asarray(_comin._var_get_buffer(self._handle))

    @property
    def pos(self):
        return _comin._var_get_pos(self._handle)

    @property
    def ncontained(self):
        return _comin._var_get_ncontained(self._handle)

    @property
    def to_3d(self):
        missing_dims = {0, 1, 2, 3, 4}.difference({*self.pos})
        if self.ncontained > 0:
            return _np.asarray(self).transpose(*self.pos, *missing_dims)[..., 0, 0]
        else:
            return _np.asarray(self).transpose(*self.pos[0:3], *missing_dims)[..., 0, 0]


def var_get(context, var_descriptor, flag=-1):
    """get variable object, arguments: entry point, name string, domain id, access flag)"""
    return _variable(_comin._var_get(context, var_descriptor, flag))


for ep in range(1, _comin._EP_DESTRUCTOR()+1):
    name = _comin.callback_get_ep_name(ep)
    vars()[name] = ep


@dataclass
class plugin_info:
    id: int
    name: str
    options: str
    comm: str


def current_get_plugin_info():
    """returns object describing the current plugin"""
    return plugin_info(**_comin._current_get_plugin_info())


class _descrdata:
    def __init__(self, properties, jg=0):
        self.properties = properties
        self.jg = jg

    def __dir__(self):
        return self.properties.keys()

    def __getattr__(self, key):
        val = _comin._descrdata_eval_property(self.properties[key], jg=self.jg)
        if isinstance(val, dict):
            return _descrdata(val, jg=self.jg)
        else:
            return val


def descrdata_get_domain(jg):
    """returns descriptive data for a given domain, arguments: jg"""
    return _descrdata(_comin._descrdata_get_domain(), jg=jg)

def descrdata_get_global():
    """returns global descriptive data object"""
    return _descrdata(_comin._descrdata_get_global())

def var_descr_list():
    """List of exposed variables (descriptors)"""
    current = _comin._var_get_descr_list_head()
    while current is not None:
        yield _comin._var_get_descr_list_var_desc(current)
        current = _comin._var_get_descr_list_next(current)


def metadata_set(var_descriptor, **kwargs):
    """sets metadata for a requested field, arguments: name string, domain id, metadata key, metadata value"""
    for n, v in kwargs.items():
        _comin._metadata_set(var_descriptor, n, v)


@dataclass
class simulation_interval:
    exp_start : str
    exp_stop  : str
    run_start : str
    run_stop  : str


def descrdata_get_simulation_interval():
    """"returns simulation intervals: exp_start, exp_stop, run_start, run_stop"""
    return simulation_interval(**_comin._descrdata_get_simulation_interval())


COMIN_FLAG_NONE = 0
COMIN_FLAG_READ = 1 << 1
COMIN_FLAG_WRITE = 1 << 2

COMIN_HGRID_UNSTRUCTURED_CELL = 1
