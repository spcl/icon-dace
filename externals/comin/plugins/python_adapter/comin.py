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
import shlex as _shlex
import sys
from collections.abc import Mapping


def register_callback(ep):
    def __callback(fun):
        _comin._callback_register(ep, fun)
        return fun
    return __callback


COMIN_ZAXIS_UNDEF   = -1
COMIN_ZAXIS_NONE    =  0
COMIN_ZAXIS_2D      =  1
COMIN_ZAXIS_3D      =  2
COMIN_ZAXIS_3D_HALF =  3


class _variable:
    def __init__(self, handle):
        self._handle = handle
        try:
            import numpy as _np
            self.np = _np
        except ImportError:
            print("Warning: cant import numpy", file=sys.stderr)

    def __array__(self):
        return self.np.asarray(_comin._var_get_buffer(self._handle))

    @property
    def __cuda_array_interface__(self):
        import numpy as _np
        host_buf = _comin._var_get_buffer(self._handle)
        return {
            "shape": host_buf.shape,
            "typestr": _np.dtype(host_buf.format).str,
            "data": (_comin._var_get_device_ptr(self._handle), False),
            "version": 3,
            "strides": host_buf.strides
        }

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
            return self.np.asarray(self).transpose(*self.pos, *missing_dims)[..., 0, 0]
        else:
            return self.np.asarray(self).transpose(*self.pos[0:3], *missing_dims)[..., 0, 0]

    @property
    def descriptor(self):
        return _comin._var_get_descriptor(self._handle)


def var_get(context, var_descriptor, flag):
    """get variable object, arguments: [entry point], (name string, domain id), access flag)"""
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

    @property
    def args(self):
        """
        Extract the argument from the options string like as the
        script was called from the command line. This is supposed to
        be passed to `argparse.ArgumentParser.parse_args`.
        """
        return _shlex.split(self.options)[1:]


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


class metadata(Mapping):
    def __init__(self, var_descr):
        self.descr = var_descr

    def __getitem__(self, key):
        return metadata_get(self.descr, key)

    def __iter__(self):
        it = _comin._metadata_get_iterator_begin(self.descr)
        end = _comin._metadata_get_iterator_end(self.descr)
        while  not _comin._metadata_iterator_compare(it, end):
            yield _comin._metadata_iterator_get_key(it)
            _comin._metadata_iterator_next(it)
        _comin._metadata_iterator_delete(it)
        _comin._metadata_iterator_delete(end)

    def __len__(self):
        # impl. could be improved
        return len(self.__iter__)

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
COMIN_FLAG_DEVICE = 1 << 4

COMIN_HGRID_UNSTRUCTURED_CELL   = 1
COMIN_HGRID_UNSTRUCTURED_EDGE   = 2
COMIN_HGRID_UNSTRUCTURED_VERTEX = 3
