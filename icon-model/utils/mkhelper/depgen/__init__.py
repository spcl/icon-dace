# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import os
import sys

try:
    from itertools import zip_longest as zip_longest23
except ImportError:
    # noinspection PyUnresolvedReferences
    from itertools import izip_longest as zip_longest23


def open23(name, mode="r"):
    if sys.version_info < (3, 0, 0):
        return open(name, mode)
    else:
        # noinspection PyArgumentList
        return open(
            name,
            mode,
            encoding=(None if "b" in mode else "UTF-8"),
            errors=(None if "b" in mode else "surrogateescape"),
        )


def map23(foo, iterable):
    if sys.version_info < (3, 0, 0):
        return map(foo, iterable)
    else:
        return list(map(foo, iterable))


if hasattr(sys, "implementation") and sys.implementation.name == "cpython":
    # see https://docs.python.org/3/library/itertools.html#itertools-recipes
    import collections

    def exhaust(it):
        return collections.deque(it, maxlen=0)

else:

    def exhaust(it):
        for _ in it:
            pass


def file_in_dir(f, d):
    if d:
        return os.path.abspath(f).startswith(os.path.abspath(d) + os.path.sep)
    else:
        return True


def find_unquoted_string(string, line, quotes="'\""):
    skip = 0
    quote = None
    while 1:
        idx = line.find(string, skip)
        if idx < 0:
            return idx

        escaped = False
        for c in line[skip:idx]:
            if escaped:
                escaped = False
            elif c in quotes:
                if quote is None:
                    quote = c
                elif quote == c:
                    quote = None
            elif c == "\\" and quote:
                escaped = True

        if quote:
            skip = idx + len(string)
        else:
            return idx


class IncludeFinder:
    def __init__(self, include_order=None, include_dirs=None):
        self.include_order = include_order
        self.include_dirs = include_dirs

    def find(self, filename, root_includer=None, current_includer=None):
        if os.path.isabs(filename) and os.path.isfile(filename):
            return filename
        elif self.include_order:
            for inc_type in self.include_order:
                if inc_type == "cwd" and os.path.isfile(filename):
                    return filename
                elif inc_type == "src" and root_includer:
                    candidate = os.path.join(
                        os.path.dirname(root_includer), filename
                    )
                    if os.path.isfile(candidate):
                        return candidate
                elif inc_type == "inc" and current_includer:
                    candidate = os.path.join(
                        os.path.dirname(current_includer), filename
                    )
                    if os.path.isfile(candidate):
                        return candidate
                elif inc_type == "flg" and self.include_dirs:
                    for d in self.include_dirs:
                        candidate = os.path.join(d, filename)
                        if os.path.isfile(candidate):
                            return candidate
        return None


class StreamStack(object):
    __slots__ = ["_stream_stack", "_close_stack", "_name_stack"]

    def __init__(self):
        # Stack of file-like objects (i.e. string iterators with the close
        # method:
        self._stream_stack = []
        self._close_stack = []
        self._name_stack = []

    @property
    def root_name(self):
        return self._name_stack[0] if self._name_stack else None

    @property
    def current_name(self):
        return self._name_stack[-1] if self._name_stack else None

    def push(self, stream, name=None, close=True):
        self._stream_stack.append(stream)
        self._close_stack.append(close)
        self._name_stack.append(name)

    def clear(self):
        for stream, close in zip(self._stream_stack, self._close_stack):
            if close:
                stream.close()
        self._stream_stack *= 0
        self._close_stack *= 0
        self._name_stack *= 0

    def __iter__(self):
        return self

    def __next__(self):
        while self._stream_stack:
            try:
                return next(self._stream_stack[-1])
            except StopIteration:
                self._name_stack.pop()
                stream = self._stream_stack.pop()
                if self._close_stack.pop():
                    stream.close()
        raise StopIteration

    if sys.version_info < (3,):
        next = __next__
