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

import re

from depgen import (
    IncludeFinder,
    StreamStack,
    file_in_dir,
    open23,
    find_unquoted_string,
)


class FortranParser:
    _re_include = re.compile(r'^\s*include\s+(\'|")(.*?)\1', re.I)
    _re_line_continue_start = re.compile(r"^(.*)&\s*$")
    _re_line_continue_end = re.compile(r"^\s*&")
    _re_module_provide = re.compile(r"^\s*module\s+(?!procedure\s)(\w+)", re.I)
    _re_module_require = re.compile(
        r"^\s*use(?:\s+|(?:\s*,\s*((?:non_)?intrinsic))?\s*::\s*)(\w+)", re.I
    )

    def __init__(
        self,
        include_order=None,
        include_dirs=None,
        include_roots=None,
        intrinsic_mods=None,
        external_mods=None,
    ):
        self.include_roots = include_roots

        if intrinsic_mods:
            self.intrinsic_mods = set(intrinsic_mods)
        else:
            self.intrinsic_mods = set()

        if external_mods:
            self.external_mods = set(external_mods)
        else:
            self.external_mods = set()

        # Callbacks:
        self.include_callback = None
        self.module_callback = None
        self.use_module_callback = None
        self.debug_callback = None

        self._include_finder = IncludeFinder(include_order, include_dirs)

    def parse(self, stream):
        with StreamStack() as s:
            s.add(stream)
            while 1:
                line = s.readline()
                if not line:
                    break

                # delete comments
                line = FortranParser._delete_comments(line)
                if line.isspace():
                    continue

                # line continuation
                match = FortranParser._re_line_continue_start.match(line)
                while match:
                    next_line = s.readline()
                    if not next_line:
                        break

                    next_line = FortranParser._delete_comments(next_line)

                    # If the line contains only comments, we need the next one
                    # TODO: implement a separate class FortranPrepcocessor
                    if next_line.isspace():
                        continue

                    line = match.group(1) + re.sub(
                        FortranParser._re_line_continue_end, "", next_line
                    )

                    match = FortranParser._re_line_continue_start.match(line)

                for line in FortranParser._split_semicolons(line):
                    # module provided
                    match = FortranParser._re_module_provide.match(line)
                    if match:
                        module_name = match.group(1).lower()
                        if self.module_callback:
                            self.module_callback(module_name)
                        if self.debug_callback:
                            self.debug_callback(
                                line, "declared module '%s'" % module_name
                            )
                        continue

                    # module required
                    match = FortranParser._re_module_require.match(line)
                    if match:
                        module_nature = (
                            match.group(1).lower()
                            if match.group(1) is not None
                            else ""
                        )
                        module_name = match.group(2).lower()
                        if module_nature == "intrinsic":
                            if self.debug_callback:
                                self.debug_callback(
                                    line,
                                    "ignored module usage ('%s' "
                                    "is explicitly intrinsic)" % module_name,
                                )
                        elif (
                            module_name in self.intrinsic_mods
                            and module_nature != "non_intrinsic"
                        ):
                            if self.debug_callback:
                                self.debug_callback(
                                    line,
                                    "ignored module usage ('%s' "
                                    "is implicitly intrinsic)" % module_name,
                                )
                        elif module_name in self.external_mods:
                            if self.debug_callback:
                                self.debug_callback(
                                    line,
                                    "ignored module usage ('%s' "
                                    "is external)" % module_name,
                                )
                        else:
                            if self.use_module_callback:
                                self.use_module_callback(module_name)
                            if self.debug_callback:
                                self.debug_callback(
                                    line, "used module '%s'" % module_name
                                )
                        continue

                    # include statement
                    match = FortranParser._re_include.match(line)
                    if match:
                        filename = match.group(2)
                        filepath = self._include_finder.find(
                            filename, s.root_name, s.current_name
                        )
                        if filepath:
                            if not self.include_roots or any(
                                [
                                    file_in_dir(filepath, d)
                                    for d in self.include_roots
                                ]
                            ):
                                s.add(open23(filepath, "r"))
                                if self.include_callback:
                                    self.include_callback(filepath)
                                if self.debug_callback:
                                    self.debug_callback(
                                        line, "included file '%s'" % filepath
                                    )
                            elif self.debug_callback:
                                self.debug_callback(
                                    line,
                                    "ignored (file '%s' "
                                    "is not in the source roots)" % filepath,
                                )
                        elif self.debug_callback:
                            self.debug_callback(
                                line, "ignored (file not found)"
                            )
                        continue

    @staticmethod
    def _split_semicolons(line):
        while 1:
            idx = find_unquoted_string(";", line)
            if idx < 0:
                if line and not line.isspace():
                    yield line
                break
            else:
                prefix = line[:idx]
                if prefix and not prefix.isspace():
                    yield prefix + "\n"
                line = line[idx + 1 :]

    @staticmethod
    def _delete_comments(line):
        comment_idx = find_unquoted_string("!", line)
        if comment_idx >= 0:
            line = line[:comment_idx]
        return line
