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
    find_unquoted_string,
    open23,
)


class Parser:
    _re_include = re.compile(r'^\s*include\s+([\'"])(.*?)\1\s*$', re.I)
    _re_line_continue_start = re.compile(r"^(.*)&\s*$")
    _re_line_continue_end = re.compile(r"^\s*&")
    _re_module_start = re.compile(
        r"^\s*module\s+(?!(?:procedure|subroutine|function)\s)(\w+)\s*$", re.I
    )
    _re_submodule_start = re.compile(
        r"^\s*submodule\s+\(\s*(\w+)(?:\s*:\s*(\w+)\s*)?\s*\)\s+(\w+)\s*$", re.I
    )
    _re_module_use = re.compile(
        r"^\s*use(?:\s+|(?:\s*,\s*((?:non_)?intrinsic))?\s*::\s*)(\w+)", re.I
    )
    _re_module_end = re.compile(r"^\s*end\s+module(?:\s+(\w+))?\s*$", re.I)
    _re_module_prefixed_procedure = re.compile(
        r"^\s*(?:\w+(?:\(.*\))?\s+)*"
        r"module\s+"
        r"(?:\w+(?:\(.*\))?\s+)*(?:function|subroutine)",
        re.I,
    )

    def __init__(
        self,
        include_order=None,
        include_dirs=None,
        include_roots=None,
        intrinsic_mods=None,
        external_mods=None,
        subparser=None,
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

        self._get_stream_iterator = (
            subparser.parse if subparser else lambda x, *_: x
        )

        # Callbacks:
        self.include_callback = None
        self.module_start_callback = None
        self.submodule_start_callback = None
        self.module_use_callback = None
        self.extendable_module_callback = None
        self.debug_callback = None

        self._include_finder = IncludeFinder(include_order, include_dirs)

    def parse(self, stream, stream_name):
        stream = self._get_stream_iterator(stream, stream_name)

        include_stack = StreamStack()
        include_stack.push(stream, stream_name)

        current_module = None
        current_module_is_extendable = False

        for line in Parser.streamline_input(include_stack):
            # module definition start
            match = Parser._re_module_start.match(line)
            if match:
                module_name = match.group(1).lower()
                current_module = module_name
                if self.module_start_callback:
                    self.module_start_callback(module_name)
                if self.debug_callback:
                    self.debug_callback(
                        line, "module '{0}' (start)".format(module_name)
                    )
                continue

            # submodule definition start
            match = Parser._re_submodule_start.match(line)
            if match:
                module_name = match.group(1).lower()
                parent_name = match.group(2)
                if parent_name:
                    parent_name = parent_name.lower()
                submodule_name = match.group(3).lower()
                if self.submodule_start_callback:
                    self.submodule_start_callback(
                        submodule_name, parent_name, module_name
                    )
                if self.debug_callback:
                    self.debug_callback(
                        line,
                        "submodule '{0}'{1} of module '{2}' (start)".format(
                            submodule_name,
                            (
                                " with parent '{0}'".format(parent_name)
                                if parent_name
                                else ""
                            ),
                            module_name,
                        ),
                    )
                continue

            # module used
            match = Parser._re_module_use.match(line)
            if match:
                module_nature = match.group(1)
                if module_nature:
                    module_nature = module_nature.lower()
                module_name = match.group(2).lower()
                if module_nature == "intrinsic":
                    if self.debug_callback:
                        self.debug_callback(
                            line,
                            "ignored module usage "
                            "('{0}' is explicitly intrinsic)".format(
                                module_name
                            ),
                        )
                elif (
                    module_name in self.intrinsic_mods
                    and module_nature != "non_intrinsic"
                ):
                    if self.debug_callback:
                        self.debug_callback(
                            line,
                            "ignored module usage "
                            "('{0}' is implicitly intrinsic)".format(
                                module_name
                            ),
                        )
                elif module_name in self.external_mods:
                    if self.debug_callback:
                        self.debug_callback(
                            line,
                            "ignored module usage "
                            "('{0}' is external)".format(module_name),
                        )
                else:
                    if self.module_use_callback:
                        self.module_use_callback(module_name)
                    if self.debug_callback:
                        self.debug_callback(
                            line,
                            "used module '{0}'".format(module_name),
                        )
                continue

            # include statement
            match = Parser._re_include.match(line)
            if match:
                filename = match.group(2)
                filepath = self._include_finder.find(
                    filename,
                    include_stack.root_name,
                    include_stack.current_name,
                )
                if filepath:
                    if not self.include_roots or any(
                        [file_in_dir(filepath, d) for d in self.include_roots]
                    ):
                        include_stack.push(open23(filepath, "r"), filepath)
                        if self.include_callback:
                            self.include_callback(filepath)
                        if self.debug_callback:
                            self.debug_callback(
                                line,
                                "included file '{0}'".format(filepath),
                            )
                    elif self.debug_callback:
                        self.debug_callback(
                            line,
                            "ignored (file '{0}' is not "
                            "in the source roots)".format(filepath),
                        )
                elif self.debug_callback:
                    self.debug_callback(line, "ignored (file not found)")
                continue

            if not (self.extendable_module_callback or self.debug_callback):
                continue

            # procedure with the module prefix
            match = Parser._re_module_prefixed_procedure.match(line)
            if match:
                if current_module:
                    if current_module_is_extendable:
                        if self.debug_callback:
                            self.debug_callback(
                                line,
                                "ignored module subroutine/function "
                                "(module '{0}' is already known "
                                "to be extendable)".format(current_module),
                            )
                    else:
                        current_module_is_extendable = True
                        if self.extendable_module_callback:
                            self.extendable_module_callback(current_module)
                        if self.debug_callback:
                            self.debug_callback(
                                line,
                                "module subroutine/function "
                                "(module '{0}' is extendable)".format(
                                    current_module
                                ),
                            )
                elif self.debug_callback:
                    self.debug_callback(
                        line,
                        "ignored module subroutine/function "
                        "(not in the module scope)".format(current_module),
                    )
                continue

            # module definition end
            match = Parser._re_module_end.match(line)
            if match:
                if self.debug_callback:
                    module_name = match.group(1)
                    if module_name:
                        module_name = module_name.lower()
                    self.debug_callback(
                        line,
                        "module '{0}' (end, {1})".format(
                            str(module_name),
                            (
                                (
                                    (
                                        "as expected"
                                        if module_name == current_module
                                        else "expected '{0}'".format(
                                            current_module
                                        )
                                    )
                                    if module_name
                                    else "assumed '{0}'".format(current_module)
                                )
                                if current_module
                                else "unexpected"
                            ),
                        ),
                    )
                current_module = None
                current_module_is_extendable = False
                continue

        include_stack.clear()

        # return an empty iterator
        return
        # noinspection PyUnreachableCode
        yield

    @staticmethod
    def streamline_input(stream):
        stream = Parser.drop_comments_and_empty_lines(stream)
        for line in stream:
            # concatenate lines
            while 1:
                match = Parser._re_line_continue_start.match(line)
                if not match:
                    break

                next_line = next(stream, None)
                if next_line is None:
                    break

                line = match.group(1) + re.sub(
                    Parser._re_line_continue_end, "", next_line
                )

            # split semicolons
            while 1:
                idx = find_unquoted_string(";", line)
                if idx < 0:
                    break

                yield line[:idx] + "\n"
                line = line[idx + 1 :]

            yield line

    @staticmethod
    def drop_comments_and_empty_lines(stream):
        for line in stream:
            comment_idx = find_unquoted_string("!", line)
            if comment_idx >= 0:
                line = line[:comment_idx]
            if line and not line.isspace():
                yield line
