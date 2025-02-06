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
    _re_ifdef = re.compile(r"^\s*#\s*if(n)?def\s+([a-zA-Z_]\w*)")
    _re_if_expr = re.compile(r"^\s*#\s*if((?:\s|\().*)")

    _re_elif = re.compile(r"^\s*#\s*elif((?:\s|\().*)")
    _re_else = re.compile(r"^\s*#\s*else(?:\s.*)")
    _re_endif = re.compile(r"^\s*#\s*endif(?:\s.*)")

    _re_include = re.compile(r'^\s*#\s*include\s+(?:"(.*?)"|<(.*?)>)')
    _re_define = re.compile(r"^\s*#\s*define\s+([a-zA-Z_]\w*)(\([^)]*\))?(.*)$")
    _re_undef = re.compile(r"^\s*#\s*undef\s+([a-zA-Z_]\w*)")

    def __init__(
        self,
        include_order=None,
        include_sys_order=None,
        include_dirs=None,
        include_roots=None,
        try_eval_expr=False,
        inc_sys=False,
        predefined_macros=None,
        subparser=None,
    ):
        self.include_roots = include_roots
        self.try_eval_expr = try_eval_expr
        self.inc_sys = inc_sys

        # Callbacks:
        self.include_callback = None
        self.debug_callback = None

        self._include_finder = IncludeFinder(include_order, include_dirs)
        self._include_sys_finder = IncludeFinder(
            include_sys_order, include_dirs
        )

        self._predefined_macros = predefined_macros

        self._get_stream_iterator = (
            subparser.parse if subparser else lambda x, *_: x
        )

    def parse(self, stream, stream_name):
        stream = self._get_stream_iterator(stream, stream_name)

        include_stack = StreamStack()
        include_stack.push(stream, stream_name)

        branch_state = BranchState()
        macro_handler = MacroHandler(self._predefined_macros)

        for line in Parser.streamline_input(include_stack):
            # if(n)def directive
            match = Parser._re_ifdef.match(line)
            if match:
                macro, negate, state = match.group(2), bool(match.group(1)), 0
                if not branch_state.is_dead():
                    state = macro_handler.eval_defined(macro, negate)
                    if self.debug_callback:
                        self.debug_callback(
                            line, "evaluated to {0}".format(state > 0)
                        )
                elif self.debug_callback:
                    self.debug_callback(line, "was not evaluated (dead branch)")
                branch_state.switch_if(state)
                continue

            # if directive
            match = Parser._re_if_expr.match(line)
            if match:
                expr, state = match.group(1), 0
                if not branch_state.is_dead():
                    if self.try_eval_expr:
                        state = macro_handler.eval_expression(expr)
                        if self.debug_callback:
                            self.debug_callback(
                                line,
                                "evaluated to {0}".format(
                                    "Unknown (evaluation failed)"
                                    if state == 0
                                    else state > 0
                                ),
                            )
                    elif self.debug_callback:
                        self.debug_callback(
                            line, "was not evaluated (evaluation disabled)"
                        )
                elif self.debug_callback:
                    self.debug_callback(line, "was not evaluated (dead branch)")
                branch_state.switch_if(state)
                continue

            # elif directive
            match = Parser._re_elif.match(line)
            if match:
                branch_state.switch_else()
                expr, state = match.group(1), 0
                if not branch_state.is_dead():
                    if self.try_eval_expr:
                        state = macro_handler.eval_expression(expr)
                        if self.debug_callback:
                            self.debug_callback(
                                line,
                                "evaluated to {0}".format(
                                    "Unknown (evaluation failed)"
                                    if state == 0
                                    else state > 0
                                ),
                            )
                    elif self.debug_callback:
                        self.debug_callback(
                            line, "was not evaluated (evaluation disabled)"
                        )
                elif self.debug_callback:
                    self.debug_callback(line, "was not evaluated (dead branch)")
                branch_state.switch_elif(state)
                continue

            # else directive
            match = Parser._re_else.match(line)
            if match:
                branch_state.switch_else()
                continue

            # endif directive
            match = Parser._re_endif.match(line)
            if match:
                branch_state.switch_endif()
                continue

            if branch_state.is_dead() and self.debug_callback is None:
                continue

            # define directive
            match = Parser._re_define.match(line)
            if match:
                if not branch_state.is_dead():
                    macro_handler.define(*match.group(1, 2, 3))
                    if self.debug_callback:
                        self.debug_callback(line, "accepted")
                elif self.debug_callback:
                    self.debug_callback(line, "ignored (dead branch)")
                continue

            # undef directive
            match = Parser._re_undef.match(line)
            if match:
                if not branch_state.is_dead():
                    macro_handler.undefine(match.group(1))
                    if self.debug_callback:
                        self.debug_callback(line, "accepted")
                elif self.debug_callback:
                    self.debug_callback(line, "ignored (dead branch)")
                continue

            # include directive
            match = Parser._re_include.match(line)
            if match:
                if not branch_state.is_dead():
                    if match.lastindex == 1:  # quoted form
                        filepath = self._include_finder.find(
                            match.group(1),
                            include_stack.root_name,
                            include_stack.current_name,
                        )
                    elif match.lastindex == 2:  # angle-bracket form
                        if self.inc_sys:
                            filepath = self._include_sys_finder.find(
                                match.group(2),
                                include_stack.root_name,
                                include_stack.current_name,
                            )
                        else:
                            if self.debug_callback:
                                self.debug_callback(
                                    line, "ignored (system header)"
                                )
                            continue
                    else:
                        if self.debug_callback:
                            self.debug_callback(
                                line, "ignored (internal error)"
                            )
                        continue

                    if filepath:
                        if not self.include_roots or any(
                            [
                                file_in_dir(filepath, d)
                                for d in self.include_roots
                            ]
                        ):
                            include_stack.push(open23(filepath, "r"), filepath)
                            if self.include_callback:
                                self.include_callback(filepath)
                            if self.debug_callback:
                                self.debug_callback(
                                    line, "included file '{0}'".format(filepath)
                                )
                        elif self.debug_callback:
                            self.debug_callback(
                                line,
                                "ignored (file '{0}' "
                                "is not in the source roots)".format(filepath),
                            )
                    elif self.debug_callback:
                        self.debug_callback(line, "ignored (file not found)")
                elif self.debug_callback:
                    self.debug_callback(line, "ignored (dead branch)")
                continue

            if branch_state.is_dead():
                continue

            yield line

    @staticmethod
    def streamline_input(stream):
        for line in stream:
            # concatenate lines
            while line.endswith("\\\n"):
                line = line[:-2] + next(stream, "")

            # remove block comments
            while 1:
                # Check whether the line contains an unquoted block comment
                # initiator '/*':
                start_idx = find_unquoted_string("/*", line)
                if start_idx < 0:
                    break
                # Check whether the line contains a block comment
                # terminator '*/' (even if it is quoted):
                term_idx = line.find("*/", start_idx + 2)
                while term_idx < 0:
                    # The block is not terminated yet, read the next line:
                    term_idx = len(line)
                    try:
                        line += next(stream)
                        term_idx = line.find("*/", term_idx)
                    except StopIteration:
                        pass
                else:
                    # Replace the block of comments with a single
                    # space:
                    line = "{0} {1}".format(
                        line[:start_idx], line[term_idx + 2 :]
                    )

            if not line or line.isspace():
                continue

            yield line


class BranchState(object):
    __slots__ = ["_if_state_stack", "_states_per_endif_stack"]

    def __init__(self):
        # Stack of #if-#else blocks holds one of the following:
        # 1 - keep current branch and ignore another
        # -1 - ignore current branch and keep another
        # 0 - keep both branches (if failed to evaluate expression)
        self._if_state_stack = []

        # Each #elif is interpreted as a combination of #else and #if.
        # Thus, each #elif increments the number of #if blocks that are
        # closed with #endif statement. This numbers are stored in a separate
        # stack:
        self._states_per_endif_stack = []

    def switch_if(self, state):
        self._if_state_stack.append(state)
        self._states_per_endif_stack.append(1)

    def switch_else(self):
        if self._if_state_stack:
            self._if_state_stack[-1] = -self._if_state_stack[-1]

    def switch_elif(self, state):
        self._if_state_stack.append(state)
        self._states_per_endif_stack[-1] += 1

    def switch_endif(self):
        if self._if_state_stack:
            pop_count = self._states_per_endif_stack.pop()
            for _ in range(pop_count):
                self._if_state_stack.pop()

    def is_dead(self):
        return any(state < 0 for state in self._if_state_stack)


class MacroHandler(object):
    # matches "defined MACRO_NAME" and "defined (MACRO_NAME)"
    _re_defined_call = re.compile(
        r"(defined\s*(\(\s*)?([a-zA-Z_]\w*)(?(2)\s*\)))"
    )

    # matches object-like and function-like macro identifiers
    _re_identifier = re.compile(
        r"(([a-zA-Z_]\w*)\s*(\(\s*(?:\w+(?:\s*,\s*\w+)*\s*)?\))?)"
    )

    __slots__ = ["_macros"]

    def __init__(self, predefined_macros=None):
        self._macros = dict(predefined_macros or [])

    def define(self, macro_name, macro_args=None, macro_body=None):
        if macro_name != "defined":
            self._macros[macro_name] = (macro_args, macro_body or "")

    def undefine(self, macro_name):
        self._macros.pop(macro_name, None)

    def eval_defined(self, macro_name, negate=False):
        return 1 if bool(macro_name in self._macros) ^ negate else -1

    def eval_expression(self, expr):
        prev_expr = None
        while expr != prev_expr:
            # replace calls to function "defined"
            defined_calls = re.findall(MacroHandler._re_defined_call, expr)
            for call in defined_calls:
                expr = expr.replace(
                    call[0], "1" if call[2] in self._macros else "0"
                )

            identifiers = re.findall(MacroHandler._re_identifier, expr)

            for identifier, identifier_name, identifier_args in identifiers:

                macro = self._macros.get(identifier_name, None)

                if identifier_args:
                    # expansion of a function-like identifier

                    if macro is None:
                        # the respective macro is not defined, which normally
                        # results in a preprocessor error
                        return 0
                    elif macro[0] is None:
                        # the respective macro is an object-like one: replace
                        # the name of the function-like identifier with the body
                        # (value) of the macro
                        identifier_repl = macro[1] + identifier_args
                    elif "".join(identifier_args.split()) == "()":
                        # the argument list of the function-like identifier is
                        # empty: replace the identifier with the body of the
                        # function-like macro
                        identifier_repl = macro[1]
                    else:
                        # we cannot expand function-like macros but we might not
                        # have to due to the short-circuit evaluation
                        continue
                else:
                    # expansion of an object-like identifier

                    if macro is None or macro[0] is not None:
                        # the respective macro is not defined or defined as a
                        # function-like one: the identifier evaluates to False
                        identifier_repl = "0"
                    else:
                        # the respective macro is defined as an object-like one:
                        # replace the identifier with the body (value) of the
                        # macro
                        identifier_repl = macro[1]

                expr = expr.replace(identifier, identifier_repl)

            prev_expr = expr

        expr = expr.replace("||", " or ")
        expr = expr.replace("&&", " and ")
        expr = expr.replace("!", "not ")

        try:
            result = bool(eval(expr, {"__builtins__": None}))
            return 1 if result else -1
        except Exception:
            return 0
