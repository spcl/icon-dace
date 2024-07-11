#!/bin/sh

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

""":"
for cmd in python3 python; do
  if command -v > /dev/null "${cmd}"; then
    exec "${cmd}" "$0" "$@"
  fi
done
echo "Error: could not find a python interpreter!" >&2
exit 1
":"""

import argparse
import collections
import fnmatch
import itertools
import os
import re
import sys

_re_rule = re.compile(
    r"^[ ]*(?P<targets>[^:|#\s]+(?:[ ]+[^:|#\s]+)*)[ ]*"
    r":[ ]*(?P<normal>[^:|#\s]+(?:[ ]+[^:|#\s]+)*)?[ ]*"
    r"\|?[ ]*(?P<order_only>[^:|#\s]+(?:[ ]+[^:|#\s]+)*)?[ ]*"
    r"(?:#-hint)?[ ]*(?P<hint>[^:|#\s]+(?:[ ]+[^:|#\s]+)*)?[ ]*"
)

_meta_root = 0
_term_colors = {
    "black": 90,
    "red": 91,
    "green": 92,
    "yellow": 93,
    "blue": 94,
    "magenta": 95,
    "cyan": 96,
    "white": 97,
}


def parse_args():
    class ArgumentParser(argparse.ArgumentParser):
        def convert_arg_line_to_args(self, arg_line):
            try:
                # Drop everything after the first occurrence of #:
                arg_line = arg_line[: arg_line.index("#")]
            except ValueError:
                pass

            result = []
            # Do not regard consecutive whitespaces as a single separator:
            for arg in arg_line.split(" "):
                if arg:
                    result.append(arg)
                elif result:
                    # The previous argument has a significant space:
                    result[-1] += " "
            return result

    parser = ArgumentParser(
        fromfile_prefix_chars="@",
        description="Reads a set of MAKEFILEs and prints a topologically "
        "sorted list of TARGETs (PREREQuisites) together with their "
        "dependencies (dependents).",
    )

    parser.add_argument(
        "-d", "--debug-file", help="dump debug information to DEBUG_FILE"
    )
    parser.add_argument(
        "-t",
        "--target",
        nargs="*",
        help="names of the makefile targets to be printed together with their "
        "dependencies; mutually exclusive with the argument '-p/--prereq'; if "
        "neither of the arguments is specified, all targets and prerequisites "
        "found in the makefiles are sent to the output",
    )
    parser.add_argument(
        "-p",
        "--prereq",
        nargs="*",
        help="names of the makefile prerequisites to be printed together with "
        "their dependents; mutually exclusive with the argument '-t/--target'; "
        "if neither of the arguments is specified, all targets and "
        "prerequisites found in the makefiles are sent to the output",
    )
    parser.add_argument(
        "--max-depth",
        metavar="MAX_DEPTH",
        type=int,
        help="print dependencies (dependents) that are at most MAX_DEPTH "
        "levels from the requested targets (prerequisites)",
    )
    parser.add_argument(
        "--ignore-order-only",
        "--no-oo",
        action="store_true",
        help="ignore order-only prerequisites",
    )
    parser.add_argument(
        "--ignore-hints",
        "--no-hints",
        action="store_true",
        help="ignore #-hint prerequisites",
    )
    parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        help="print the output list in the reversed order",
    )
    parser.add_argument(
        "--check-unique-prereq",
        action="append",
        # Unfortunately, we cannot set nargs to 'two or more', therefore we
        # set nargs to 'one or more':
        nargs="+",
        metavar="PATTERN",
        help="list of two or more shell-like wildcards; the option enables "
        "additional consistency checks of the dependency graph: each target "
        "that matches the first pattern of the list is checked whether it has "
        "no more than one prerequisite matching any of the rest of the "
        "patterns; if the check fails, a warning message is emitted to the "
        "standard error stream",
    )
    parser.add_argument(
        "--check-unique-basename",
        action="append",
        nargs="+",
        metavar="PATTERN",
        help="list of shell-like wildcards; the option enables additional "
        "consistency checks of the dependency graph; all targets that match at "
        "least one of the patterns are checked whether none of them have the "
        "same basename; if the check fails, a warning message is emitted to "
        "the standard error stream",
    )
    parser.add_argument(
        "--check-exists-prereq",
        action="append",
        # Unfortunately, we cannot set nargs to 'two or more', therefore we
        # set nargs to 'one or more':
        nargs="+",
        metavar="PATTERN",
        help="list of two or more shell-like wildcards; the option enables "
        "additional consistency checks of the dependency graph: each target "
        "that matches the first pattern of the list is checked whether it has "
        "at least one prerequisite matching any of the rest of the patterns; "
        "if the check fails, a warning message is emitted to the standard "
        "error stream",
    )
    parser.add_argument(
        "--check-cycles",
        action="store_true",
        help="check whether the dependency graph is acyclic, e.g. there is no "
        "circular dependencies; if a cycle is found, a warning message is "
        "emitted to the standard output",
    )
    parser.add_argument(
        "--check-colour",
        choices=_term_colors.keys(),
        help="colour the message output of the checks using ANSI escape "
        "sequences; the argument is ignored if the standard error stream is "
        "not associated with a terminal device",
    )
    parser.add_argument(
        "-f",
        "--makefile",
        nargs="*",
        help="paths to makefiles; a single dash (-) triggers reading from the "
        "standard input stream",
    )

    args = parser.parse_args()

    if args.target is not None and args.prereq is not None:
        parser.error(
            "arguments -t/--target and -p/--prereq are mutually " "exclusive"
        )

    if args.max_depth is not None and args.max_depth < 0:
        args.max_depth = None

    if args.check_unique_prereq:
        for pattern_list in args.check_unique_prereq:
            if len(pattern_list) < 2:
                parser.error(
                    "argument --check-unique-prereq: expected 2 or "
                    "more arguments"
                )

    if args.check_exists_prereq:
        for pattern_list in args.check_exists_prereq:
            if len(pattern_list) < 2:
                parser.error(
                    "argument --check-exists-prereq: expected 2 or "
                    "more arguments"
                )

    if not sys.stderr.isatty():
        args.check_colour = None

    return args


def read_makefiles(makefiles, ignore_order_only, ignore_hints):
    dep_graph = collections.defaultdict(list)
    extra_edges = collections.defaultdict(list)

    for mkf in makefiles:
        if mkf == "-":
            stream = sys.stdin
        elif not os.path.isfile(mkf):
            continue
        else:
            stream = open(mkf, "r")

        it = iter(stream)

        for line in it:
            while line.endswith("\\\n"):
                line = line[:-2]
                try:
                    line += next(it)
                except StopIteration:
                    break

            match = _re_rule.match(line)
            if match:
                targets = set(match.group("targets").split())
                prereqs = []

                prereqs_string = match.group("normal")
                if prereqs_string:
                    prereqs.extend(prereqs_string.split())

                if not ignore_order_only:
                    prereqs_string = match.group("order_only")
                    if prereqs_string:
                        prereqs.extend(prereqs_string.split())

                for target in targets:
                    dep_graph[target].extend(prereqs)

                if not ignore_hints:
                    prereqs_string = match.group("hint")
                    if prereqs_string:
                        for target in targets:
                            extra_edges[target].extend(prereqs_string.split())

        stream.close()
    return dep_graph, extra_edges


def visit_dfs(
    dep_graph,
    vertex,
    current_depth=0,
    max_depth=None,
    visited=None,
    start_visit_cb_list=None,
    finish_visit_cb_list=None,
    skip_visit_cb_list=None,
):
    if max_depth is not None and current_depth > max_depth:
        return

    if visited is None:
        visited = dict()

    if vertex in visited:
        if skip_visit_cb_list:
            for skip_visit_cb in skip_visit_cb_list:
                skip_visit_cb(vertex)
        return

    if start_visit_cb_list:
        for start_visit_cb in start_visit_cb_list:
            start_visit_cb(vertex)

    visited[vertex] = current_depth

    if vertex in dep_graph:
        for child in dep_graph[vertex]:
            visit_dfs(
                dep_graph,
                child,
                current_depth + 1,
                max_depth,
                visited,
                start_visit_cb_list,
                finish_visit_cb_list,
                skip_visit_cb_list,
            )

    if finish_visit_cb_list:
        for finish_visit_cb in finish_visit_cb_list:
            finish_visit_cb(vertex)


def dedupe(sequence):
    seen = set()
    for x in sequence:
        if x not in seen:
            yield x
            seen.add(x)


def sanitize_graph(graph):
    # Remove duplicates (we do not use sets as values of the dictionary to keep
    # the order of prerequisites):
    for target in graph.keys():
        graph[target] = graph.default_factory(
            dedupe(prereq for prereq in graph[target])
        )

    # Make leaves (i.e. prerequisites without any prerequisites) explicit nodes
    # of the graph:
    leaves = set(
        prereq
        for prereqs in graph.values()
        for prereq in prereqs
        if prereq not in graph
    )
    graph.update((prereq, graph.default_factory()) for prereq in leaves)


def flip_edges(graph):
    result = collections.defaultdict(list)
    for parent, children in graph.items():
        for child in children:
            result[child].append(parent)
        else:
            _ = result[parent]
    return result


def warn(msg, colour=None):
    sys.stderr.write(
        "{0}{1}: WARNING: {2}{3}\n".format(
            ("\033[{0}m".format(_term_colors[colour])) if colour else "",
            os.path.basename(__file__),
            msg,
            "\033[0m" if colour else "",
        )
    )


def main():
    args = parse_args()

    if args.debug_file:
        with open(args.debug_file, "w") as debug_file:
            debug_file.writelines(
                [
                    "# Python version: ",
                    sys.version.replace("\n", " "),
                    "\n",
                    "#\n",
                    "# Command:\n",
                    "#  ",
                    " ".join(sys.argv),
                    "\n",
                    "#\n",
                    "# Parsed arguments:\n",
                    "#  ",
                    "\n#  ".join(
                        [k + "=" + str(v) for k, v in vars(args).items()]
                    ),
                    "\n",
                ]
            )

    if args.makefile is None:
        return

    dep_graph, extra_edges = read_makefiles(
        args.makefile, args.ignore_order_only, args.ignore_hints
    )

    if not dep_graph:
        return

    sanitize_graph(dep_graph)
    sanitize_graph(extra_edges)

    if args.prereq is None:
        traversed_graph = dep_graph
        start_nodes = args.target
    else:
        traversed_graph = flip_edges(dep_graph)
        start_nodes = args.prereq
        extra_edges = flip_edges(extra_edges)

    # Insert _meta_root, which will be the starting-point for the dependency
    # graph traverse:
    if start_nodes is None:
        traversed_graph[_meta_root] = sorted(traversed_graph.keys())
    else:
        traversed_graph[_meta_root] = [
            t for t in start_nodes if t in traversed_graph
        ]

    # Visitor callbacks:
    start_visit_cb_list = []
    finish_visit_cb_list = []
    skip_visit_cb_list = []

    # Callbacks that are called once the graph is traversed:
    postprocess_cb_list = []

    if args.check_unique_prereq:

        def check_unique_prereq_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root or does not have descendants:
            if vertex == _meta_root:
                return
            for pattern_list in args.check_unique_prereq:
                if fnmatch.fnmatch(vertex, pattern_list[0]):
                    vertex_prereqs = dep_graph[vertex]
                    prereq_patterns = pattern_list[1:]
                    matching_prereqs = [
                        prereq
                        for prereq_pattern in prereq_patterns
                        for prereq in fnmatch.filter(
                            vertex_prereqs, prereq_pattern
                        )
                    ]
                    if len(matching_prereqs) > 1:
                        warn(
                            "target '{0}' has more than one immediate "
                            "prerequisite matching any of the patterns: "
                            "'{1}':\n\t{2}".format(
                                vertex,
                                "', '".join(prereq_patterns),
                                "\n\t".join(matching_prereqs),
                            ),
                            args.check_colour,
                        )

        start_visit_cb_list.append(check_unique_prereq_start_visit_cb)

    if args.check_unique_basename:
        basenames = [
            collections.defaultdict(set)
            for _ in range(len(args.check_unique_basename))
        ]

        def check_unique_basename_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root:
            if vertex == _meta_root:
                return
            for i, pattern_list in enumerate(args.check_unique_basename):
                for pattern in pattern_list:
                    if fnmatch.fnmatch(vertex, pattern):
                        basenames[i][os.path.basename(vertex)].add(vertex)

        start_visit_cb_list.append(check_unique_basename_start_visit_cb)

        def check_unique_basename_postprocess_cb():
            for basename_group in basenames:
                for basename, paths in basename_group.items():
                    if len(paths) > 1 and basename:
                        warn(
                            "the dependency graph contains more than one "
                            "target with basename '{0}':\n\t{1}".format(
                                basename, "\n\t".join(paths)
                            ),
                            args.check_colour,
                        )

        postprocess_cb_list.append(check_unique_basename_postprocess_cb)

    if args.check_exists_prereq:

        def check_exists_prereq_start_visit_cb(vertex):
            # Skip if the vertex is _meta_root:
            if vertex == _meta_root:
                return
            for pattern_list in args.check_exists_prereq:
                if fnmatch.fnmatch(vertex, pattern_list[0]):
                    vertex_prereqs = dep_graph.get(vertex, set())
                    prereq_patterns = pattern_list[1:]
                    if not any(
                        fnmatch.filter(vertex_prereqs, prereq_pattern)
                        for prereq_pattern in prereq_patterns
                    ):
                        warn(
                            "target '{0}' does not have an immediate "
                            "prerequisite matching any of the patterns: "
                            "'{1}'".format(
                                vertex, "', '".join(prereq_patterns)
                            ),
                            args.check_colour,
                        )

        start_visit_cb_list.append(check_exists_prereq_start_visit_cb)

    if args.check_cycles:
        path = []

        def check_cycles_start_visit_cb(vertex):
            path.append(vertex)

        def check_cycles_skip_visit_cb(vertex):
            if vertex in path:
                start_cycle_idx = path.index(vertex)

                if args.prereq is None:
                    msg_lines = (
                        path[1:start_cycle_idx]
                        + [path[start_cycle_idx] + " <- start of cycle"]
                        + path[start_cycle_idx + 1 :]
                        + [vertex + " <- end of cycle"]
                    )
                else:
                    msg_lines = (
                        [vertex + " <- start of cycle"]
                        + path[-1:start_cycle_idx:-1]
                        + [path[start_cycle_idx] + " <- end of cycle"]
                        + path[start_cycle_idx - 1 : 0 : -1]
                    )

                warn(
                    "the dependency graph has a cycle:\n"
                    "\t{0}".format("\n\t".join(msg_lines)),
                    args.check_colour,
                )

        def check_cycles_finish_visit_cb(_):
            path.pop()

        start_visit_cb_list.append(check_cycles_start_visit_cb)
        skip_visit_cb_list.append(check_cycles_skip_visit_cb)
        finish_visit_cb_list.append(check_cycles_finish_visit_cb)

    toposort = []

    def toposort_finish_visit_cb(vertex):
        toposort.append(vertex)

    def toposort_postprocess_cb():
        # The last element of toposort is _meta_root:
        toposort.pop()

    finish_visit_cb_list.append(toposort_finish_visit_cb)
    postprocess_cb_list.append(toposort_postprocess_cb)

    visited_vertices = dict()

    def traverse(start_depth=-1):
        visit_dfs(
            traversed_graph,
            _meta_root,
            current_depth=start_depth,
            max_depth=args.max_depth,
            visited=visited_vertices,
            start_visit_cb_list=start_visit_cb_list,
            finish_visit_cb_list=finish_visit_cb_list,
            skip_visit_cb_list=skip_visit_cb_list,
        )

        for postprocess_cb in postprocess_cb_list:
            postprocess_cb()

    traverse()

    # Add the extra prerequisites to the graph:
    for target, prereqs in extra_edges.items():
        traversed_graph[target] = traversed_graph.default_factory(
            dedupe(itertools.chain(traversed_graph[target], prereqs))
        )

    for target, prereqs in extra_edges.items():
        target_depth = visited_vertices.get(target, None)
        if target_depth is None:
            continue

        # Reset the _meta_root and traverse the graph:
        visited_vertices.pop(_meta_root)
        traversed_graph[_meta_root] = prereqs
        traverse(target_depth)

    if args.reverse ^ (args.prereq is not None):
        toposort.reverse()

    print("\n".join(toposort))


if __name__ == "__main__":
    main()
