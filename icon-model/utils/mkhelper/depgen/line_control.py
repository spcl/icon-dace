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
import re

from depgen import file_in_dir


class Parser:
    _re_lc = re.compile(r'^#\s*[1-9]\d*\s*"(.*?)"\s*(?:[1-9]\d*)?')

    def __init__(self, include_roots=None, subparser=None):
        self.include_roots = include_roots

        self._get_stream_iterator = (
            subparser.parse if subparser else lambda x, *_: x
        )

        # Callbacks:
        self.lc_callback = None
        self.debug_callback = None

    def parse(self, stream, stream_name):
        stream = self._get_stream_iterator(stream, stream_name)

        for line in stream:
            match = Parser._re_lc.match(line)
            if match:
                filepath = match.group(1)
                if os.path.isfile(filepath):
                    if not self.include_roots or any(
                        [file_in_dir(filepath, d) for d in self.include_roots]
                    ):
                        if self.lc_callback:
                            self.lc_callback(filepath)
                        if self.debug_callback:
                            self.debug_callback(
                                line, "accepted file '{0}'".format(filepath)
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

            yield line
