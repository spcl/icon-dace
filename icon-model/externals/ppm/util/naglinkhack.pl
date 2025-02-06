#! /usr/bin/env perl
#
# naglinkhack.pl --- translation of linker flags for NAG Fortran compiler
#
# Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
#
# Version: 1.0
# Keywords: NAG compiler -Wl fix
# Author: Thomas Jahns <jahns@dkrz.de>
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: https://www.dkrz.de/redmine/projects/scales-ppm
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
# The NAG Fortran compiler has its own definition of -Wl: -Wl
# arguments are passed to the underlying C compiler, to pass them to
# the link editor, they must be escaped by prepending another -Wl and
# doubling every comma.
use strict;
use warnings;

while (<>)
{
  my $output = $_;
  # generate -Wl,-Wl,,x,,y from -Wl,x,y
  my ($suffix, $prefix) = ($output, '');
  while (length($suffix))
  {
    my $prefix_add;
    ($prefix_add, $suffix) = $suffix =~ m{^(.*?)(-Wl,.*)?$};
    $prefix_add = '' if !defined($prefix_add);
    $suffix = '' if !defined($suffix);
    $prefix .= $prefix_add;
    if (length($suffix))
    {
      my ($arg);
      ($arg, $suffix) = $suffix =~ m{^(-Wl,\S*)(.*)};
      $arg =~ s/,/,,/g;
      $arg = '-Wl,' . $arg;
      $prefix .= $arg;
    }
  }
  $output = $prefix;
  # remove -pthread
  $output =~ s{(?:^|(?<=\s))-pthread(?=\s|$)}{}gx;
  print $output, "\n";
}

#
# Local Variables:
# license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:
#
