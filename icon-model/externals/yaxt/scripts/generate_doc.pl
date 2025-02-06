#! /usr/bin/env perl
#
# scripts/generate_doc.pl --- script for yaxt tests
#
# Copyright  (C)  2022 Jörg Behrens <behrens@dkrz.de>
#                      Moritz Hanke <hanke@dkrz.de>
#                      Thomas Jahns <jahns@dkrz.de>
#
# Author: Jörg Behrens <behrens@dkrz.de>
#         Moritz Hanke <hanke@dkrz.de>
#         Thomas Jahns <jahns@dkrz.de>
#
# Maintainer: Jörg Behrens <behrens@dkrz.de>
#             Moritz Hanke <hanke@dkrz.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: http://https://www.dkrz.de/redmine/projects/scales-ppm
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

use strict;
use warnings;

# find package version in AC_INIT macro of configure.ac
sub match_version_from_ac_init {
  my $file = shift || die "sub match_version_from_ac_init: missing argument";
  my $fh;
  open($fh, '<', $file) || die "Cannot open $file.";
  while(my $line = <$fh>) {
    if ($line =~ /^\s*AC_INIT/) {
      $line =~ s/\s+//g;
      if ($line =~ /^AC_INIT\(\[yaxt\],\[([\d\.\w]+)\]/) {
        return $1;
      }
    }
  }
  close($fh);
  return undef;
}

sub create_doxyfile {
  my ($version, $in, $out) = @_;
  die "sub create_doxyfile: missing argument"
    unless defined($version) and defined($in) and defined($out);
  my ($fh_in, $fh_out);
  open($fh_in, '<', $in) || die "Cannot open $in.";
  open($fh_out, '>', $out) || die "Cannot open $out.";
  my $match = 0;
  while(my $line = <$fh_in>) {
    if ($line =~ /^(\s*PROJECT_NUMBER\s*=\s*)/) {
      $line = "PROJECT_NUMBER = $version\n";
      $match++;
    }
    print $fh_out $line;
  }
  close($fh_out);
  close($fh_in);
  return $match;
}

my ($conf, $doxy_in, $doxy_out) = @ARGV;

die "$0 usage: your_configure.ac Doxygen_template outputfile\n"
  unless defined($conf) and defined($doxy_in) and defined($doxy_out);


my $version = match_version_from_ac_init($conf);
die "Cannot find yaxt version in $conf." unless defined $version;

create_doxyfile($version, $doxy_in, $doxy_out) || die "Found no match in $doxy_in\n";


