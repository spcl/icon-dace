#! /bin/bash
#
# recreate-testsuite.sh --- regenerate configure-dependent variables
#                           for autotest runs
#
#
# Copyright  (C)  2021  Thomas Jahns <jahns@dkrz.de>
#
# Version: 1.0
# Author: Thomas Jahns <jahns@dkrz.de>
# Keywords:
# Maintainer: Thomas Jahns <jahns@dkrz.de>
# URL: http://
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
set -e
if [[ ! -r configure || ! -d test || ! -r test/testsuite.at ]]; then
  echo "error: change to source directory first!" >&2
  exit 1
fi
cd test
coproc grep -E '^PACKAGE_(NAME|TARNAME|VERSION|STRING|BUGREPORT|URL)=' \
     ../configure
exec 5>package.m4.temp
trap '\rm package.m4.temp' ERR
echo '# Signature of the current package.' >&5
while IFS='=' read sym val
do
  val=${val%\'}
  val=${val#\'}
  echo 'm4_define([AT_'"$sym"'],' >&5
  echo '  ['"$val"'])' >&5
done <&${COPROC[0]}
mv package.m4.temp package.m4
exec autom4te --language=autotest -I . -o testsuite testsuite.at
#
# Local Variables:
# mode: sh
# license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
# license-default: "bsd"
# End:
#
