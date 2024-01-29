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

# ACX_HOST_FQDN([ACTION-IF-SUCCESS],
#               [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Detects the fully qualified domain name (FQDN) of the host. The result is
# either "unknown" or a string with the fully qualified domain name.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_host_fqdn variable.
#
AC_DEFUN([ACX_HOST_FQDN],
  [AC_CACHE_CHECK([for the fully qualified domain name of the host],
     [acx_cv_host_fqdn],
     [dnl The easiest way to get the result:
      acx_cv_host_fqdn=`hostname -f 2>/dev/null`
      AS_IF([test $? -ne 0 || test -z "$acx_cv_host_fqdn"],
        [acx_cv_host_fqdn=unknown
dnl If the previous attempt didn't work, try to find the hostname first:
dnl   use the 'hostname' and 'uname -n' commands:
         acx_hostname=`(hostname || uname -n) 2>/dev/null | sed 1q`
dnl   or use the value of the HOSTNAME environment variable:
         test -n "$acx_hostname" || acx_hostname=$HOSTNAME
dnl   or read '/proc/sys/kernel/hostname' (not portable):
         test -n "$acx_hostname" || acx_hostname=`cat /proc/sys/kernel/hostname 2>/dev/null`
dnl   Try to call to the 'host' command if the hostname is detected:
         test -n "$acx_hostname" && dnl
acx_cv_host_fqdn=`host "$acx_hostname" 2>/dev/null | dnl
awk '/has address/ {print $[1]}' 2>/dev/null`])
dnl Set the result to unknown if it is empty:
      test -n "$acx_cv_host_fqdn" || acx_cv_host_fqdn=unknown])
   AS_VAR_IF([acx_cv_host_fqdn], [unknown], [m4_default([$2],
     [AC_MSG_FAILURE([unable to detect the fully qualified domain name of dnl
the host])])], [$1])])
