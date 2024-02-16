#
# @file Makefile.am
#
# @copyright Copyright  (C)  2020 Moritz Hanke <hanke@dkrz.de>
#                                 Rene Redler <rene.redler@mpimet.mpg.de>
#                                 Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#                                 Teresa Holfeld <teresa.holfeld@zmaw.de>
#                                 Thomas Jahns <jahns@dkrz.de>
#
# @version 1.0
# @author Moritz Hanke <hanke@dkrz.de>
#         Rene Redler <rene.redler@mpimet.mpg.de>
#         Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#         Teresa Holfeld <teresa.holfeld@zmaw.de>
#         Thomas Jahns <jahns@dkrz.de>
#
# Keywords:
# Maintainer: Moritz Hanke <hanke@dkrz.de>
#             Rene Redler <rene.redler@mpimet.mpg.de>
#             Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#             Teresa Holfeld <teresa.holfeld@zmaw.de>
#             Thomas Jahns <jahns@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
#
# This file is part of YAC.
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

lib_LIBRARIES = libyac_mtime.a

libyac_mtime_a_SOURCES =                  \
        include/mtime_calendar.h          \
        include/mtime_calendar360day.h    \
        include/mtime_calendar365day.h    \
        include/mtime_calendarGregorian.h \
        include/mtime_date.h              \
        include/mtime_datetime.h          \
        include/mtime_eventHandling.h     \
        include/mtime_eventList.h         \
        include/mtime_iso8601.h           \
        include/mtime_julianDay.h         \
        include/mtime_time.h              \
        include/mtime_timedelta.h         \
        src/mtime_calendar.c              \
        src/mtime_calendar360day.c        \
        src/mtime_calendar365day.c        \
        src/mtime_calendarGregorian.c     \
        src/mtime_date.c                  \
        src/mtime_datetime.c              \
        src/mtime_eventHandling.c         \
        src/mtime_eventList.c             \
        src/mtime_iso8601.c               \
        src/mtime_julianDay.c             \
        src/mtime_time.c                  \
        src/mtime_timedelta.c

AM_CPPFLAGS = -I$(top_srcdir)/mtime/include
AM_CFLAGS = $(PIC_CFLAGS)
