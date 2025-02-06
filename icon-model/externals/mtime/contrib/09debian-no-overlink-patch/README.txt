The patches are taken the Debian project:
https://sources.debian.org/patches/libtool/2.4.7-7/deplib_binary.patch
https://sources.debian.org/patches/libtool/2.4.7-7/link_all_deplibs.patch

The patches prevent overlinking.

The comment in link_all_deplibs.patch says:
## Do not link against deplibs.  This is not needed for shared libs
## on atleast ELF systems since those already know which libs they
## need themself.  This seems to break a few things and will be fixed
## in a better way in a future upstream version.
