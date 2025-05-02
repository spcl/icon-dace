This patch is needed so that -Wl,-rpath option arguments make it into
the dependent libraries of a library and also ensures that -Wl,-Wl,,
options are passed to the compiler unaltered.

It is intended for libtool 2.4.2 and can be applied after
e.g. autoreconf by

$ patch -p1 <contrib/nagfor-libtool-patch/nagfor-libtool.patch

