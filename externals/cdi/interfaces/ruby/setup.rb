require 'mkmf'

$INCFLAGS = "$(CDI_INCFLAGS) #{$INCFLAGS}"
$CPPFLAGS = "$(CDI_CPPFLAGS) #{$CPPFLAGS}"
$CXXFLAGS = "$(CDI_CXXFLAGS) #{$CXXFLAGS}"
$LDFLAGS  = "$(CDI_LDFLAGS) #{$LDFLAGS}"
$LIBS     = "#{$LIBS} $(CDI_LIBS)"

$libs = "../libcdipp.la #{$libs}"
$srcs = %w[cdi_wrapper.cpp]
$objs = %w[cdi_wrapper.lo]
create_makefile('Cdi')
