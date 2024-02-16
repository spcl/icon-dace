lib_LIBRARIES = libyac_clapack.a

BUILT_SOURCES = sysdep1.h
CLEANFILES = sysdep1.h
EXTRA_DIST = F2CLIBS/libf2c/sysdep1.h0

sysdep1.h: F2CLIBS/libf2c/sysdep1.h0
	$(AM_V_GEN)cp $< $@

AM_CPPFLAGS = -I$(srcdir)/INCLUDE
AM_CFLAGS = $(PIC_CFLAGS)

libyac_clapack_a_SOURCES =        \
        BLAS/SRC/dcopy.c          \
        BLAS/SRC/ddot.c           \
        BLAS/SRC/dgemm.c          \
        BLAS/SRC/dgemv.c          \
        BLAS/SRC/dger.c           \
        BLAS/SRC/dnrm2.c          \
        BLAS/SRC/dscal.c          \
        BLAS/SRC/dswap.c          \
        BLAS/SRC/dsymv.c          \
        BLAS/SRC/dsyr.c           \
        BLAS/SRC/dtrmm.c          \
        BLAS/SRC/dtrmv.c          \
        BLAS/SRC/dtrsm.c          \
        BLAS/SRC/idamax.c         \
        F2CLIBS/libf2c/close.c    \
        F2CLIBS/libf2c/d_lg10.c   \
        F2CLIBS/libf2c/d_sign.c   \
        F2CLIBS/libf2c/endfile.c  \
        F2CLIBS/libf2c/err.c      \
        F2CLIBS/libf2c/exit_.c    \
        F2CLIBS/libf2c/f77_aloc.c \
        F2CLIBS/libf2c/fio.h      \
        F2CLIBS/libf2c/fmt.c      \
        F2CLIBS/libf2c/fmt.h      \
        F2CLIBS/libf2c/fmtlib.c   \
        F2CLIBS/libf2c/fp.h       \
        F2CLIBS/libf2c/i_nint.c   \
        F2CLIBS/libf2c/open.c     \
        F2CLIBS/libf2c/pow_di.c   \
        F2CLIBS/libf2c/s_cat.c    \
        F2CLIBS/libf2c/s_cmp.c    \
        F2CLIBS/libf2c/s_copy.c   \
        F2CLIBS/libf2c/sfe.c      \
        F2CLIBS/libf2c/sig_die.c  \
        F2CLIBS/libf2c/util.c     \
        F2CLIBS/libf2c/wref.c     \
        F2CLIBS/libf2c/wrtfmt.c   \
        F2CLIBS/libf2c/wsfe.c     \
        INCLUDE/blaswrap.h        \
        INCLUDE/clapack.h         \
        INCLUDE/f2c.h             \
        INSTALL/dlamch.c          \
        INSTALL/lsame.c           \
        SRC/dgelq2.c              \
        SRC/dgelqf.c              \
        SRC/dgels.c               \
        SRC/dgeqr2.c              \
        SRC/dgeqrf.c              \
        SRC/dgesv.c               \
        SRC/dgetf2.c              \
        SRC/dgetrf.c              \
        SRC/dgetri.c              \
        SRC/dgetrs.c              \
        SRC/disnan.c              \
        SRC/dlabad.c              \
        SRC/dlaisnan.c            \
        SRC/dlange.c              \
        SRC/dlapy2.c              \
        SRC/dlarfb.c              \
        SRC/dlarf.c               \
        SRC/dlarfp.c              \
        SRC/dlarft.c              \
        SRC/dlascl.c              \
        SRC/dlaset.c              \
        SRC/dlassq.c              \
        SRC/dlaswp.c              \
        SRC/dlasyf.c              \
        SRC/dorm2r.c              \
        SRC/dorml2.c              \
        SRC/dormlq.c              \
        SRC/dormqr.c              \
        SRC/dsytf2.c              \
        SRC/dsytrf.c              \
        SRC/dsytri.c              \
        SRC/dtrti2.c              \
        SRC/dtrtri.c              \
        SRC/dtrtrs.c              \
        SRC/ieeeck.c              \
        SRC/iladlc.c              \
        SRC/iladlr.c              \
        SRC/ilaenv.c              \
        SRC/iparmq.c              \
        SRC/xerbla.c
