# this file overrides the macro from fortran.m4 until the macro is
# fixed in upstream autoconf

# _AC_PROG_FC_V_OUTPUT([FLAG = $ac_cv_prog_{f77/fc}_v])
# -------------------------------------------------
# Link a trivial Fortran program, compiling with a verbose output FLAG
# (whose default value, $ac_cv_prog_{f77/fc}_v, is computed by
# _AC_PROG_FC_V), and return the output in $ac_{f77/fc}_v_output.  This
# output is processed in the way expected by _AC_FC_LIBRARY_LDFLAGS,
# so that any link flags that are echoed by the compiler appear as
# space-separated items.
AC_DEFUN([_AC_PROG_FC_V_OUTPUT],
[_AC_FORTRAN_ASSERT()dnl
AC_LANG_CONFTEST([AC_LANG_PROGRAM([])])

# Compile and link our simple test program by passing a flag (argument
# 1 to this macro) to the Fortran compiler in order to get
# "verbose" output that we can then parse for the Fortran linker
# flags.
ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
_AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS m4_default([$1], [$ac_cv_prog_[]_AC_LANG_ABBREV[]_v])"
eval "set x $ac_link"
shift
_AS_ECHO_LOG([$[*]])
# gfortran 4.3 outputs lines setting COLLECT_GCC_OPTIONS, COMPILER_PATH,
# LIBRARY_PATH; skip all such settings.
ac_[]_AC_LANG_ABBREV[]_v_output=`eval $ac_link AS_MESSAGE_LOG_FD>&1 2>&1 |
  grep -v 'Driving:' | grep -v "^[[_$as_cr_Letters]][[_$as_cr_alnum]]*="`
AS_ECHO(["$ac_[]_AC_LANG_ABBREV[]_v_output"]) >&AS_MESSAGE_LOG_FD
_AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS

rm -rf conftest*

# On HP/UX there is a line like: "LPATH is: /foo:/bar:/baz" where
# /foo, /bar, and /baz are search directories for the Fortran linker.
# Here, we change these into -L/foo -L/bar -L/baz (and put it first):
ac_[]_AC_LANG_ABBREV[]_v_output="`echo $ac_[]_AC_LANG_ABBREV[]_v_output |
	grep 'LPATH is:' |
	sed 's|.*LPATH is\(: *[[^ ]]*\).*|\1|;s|: */| -L/|g'` $ac_[]_AC_LANG_ABBREV[]_v_output"

# FIXME: we keep getting bitten by quoted arguments; a more general fix
#        that detects unbalanced quotes in FLIBS should be implemented
#        and (ugh) tested at some point.
case $ac_[]_AC_LANG_ABBREV[]_v_output in
  # If we are using xlf then replace all the commas with spaces.
  # Also xlf appends a closing paren on its output and the ipa step adds a
  # -link argument that might be mistaken for a library
  *xlfentry*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo "$ac_[]_AC_LANG_ABBREV[]_v_output" | sed 's/,/ /g
s/) *$//g
s/ -link / /g'` ;;
  # With Intel ifc, ignore the quoted -mGLOB_options_string stuff (quoted
  # $LIBS confuse us, and the libraries appear later in the output anyway).
  *mGLOB_options_string*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/"-mGLOB[[^"]]*"/ /g'` ;;

  # Portland Group compiler has singly- or doubly-quoted -cmdline argument
  # Singly-quoted arguments were reported for versions 5.2-4 and 6.0-4.
  # Doubly-quoted arguments were reported for "PGF90/x86 Linux/x86 5.0-2".
  *-cmdline\ * | *-ignore\ * | *-def\ *)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed "\
	s/-cmdline  *'[[^']]*'/ /g; s/-cmdline  *\"[[^\"]]*\"/ /g
	s/-ignore  *'[[^']]*'/ /g; s/-ignore  *\"[[^\"]]*\"/ /g
	s/-def  *'[[^']]*'/ /g; s/-def  *\"[[^\"]]*\"/ /g"` ;;

  # If we are using fort77 (the f2c wrapper) then filter output and delete quotes.
  *fort77*f2c*gcc*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo "$ac_[]_AC_LANG_ABBREV[]_v_output" | sed -n '
        /:[[	 ]]\+Running[[	 ]]\{1,\}"gcc"/{
          /"-c"/d
          /[[.]]c"*/d
          s/^.*"gcc"/"gcc"/
          s/"//gp
        }'` ;;

  # If we are using Cray Fortran then delete quotes.
  *cft90*)
    ac_[]_AC_LANG_ABBREV[]_v_output=`echo $ac_[]_AC_LANG_ABBREV[]_v_output | sed 's/"//g'` ;;
esac

])# _AC_PROG_FC_V_OUTPUT
dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
