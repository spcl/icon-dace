!
! mo_fplume_dlsode
! This module contains Livermore Solver (used in FPLUME)
! Original code by Hindmarsh, Alan C., (LLNL)
!
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_fplume_dlsode

  USE mo_exception,                     ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::dlsode
  CONTAINS
!DECK DLSODE
    SUBROUTINE dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
      iopt, rwork, lrw, iwork, liw, jac, mf)
      EXTERNAL f, jac
      INTEGER neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      DOUBLE PRECISION y, t, tout, rtol, atol, rwork
      DIMENSION neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)
!***BEGIN PROLOGUE  DLSODE
!***PURPOSE  Livermore Solver for Ordinary Differential Equations.
!            DLSODE solves the initial-value problem for stiff or
!            nonstiff systems of first-order ODE's,
!               dy/dt = f(t,y),   or, in component form,
!               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
!***CATEGORY  I1A
!***TYPE      DOUBLE PRECISION (SLSODE-S, DLSODE-D)
!***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
!             STIFF, NONSTIFF
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!             Center for Applied Scientific Computing, L-561
!             Lawrence Livermore National Laboratory
!             Livermore, CA 94551.
!***DESCRIPTION
!
!     NOTE: The "Usage" and "Arguments" sections treat only a subset of
!           available options, in condensed fashion.  The options
!           covered and the information supplied will support most
!           standard uses of DLSODE.
!
!           For more sophisticated uses, full details on all options are
!           given in the concluding section, headed "Long Description."
!           A synopsis of the DLSODE Long Description is provided at the
!           beginning of that section; general topics covered are:
!           - Elements of the call sequence; optional input and output
!           - Optional supplemental routines in the DLSODE package
!           - internal COMMON block
!
! *Usage:
!     Communication between the user and the DLSODE package, for normal
!     situations, is summarized here.  This summary describes a subset
!     of the available options.  See "Long Description" for complete
!     details, including optional communication, nonstandard options,
!     and instructions for special situations.
!
!     A sample program is given in the "Examples" section.
!
!     Refer to the argument descriptions for the definitions of the
!     quantities that appear in the following sample declarations.
!
!     For MF = 10,
!        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20)
!     For MF = 21 or 22,
!        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ)
!     For MF = 24 or 25,
!        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ,
!       *                                         LIW = 20 + NEQ)
!
!        EXTERNAL F, JAC
!        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW),
!       *         LIW, MF
!        DOUBLE PRECISION Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW)
!
!        CALL DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
!       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
!
! *Arguments:
!     F     :EXT    Name of subroutine for right-hand-side vector f.
!                   This name must be declared EXTERNAL in calling
!                   program.  The form of F must be:
!
!                   SUBROUTINE  F (NEQ, T, Y, YDOT)
!                   INTEGER  NEQ
!                   DOUBLE PRECISION  T, Y(*), YDOT(*)
!
!                   The inputs are NEQ, T, Y.  F is to set
!
!                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)),
!                                                     i = 1, ..., NEQ .
!
!     NEQ   :IN     Number of first-order ODE's.
!
!     Y     :INOUT  Array of values of the y(t) vector, of length NEQ.
!                   Input:  For the first call, Y should contain the
!                           values of y(t) at t = T. (Y is an input
!                           variable only if ISTATE = 1.)
!                   Output: On return, Y will contain the values at the
!                           new t-value.
!
!     T     :INOUT  Value of the independent variable.  On return it
!                   will be the current value of t (normally TOUT).
!
!     TOUT  :IN     Next point where output is desired (.NE. T).
!
!     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or
!                   an array.
!
!     RTOL  :IN     Relative tolerance parameter (scalar).
!
!     ATOL  :IN     Absolute tolerance parameter (scalar or array).
!                   If ITOL = 1, ATOL need not be dimensioned.
!                   If ITOL = 2, ATOL must be dimensioned at least NEQ.
!
!                   The estimated local error in Y(i) will be controlled
!                   so as to be roughly less (in magnitude) than
!
!                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!
!                   Thus the local error test passes if, in each
!                   component, either the absolute error is less than
!                   ATOL (or ATOL(i)), or the relative error is less
!                   than RTOL.
!
!                   Use RTOL = 0.0 for pure absolute error control, and
!                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative
!                   error control.  Caution:  Actual (global) errors may
!                   exceed these local tolerances, so choose them
!                   conservatively.
!
!     ITASK :IN     Flag indicating the task DLSODE is to perform.
!                   Use ITASK = 1 for normal computation of output
!                   values of y at t = TOUT.
!
!     ISTATE:INOUT  Index used for input and output to specify the state
!                   of the calculation.
!                   Input:
!                    1   This is the first call for a problem.
!                    2   This is a subsequent call.
!                   Output:
!                    1   Nothing was done, because TOUT was equal to T.
!                    2   DLSODE was successful (otherwise, negative).
!                        Note that ISTATE need not be modified after a
!                        successful return.
!                   -1   Excess work done on this call (perhaps wrong
!                        MF).
!                   -2   Excess accuracy requested (tolerances too
!                        small).
!                   -3   Illegal input detected (see printed message).
!                   -4   Repeated error test failures (check all
!                        inputs).
!                   -5   Repeated convergence failures (perhaps bad
!                        Jacobian supplied or wrong choice of MF or
!                        tolerances).
!                   -6   Error weight became zero during problem
!                        (solution component i vanished, and ATOL or
!                        ATOL(i) = 0.).
!
!     IOPT  :IN     Flag indicating whether optional inputs are used:
!                   0   No.
!                   1   Yes.  (See "Optional inputs" under "Long
!                       Description," Part 1.)
!
!     RWORK :WORK   Real work array of length at least:
!                   20 + 16*NEQ                    for MF = 10,
!                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
!                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
!
!     LRW   :IN     Declared length of RWORK (in user's DIMENSION
!                   statement).
!
!     IWORK :WORK   Integer work array of length at least:
!                   20        for MF = 10,
!                   20 + NEQ  for MF = 21, 22, 24, or 25.
!
!                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the
!                   lower and upper Jacobian half-bandwidths ML,MU.
!
!                   On return, IWORK contains information that may be
!                   of interest to the user:
!
!            Name   Location   Meaning
!            -----  ---------  -----------------------------------------
!            NST    IWORK(11)  Number of steps taken for the problem so
!                              far.
!            NFE    IWORK(12)  Number of f evaluations for the problem
!                              so far.
!            NJE    IWORK(13)  Number of Jacobian evaluations (and of
!                              matrix LU decompositions) for the problem
!                              so far.
!            NQU    IWORK(14)  Method order last used (successfully).
!            LENRW  IWORK(17)  Length of RWORK actually required.  This
!                              is defined on normal returns and on an
!                              illegal input return for insufficient
!                              storage.
!            LENIW  IWORK(18)  Length of IWORK actually required.  This
!                              is defined on normal returns and on an
!                              illegal input return for insufficient
!                              storage.
!
!     LIW   :IN     Declared length of IWORK (in user's DIMENSION
!                   statement).
!
!     JAC   :EXT    Name of subroutine for Jacobian matrix (MF =
!                   21 or 24).  If used, this name must be declared
!                   EXTERNAL in calling program.  If not used, pass a
!                   dummy name.  The form of JAC must be:
!
!                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!                   INTEGER  NEQ, ML, MU, NROWPD
!                   DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
!
!                   See item c, under "Description" below for more
!                   information about JAC.
!
!     MF    :IN     Method flag.  Standard values are:
!                   10  Nonstiff (Adams) method, no Jacobian used.
!                   21  Stiff (BDF) method, user-supplied full Jacobian.
!                   22  Stiff method, internally generated full
!                       Jacobian.
!                   24  Stiff method, user-supplied banded Jacobian.
!                   25  Stiff method, internally generated banded
!                       Jacobian.
!
! *Description:
!     DLSODE solves the initial value problem for stiff or nonstiff
!     systems of first-order ODE's,
!
!        dy/dt = f(t,y) ,
!
!     or, in component form,
!
!        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ))
!                                                  (i = 1, ..., NEQ) .
!
!     DLSODE is a package based on the GEAR and GEARB packages, and on
!     the October 23, 1978, version of the tentative ODEPACK user
!     interface standard, with minor modifications.
!
!     The steps in solving such a problem are as follows.
!
!     a. First write a subroutine of the form
!
!           SUBROUTINE  F (NEQ, T, Y, YDOT)
!           INTEGER  NEQ
!           DOUBLE PRECISION  T, Y(*), YDOT(*)
!
!        which supplies the vector function f by loading YDOT(i) with
!        f(i).
!
!     b. Next determine (or guess) whether or not the problem is stiff.
!        Stiffness occurs when the Jacobian matrix df/dy has an
!        eigenvalue whose real part is negative and large in magnitude
!        compared to the reciprocal of the t span of interest.  If the
!        problem is nonstiff, use method flag MF = 10.  If it is stiff,
!        there are four standard choices for MF, and DLSODE requires the
!        Jacobian matrix in some form.  This matrix is regarded either
!        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the
!        banded case, DLSODE requires two half-bandwidth parameters ML
!        and MU. These are, respectively, the widths of the lower and
!        upper parts of the band, excluding the main diagonal.  Thus the
!        band consists of the locations (i,j) with
!
!           i - ML <= j <= i + MU ,
!
!        and the full bandwidth is ML + MU + 1 .
!
!     c. If the problem is stiff, you are encouraged to supply the
!        Jacobian directly (MF = 21 or 24), but if this is not feasible,
!        DLSODE will compute it internally by difference quotients (MF =
!        22 or 25).  If you are supplying the Jacobian, write a
!        subroutine of the form
!
!           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!           INTEGER  NEQ, ML, MU, NRWOPD
!           DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
!
!        which provides df/dy by loading PD as follows:
!        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
!          the partial derivative of f(i) with respect to y(j).  (Ignore
!          the ML and MU arguments in this case.)
!        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
!          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the
!          rows of PD from the top down.
!        - In either case, only nonzero elements need be loaded.
!
!     d. Write a main program that calls subroutine DLSODE once for each
!        point at which answers are desired.  This should also provide
!        for possible use of logical unit 6 for output of error messages
!        by DLSODE.
!
!        Before the first call to DLSODE, set ISTATE = 1, set Y and T to
!        the initial values, and set TOUT to the first output point.  To
!        continue the integration after a successful return, simply
!        reset TOUT and call DLSODE again.  No other parameters need be
!        reset.
!
! *Examples:
!     The following is a simple example problem, with the coding needed
!     for its solution by DLSODE. The problem is from chemical kinetics,
!     and consists of the following three rate equations:
!
!        dy1/dt = -.04*y1 + 1.E4*y2*y3
!        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
!        dy3/dt = 3.E7*y2**2
!
!     on the interval from t = 0.0 to t = 4.E10, with initial conditions
!     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
!
!     The following coding solves this problem with DLSODE, using
!     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses
!     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2
!     has much smaller values.  At the end of the run, statistical
!     quantities of interest are printed.
!
!        EXTERNAL  FEX, JEX
!        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
!       *         MF, NEQ
!        DOUBLE PRECISION  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
!        NEQ = 3
!        Y(1) = 1.D0
!        Y(2) = 0.D0
!        Y(3) = 0.D0
!        T = 0.D0
!        TOUT = .4D0
!        ITOL = 2
!        RTOL = 1.D-4
!        ATOL(1) = 1.D-6
!        ATOL(2) = 1.D-10
!        ATOL(3) = 1.D-6
!        ITASK = 1
!        ISTATE = 1
!        IOPT = 0
!        LRW = 58
!        LIW = 23
!        MF = 21
!        DO 40 IOUT = 1,12
!          CALL DLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
!       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
!          WRITE(6,20)  T, Y(1), Y(2), Y(3)
!    20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
!          IF (ISTATE .LT. 0)  GO TO 80
!    40    TOUT = TOUT*10.D0
!        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
!    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
!        STOP
!    80  WRITE(6,90)  ISTATE
!    90  FORMAT(///' Error halt.. ISTATE =',I3)
!        STOP
!        END
!
!        SUBROUTINE  FEX (NEQ, T, Y, YDOT)
!        INTEGER  NEQ
!        DOUBLE PRECISION  T, Y(3), YDOT(3)
!        YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
!        YDOT(3) = 3.D7*Y(2)*Y(2)
!        YDOT(2) = -YDOT(1) - YDOT(3)
!        RETURN
!        END
!
!        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
!        INTEGER  NEQ, ML, MU, NRPD
!        DOUBLE PRECISION  T, Y(3), PD(NRPD,3)
!        PD(1,1) = -.04D0
!        PD(1,2) = 1.D4*Y(3)
!        PD(1,3) = 1.D4*Y(2)
!        PD(2,1) = .04D0
!        PD(2,3) = -PD(1,3)
!        PD(3,2) = 6.D7*Y(2)
!        PD(2,2) = -PD(1,2) - PD(3,2)
!        RETURN
!        END
!
!     The output from this program (on a Cray-1 in single precision)
!     is as follows.
!
!     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
!     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
!     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
!     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
!     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
!     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
!     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
!     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
!     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
!     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01
!     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
!     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00
!
!     No. steps = 330,  No. f-s = 405,  No. J-s = 69
!
! *Accuracy:
!     The accuracy of the solution depends on the choice of tolerances
!     RTOL and ATOL.  Actual (global) errors may exceed these local
!     tolerances, so choose them conservatively.
!
! *Cautions:
!     The work arrays should not be altered between calls to DLSODE for
!     the same problem, except possibly for the conditional and optional
!     inputs.
!
! *Portability:
!     Since NEQ is dimensioned inside DLSODE, some compilers may object
!     to a call to DLSODE with NEQ a scalar variable.  In this event,
!     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL.
!
!     Note to Cray users:
!     For maximum efficiency, use the CFT77 compiler.  Appropriate
!     compiler optimization directives have been inserted for CFT77.
!
! *Reference:
!     Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
!     Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
!     (North-Holland, Amsterdam, 1983), pp. 55-64.
!
! *Long Description:
!     The following complete description of the user interface to
!     DLSODE consists of four parts:
!
!     1.  The call sequence to subroutine DLSODE, which is a driver
!         routine for the solver.  This includes descriptions of both
!         the call sequence arguments and user-supplied routines.
!         Following these descriptions is a description of optional
!         inputs available through the call sequence, and then a
!         description of optional outputs in the work arrays.
!
!     2.  Descriptions of other routines in the DLSODE package that may
!         be (optionally) called by the user.  These provide the ability
!         to alter error message handling, save and restore the internal
!         COMMON, and obtain specified derivatives of the solution y(t).
!
!     3.  Descriptions of COMMON block to be declared in overlay or
!         similar environments, or to be saved when doing an interrupt
!         of the problem and continued solution later.
!
!     4.  Description of two routines in the DLSODE package, either of
!         which the user may replace with his own version, if desired.
!         These relate to the measurement of errors.
!
!
!                         Part 1.  Call Sequence
!                         ----------------------
!
!     Arguments
!     ---------
!     The call sequence parameters used for input only are
!
!        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
!
!     and those used for both input and output are
!
!        Y, T, ISTATE.
!
!     The work arrays RWORK and IWORK are also used for conditional and
!     optional inputs and optional outputs.  (The term output here
!     refers to the return from subroutine DLSODE to the user's calling
!     program.)
!
!     The legality of input parameters will be thoroughly checked on the
!     initial call for the problem, but not checked thereafter unless a
!     change in input parameters is flagged by ISTATE = 3 on input.
!
!     The descriptions of the call arguments are as follows.
!
!     F        The name of the user-supplied subroutine defining the ODE
!              system.  The system must be put in the first-order form
!              dy/dt = f(t,y), where f is a vector-valued function of
!              the scalar t and the vector y. Subroutine F is to compute
!              the function f. It is to have the form
!
!                 SUBROUTINE F (NEQ, T, Y, YDOT)
!                 DOUBLE PRECISION  T, Y(*), YDOT(*)
!
!              where NEQ, T, and Y are input, and the array YDOT =
!              f(T,Y) is output.  Y and YDOT are arrays of length NEQ.
!              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be
!              declared EXTERNAL in the calling program.
!
!              Subroutine F may access user-defined quantities in
!              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array
!              (dimensioned in F) and/or Y has length exceeding NEQ(1).
!              See the descriptions of NEQ and Y below.
!
!              If quantities computed in the F routine are needed
!              externally to DLSODE, an extra call to F should be made
!              for this purpose, for consistent and accurate results.
!              If only the derivative dy/dt is needed, use DINTDY
!              instead.
!
!     NEQ      The size of the ODE system (number of first-order
!              ordinary differential equations).  Used only for input.
!              NEQ may be decreased, but not increased, during the
!              problem.  If NEQ is decreased (with ISTATE = 3 on input),
!              the remaining components of Y should be left undisturbed,
!              if these are to be accessed in F and/or JAC.
!
!              Normally, NEQ is a scalar, and it is generally referred
!              to as a scalar in this user interface description.
!              However, NEQ may be an array, with NEQ(1) set to the
!              system size.  (The DLSODE package accesses only NEQ(1).)
!              In either case, this parameter is passed as the NEQ
!              argument in all calls to F and JAC.  Hence, if it is an
!              array, locations NEQ(2),... may be used to store other
!              integer data and pass it to F and/or JAC.  Subroutines
!              F and/or JAC must include NEQ in a DIMENSION statement
!              in that case.
!
!     Y        A real array for the vector of dependent variables, of
!              length NEQ or more.  Used for both input and output on
!              the first call (ISTATE = 1), and only for output on
!              other calls.  On the first call, Y must contain the
!              vector of initial values.  On output, Y contains the
!              computed solution vector, evaluated at T. If desired,
!              the Y array may be used for other purposes between
!              calls to the solver.
!
!              This array is passed as the Y argument in all calls to F
!              and JAC.  Hence its length may exceed NEQ, and locations
!              Y(NEQ+1),... may be used to store other real data and
!              pass it to F and/or JAC.  (The DLSODE package accesses
!              only Y(1),...,Y(NEQ).)
!
!     T        The independent variable.  On input, T is used only on
!              the first call, as the initial point of the integration.
!              On output, after each call, T is the value at which a
!              computed solution Y is evaluated (usually the same as
!              TOUT).  On an error return, T is the farthest point
!              reached.
!
!     TOUT     The next value of T at which a computed solution is
!              desired.  Used only for input.
!
!              When starting the problem (ISTATE = 1), TOUT may be equal
!              to T for one call, then should not equal T for the next
!              call.  For the initial T, an input value of TOUT .NE. T
!              is used in order to determine the direction of the
!              integration (i.e., the algebraic sign of the step sizes)
!              and the rough scale of the problem.  Integration in
!              either direction (forward or backward in T) is permitted.
!
!              If ITASK = 2 or 5 (one-step modes), TOUT is ignored
!              after the first call (i.e., the first call with
!              TOUT .NE. T).  Otherwise, TOUT is required on every call.
!
!              If ITASK = 1, 3, or 4, the values of TOUT need not be
!              monotone, but a value of TOUT which backs up is limited
!              to the current internal T interval, whose endpoints are
!              TCUR - HU and TCUR.  (See "Optional Outputs" below for
!              TCUR and HU.)
!
!
!     ITOL     An indicator for the type of error control.  See
!              description below under ATOL.  Used only for input.
!
!     RTOL     A relative error tolerance parameter, either a scalar or
!              an array of length NEQ.  See description below under
!              ATOL.  Input only.
!
!     ATOL     An absolute error tolerance parameter, either a scalar or
!              an array of length NEQ.  Input only.
!
!              The input parameters ITOL, RTOL, and ATOL determine the
!              error control performed by the solver.  The solver will
!              control the vector e = (e(i)) of estimated local errors
!              in Y, according to an inequality of the form
!
!                 rms-norm of ( e(i)/EWT(i) ) <= 1,
!
!              where
!
!                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!
!              and the rms-norm (root-mean-square norm) here is
!
!                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ).
!
!              Here EWT = (EWT(i)) is a vector of weights which must
!              always be positive, and the values of RTOL and ATOL
!              should all be nonnegative.  The following table gives the
!              types (scalar/array) of RTOL and ATOL, and the
!              corresponding form of EWT(i).
!
!              ITOL    RTOL      ATOL      EWT(i)
!              ----    ------    ------    -----------------------------
!              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
!              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
!              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
!              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!              When either of these parameters is a scalar, it need not
!              be dimensioned in the user's calling program.
!
!              If none of the above choices (with ITOL, RTOL, and ATOL
!              fixed throughout the problem) is suitable, more general
!              error controls can be obtained by substituting
!              user-supplied routines for the setting of EWT and/or for
!              the norm calculation.  See Part 4 below.
!
!              If global errors are to be estimated by making a repeated
!              run on the same problem with smaller tolerances, then all
!              components of RTOL and ATOL (i.e., of EWT) should be
!              scaled down uniformly.
!
!     ITASK    An index specifying the task to be performed.  Input
!              only.  ITASK has the following values and meanings:
!              1   Normal computation of output values of y(t) at
!                  t = TOUT (by overshooting and interpolating).
!              2   Take one step only and return.
!              3   Stop at the first internal mesh point at or beyond
!                  t = TOUT and return.
!              4   Normal computation of output values of y(t) at
!                  t = TOUT but without overshooting t = TCRIT.  TCRIT
!                  must be input as RWORK(1).  TCRIT may be equal to or
!                  beyond TOUT, but not behind it in the direction of
!                  integration.  This option is useful if the problem
!                  has a singularity at or beyond t = TCRIT.
!              5   Take one step, without passing TCRIT, and return.
!                  TCRIT must be input as RWORK(1).
!
!              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!              (within roundoff), it will return T = TCRIT (exactly) to
!              indicate this (unless ITASK = 4 and TOUT comes before
!              TCRIT, in which case answers at T = TOUT are returned
!              first).
!
!     ISTATE   An index used for input and output to specify the state
!              of the calculation.
!
!              On input, the values of ISTATE are as follows:
!              1   This is the first call for the problem
!                  (initializations will be done).  See "Note" below.
!              2   This is not the first call, and the calculation is to
!                  continue normally, with no change in any input
!                  parameters except possibly TOUT and ITASK.  (If ITOL,
!                  RTOL, and/or ATOL are changed between calls with
!                  ISTATE = 2, the new values will be used but not
!                  tested for legality.)
!              3   This is not the first call, and the calculation is to
!                  continue normally, but with a change in input
!                  parameters other than TOUT and ITASK.  Changes are
!                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
!                  ML, MU, and any of the optional inputs except H0.
!                  (See IWORK description for ML and MU.)
!
!              Note:  A preliminary call with TOUT = T is not counted as
!              a first call here, as no initialization or checking of
!              input is done.  (Such a call is sometimes useful for the
!              purpose of outputting the initial conditions.)  Thus the
!              first call for which TOUT .NE. T requires ISTATE = 1 on
!              input.
!
!              On output, ISTATE has the following values and meanings:
!               1  Nothing was done, as TOUT was equal to T with
!                  ISTATE = 1 on input.
!               2  The integration was performed successfully.
!              -1  An excessive amount of work (more than MXSTEP steps)
!                  was done on this call, before completing the
!                  requested task, but the integration was otherwise
!                  successful as far as T. (MXSTEP is an optional input
!                  and is normally 500.)  To continue, the user may
!                  simply reset ISTATE to a value >1 and call again (the
!                  excess work step counter will be reset to 0).  In
!                  addition, the user may increase MXSTEP to avoid this
!                  error return; see "Optional Inputs" below.
!              -2  Too much accuracy was requested for the precision of
!                  the machine being used.  This was detected before
!                  completing the requested task, but the integration
!                  was successful as far as T. To continue, the
!                  tolerance parameters must be reset, and ISTATE must
!                  be set to 3. The optional output TOLSF may be used
!                  for this purpose.  (Note:  If this condition is
!                  detected before taking any steps, then an illegal
!                  input return (ISTATE = -3) occurs instead.)
!              -3  Illegal input was detected, before taking any
!                  integration steps.  See written message for details.
!                  (Note:  If the solver detects an infinite loop of
!                  calls to the solver with illegal input, it will cause
!                  the run to stop.)
!              -4  There were repeated error-test failures on one
!                  attempted step, before completing the requested task,
!                  but the integration was successful as far as T.  The
!                  problem may have a singularity, or the input may be
!                  inappropriate.
!              -5  There were repeated convergence-test failures on one
!                  attempted step, before completing the requested task,
!                  but the integration was successful as far as T. This
!                  may be caused by an inaccurate Jacobian matrix, if
!                  one is being used.
!              -6  EWT(i) became zero for some i during the integration.
!                  Pure relative error control (ATOL(i)=0.0) was
!                  requested on a variable which has now vanished.  The
!                  integration was successful as far as T.
!
!              Note:  Since the normal output value of ISTATE is 2, it
!              does not need to be reset for normal continuation.  Also,
!              since a negative input value of ISTATE will be regarded
!              as illegal, a negative output value requires the user to
!              change it, and possibly other inputs, before calling the
!              solver again.
!
!     IOPT     An integer flag to specify whether any optional inputs
!              are being used on this call.  Input only.  The optional
!              inputs are listed under a separate heading below.
!              0   No optional inputs are being used.  Default values
!                  will be used in all cases.
!              1   One or more optional inputs are being used.
!
!     RWORK    A real working array (double precision).  The length of
!              RWORK must be at least
!
!                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
!
!              where
!                 NYH = the initial value of NEQ,
!              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                       smaller value is given as an optional input),
!                 LWM = 0           if MITER = 0,
!                 LWM = NEQ**2 + 2  if MITER = 1 or 2,
!                 LWM = NEQ + 2     if MITER = 3, and
!                 LWM = (2*ML + MU + 1)*NEQ + 2
!                                   if MITER = 4 or 5.
!              (See the MF description below for METH and MITER.)
!
!              Thus if MAXORD has its default value and NEQ is constant,
!              this length is:
!              20 + 16*NEQ                    for MF = 10,
!              22 + 16*NEQ + NEQ**2           for MF = 11 or 12,
!              22 + 17*NEQ                    for MF = 13,
!              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15,
!              20 +  9*NEQ                    for MF = 20,
!              22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
!              22 + 10*NEQ                    for MF = 23,
!              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
!
!              The first 20 words of RWORK are reserved for conditional
!              and optional inputs and optional outputs.
!
!              The following word in RWORK is a conditional input:
!              RWORK(1) = TCRIT, the critical value of t which the
!                         solver is not to overshoot.  Required if ITASK
!                         is 4 or 5, and ignored otherwise.  See ITASK.
!
!     LRW      The length of the array RWORK, as declared by the user.
!              (This will be checked by the solver.)
!
!     IWORK    An integer work array.  Its length must be at least
!              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25).
!              (See the MF description below for MITER.)  The first few
!              words of IWORK are used for conditional and optional
!              inputs and optional outputs.
!
!              The following two words in IWORK are conditional inputs:
!              IWORK(1) = ML   These are the lower and upper half-
!              IWORK(2) = MU   bandwidths, respectively, of the banded
!                              Jacobian, excluding the main diagonal.
!                         The band is defined by the matrix locations
!                         (i,j) with i - ML <= j <= i + MU. ML and MU
!                         must satisfy 0 <= ML,MU <= NEQ - 1. These are
!                         required if MITER is 4 or 5, and ignored
!                         otherwise.  ML and MU may in fact be the band
!                         parameters for a matrix to which df/dy is only
!                         approximately equal.
!
!     LIW      The length of the array IWORK, as declared by the user.
!              (This will be checked by the solver.)
!
!     Note:  The work arrays must not be altered between calls to DLSODE
!     for the same problem, except possibly for the conditional and
!     optional inputs, and except for the last 3*NEQ words of RWORK.
!     The latter space is used for internal scratch space, and so is
!     available for use by the user outside DLSODE between calls, if
!     desired (but not for use by F or JAC).
!
!     JAC      The name of the user-supplied routine (MITER = 1 or 4) to
!              compute the Jacobian matrix, df/dy, as a function of the
!              scalar t and the vector y.  (See the MF description below
!              for MITER.)  It is to have the form
!
!                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
!                 DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
!
!              where NEQ, T, Y, ML, MU, and NROWPD are input and the
!              array PD is to be loaded with partial derivatives
!              (elements of the Jacobian matrix) on output.  PD must be
!              given a first dimension of NROWPD.  T and Y have the same
!              meaning as in subroutine F.
!
!              In the full matrix case (MITER = 1), ML and MU are
!              ignored, and the Jacobian is to be loaded into PD in
!              columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!
!              In the band matrix case (MITER = 4), the elements within
!              the band are to be loaded into PD in columnwise manner,
!              with diagonal lines of df/dy loaded into the rows of PD.
!              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML
!              and MU are the half-bandwidth parameters (see IWORK).
!              The locations in PD in the two triangular areas which
!              correspond to nonexistent matrix elements can be ignored
!              or loaded arbitrarily, as they are overwritten by DLSODE.
!
!              JAC need not provide df/dy exactly. A crude approximation
!              (possibly with a smaller bandwidth) will do.
!
!              In either case, PD is preset to zero by the solver, so
!              that only the nonzero elements need be loaded by JAC.
!              Each call to JAC is preceded by a call to F with the same
!              arguments NEQ, T, and Y. Thus to gain some efficiency,
!              intermediate quantities shared by both calculations may
!              be saved in a user COMMON block by F and not recomputed
!              by JAC, if desired.  Also, JAC may alter the Y array, if
!              desired.  JAC must be declared EXTERNAL in the calling
!              program.
!
!              Subroutine JAC may access user-defined quantities in
!              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
!              (dimensioned in JAC) and/or Y has length exceeding
!              NEQ(1).  See the descriptions of NEQ and Y above.
!
!     MF       The method flag.  Used only for input.  The legal values
!              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24,
!              and 25.  MF has decimal digits METH and MITER:
!                 MF = 10*METH + MITER .
!
!              METH indicates the basic linear multistep method:
!              1   Implicit Adams method.
!              2   Method based on backward differentiation formulas
!                  (BDF's).
!
!              MITER indicates the corrector iteration method:
!              0   Functional iteration (no Jacobian matrix is
!                  involved).
!              1   Chord iteration with a user-supplied full (NEQ by
!                  NEQ) Jacobian.
!              2   Chord iteration with an internally generated
!                  (difference quotient) full Jacobian (using NEQ
!                  extra calls to F per df/dy value).
!              3   Chord iteration with an internally generated
!                  diagonal Jacobian approximation (using one extra call
!                  to F per df/dy evaluation).
!              4   Chord iteration with a user-supplied banded Jacobian.
!              5   Chord iteration with an internally generated banded
!                  Jacobian (using ML + MU + 1 extra calls to F per
!                  df/dy evaluation).
!
!              If MITER = 1 or 4, the user must supply a subroutine JAC
!              (the name is arbitrary) as described above under JAC.
!              For other values of MITER, a dummy argument can be used.
!
!     Optional Inputs
!     ---------------
!     The following is a list of the optional inputs provided for in the
!     call sequence.  (See also Part 2.)  For each such input variable,
!     this table lists its name as used in this documentation, its
!     location in the call sequence, its meaning, and the default value.
!     The use of any of these inputs requires IOPT = 1, and in that case
!     all of these inputs are examined.  A value of zero for any of
!     these optional inputs will cause the default value to be used.
!     Thus to use a subset of the optional inputs, simply preload
!     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively,
!     and then set those of interest to nonzero values.
!
!     Name    Location   Meaning and default value
!     ------  ---------  -----------------------------------------------
!     H0      RWORK(5)   Step size to be attempted on the first step.
!                        The default value is determined by the solver.
!     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
!                        default value is infinite.
!     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
!                        default value is 0.  (This lower bound is not
!                        enforced on the final step before reaching
!                        TCRIT when ITASK = 4 or 5.)
!     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
!                        is 12 if METH = 1, and 5 if METH = 2. (See the
!                        MF description above for METH.)  If MAXORD
!                        exceeds the default value, it will be reduced
!                        to the default value.  If MAXORD is changed
!                        during the problem, it may cause the current
!                        order to be reduced.
!     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
!                        allowed during one call to the solver.  The
!                        default value is 500.
!     MXHNIL  IWORK(7)   Maximum number of messages printed (per
!                        problem) warning that T + H = T on a step
!                        (H = step size).  This must be positive to
!                        result in a nondefault value.  The default
!                        value is 10.
!
!     Optional Outputs
!     ----------------
!     As optional additional output from DLSODE, the variables listed
!     below are quantities related to the performance of DLSODE which
!     are available to the user.  These are communicated by way of the
!     work arrays, but also have internal mnemonic names as shown.
!     Except where stated otherwise, all of these outputs are defined on
!     any successful return from DLSODE, and on any return with ISTATE =
!     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3),
!     they will be unchanged from their existing values (if any), except
!     possibly for TOLSF, LENRW, and LENIW.  On any error return,
!     outputs relevant to the error will be defined, as noted below.
!
!     Name   Location   Meaning
!     -----  ---------  ------------------------------------------------
!     HU     RWORK(11)  Step size in t last used (successfully).
!     HCUR   RWORK(12)  Step size to be attempted on the next step.
!     TCUR   RWORK(13)  Current value of the independent variable which
!                       the solver has actually reached, i.e., the
!                       current internal mesh point in t. On output,
!                       TCUR will always be at least as far as the
!                       argument T, but may be farther (if interpolation
!                       was done).
!     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0,
!                       computed when a request for too much accuracy
!                       was detected (ISTATE = -3 if detected at the
!                       start of the problem, ISTATE = -2 otherwise).
!                       If ITOL is left unaltered but RTOL and ATOL are
!                       uniformly scaled up by a factor of TOLSF for the
!                       next call, then the solver is deemed likely to
!                       succeed.  (The user may also ignore TOLSF and
!                       alter the tolerance parameters in any other way
!                       appropriate.)
!     NST    IWORK(11)  Number of steps taken for the problem so far.
!     NFE    IWORK(12)  Number of F evaluations for the problem so far.
!     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU
!                       decompositions) for the problem so far.
!     NQU    IWORK(14)  Method order last used (successfully).
!     NQCUR  IWORK(15)  Order to be attempted on the next step.
!     IMXER  IWORK(16)  Index of the component of largest magnitude in
!                       the weighted local error vector ( e(i)/EWT(i) ),
!                       on an error return with ISTATE = -4 or -5.
!     LENRW  IWORK(17)  Length of RWORK actually required.  This is
!                       defined on normal returns and on an illegal
!                       input return for insufficient storage.
!     LENIW  IWORK(18)  Length of IWORK actually required.  This is
!                       defined on normal returns and on an illegal
!                       input return for insufficient storage.
!
!     The following two arrays are segments of the RWORK array which may
!     also be of interest to the user as optional outputs.  For each
!     array, the table below gives its internal name, its base address
!     in RWORK, and its description.
!
!     Name  Base address  Description
!     ----  ------------  ----------------------------------------------
!     YH    21            The Nordsieck history array, of size NYH by
!                         (NQCUR + 1), where NYH is the initial value of
!                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of
!                         YH contains HCUR**j/factorial(j) times the jth
!                         derivative of the interpolating polynomial
!                         currently representing the solution, evaluated
!                         at t = TCUR.
!     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated
!                         corrections on each step, scaled on output to
!                         represent the estimated local error in Y on
!                         the last step.  This is the vector e in the
!                         description of the error control.  It is
!                         defined only on successful return from DLSODE.
!
!
!                    Part 2.  Other Callable Routines
!                    --------------------------------
!
!     The following are optional calls which the user may make to gain
!     additional capabilities in conjunction with DLSODE.
!
!     Form of call              Function
!     ------------------------  ----------------------------------------
!     CALL XSETUN(LUN)          Set the logical unit number, LUN, for
!                               output of messages from DLSODE, if the
!                               default is not desired.  The default
!                               value of LUN is 6. This call may be made
!                               at any time and will take effect
!                               immediately.
!     CALL XSETF(MFLAG)         Set a flag to control the printing of
!                               messages by DLSODE.  MFLAG = 0 means do
!                               not print.  (Danger:  this risks losing
!                               valuable information.)  MFLAG = 1 means
!                               print (the default).  This call may be
!                               made at any time and will take effect
!                               immediately.
!     CALL DSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the
!                               internal COMMON blocks used by DLSODE
!                               (see Part 3 below).  RSAV must be a
!                               real array of length 218 or more, and
!                               ISAV must be an integer array of length
!                               37 or more.  JOB = 1 means save COMMON
!                               into RSAV/ISAV.  JOB = 2 means restore
!                               COMMON from same.  DSRCOM is useful if
!                               one is interrupting a run and restarting
!                               later, or alternating between two or
!                               more problems solved with DLSODE.
!     CALL DINTDY(,,,,,)        Provide derivatives of y, of various
!     (see below)               orders, at a specified point t, if
!                               desired.  It may be called only after a
!                               successful return from DLSODE.  Detailed
!                               instructions follow.
!
!     Detailed instructions for using DINTDY
!     --------------------------------------
!     The form of the CALL is:
!
!           CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
!     The input parameters are:
!
!     T          Value of independent variable where answers are
!                desired (normally the same as the T last returned by
!                DLSODE).  For valid results, T must lie between
!                TCUR - HU and TCUR.  (See "Optional Outputs" above
!                for TCUR and HU.)
!     K          Integer order of the derivative desired.  K must
!                satisfy 0 <= K <= NQCUR, where NQCUR is the current
!                order (see "Optional Outputs").  The capability
!                corresponding to K = 0, i.e., computing y(t), is
!                already provided by DLSODE directly.  Since
!                NQCUR >= 1, the first derivative dy/dt is always
!                available with DINTDY.
!     RWORK(21)  The base address of the history array YH.
!     NYH        Column length of YH, equal to the initial value of NEQ.
!
!     The output parameters are:
!
!     DKY        Real array of length NEQ containing the computed value
!                of the Kth derivative of y(t).
!     IFLAG      Integer flag, returned as 0 if K and T were legal,
!                -1 if K was illegal, and -2 if T was illegal.
!                On an error return, a message is also written.
!
!
!                          Part 3.  Common Blocks
!                          ----------------------
!
!     If DLSODE is to be used in an overlay situation, the user must
!     declare, in the primary overlay, the variables in:
!     (1) the call sequence to DLSODE,
!     (2) the internal COMMON block /DLS001/, of length 255
!         (218 double precision words followed by 37 integer words).
!
!     If DLSODE is used on a system in which the contents of internal
!     COMMON blocks are not preserved between calls, the user should
!     declare the above COMMON block in his main program to insure that
!     its contents are preserved.
!
!     If the solution of a given problem by DLSODE is to be interrupted
!     and then later continued, as when restarting an interrupted run or
!     alternating between two or more problems, the user should save,
!     following the return from the last DLSODE call prior to the
!     interruption, the contents of the call sequence variables and the
!     internal COMMON block, and later restore these values before the
!     next DLSODE call for that problem.   In addition, if XSETUN and/or
!     XSETF was called for non-default handling of error messages, then
!     these calls must be repeated.  To save and restore the COMMON
!     block, use subroutine DSRCOM (see Part 2 above).
!
!
!              Part 4.  Optionally Replaceable Solver Routines
!              -----------------------------------------------
!
!     Below are descriptions of two routines in the DLSODE package which
!     relate to the measurement of errors.  Either routine can be
!     replaced by a user-supplied version, if desired.  However, since
!     such a replacement may have a major impact on performance, it
!     should be done only when absolutely necessary, and only with great
!     caution.  (Note:  The means by which the package version of a
!     routine is superseded by the user's version may be system-
!     dependent.)
!
!     DEWSET
!     ------
!     The following subroutine is called just before each internal
!     integration step, and sets the array of error weights, EWT, as
!     described under ITOL/RTOL/ATOL above:
!
!           SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
!
!     where NEQ, ITOL, RTOL, and ATOL are as in the DLSODE call
!     sequence, YCUR contains the current dependent variable vector,
!     and EWT is the array of weights set by DEWSET.
!
!     If the user supplies this subroutine, it must return in EWT(i)
!     (i = 1,...,NEQ) a positive quantity suitable for comparing errors
!     in Y(i) to.  The EWT array returned by DEWSET is passed to the
!     DVNORM routine (see below), and also used by DLSODE in the
!     computation of the optional output IMXER, the diagonal Jacobian
!     approximation, and the increments for difference quotient
!     Jacobians.
!
!     In the user-supplied version of DEWSET, it may be desirable to use
!     the current values of derivatives of y. Derivatives up to order NQ
!     are available from the history array YH, described above under
!     optional outputs.  In DEWSET, YH is identical to the YCUR array,
!     extended to NQ + 1 columns with a column length of NYH and scale
!     factors of H**j/factorial(j).  On the first call for the problem,
!     given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
!     NYH is the initial value of NEQ.  The quantities NQ, H, and NST
!     can be obtained by including in SEWSET the statements:
!           DOUBLE PRECISION RLS
!           COMMON /DLS001/ RLS(218),ILS(37)
!           NQ = ILS(33)
!           NST = ILS(34)
!           H = RLS(212)
!     Thus, for example, the current value of dy/dt can be obtained as
!     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary
!     when NST = 0).
!
!     DVNORM
!     ------
!     DVNORM is a real function routine which computes the weighted
!     root-mean-square norm of a vector v:
!
!        d = DVNORM (n, v, w)
!
!     where:
!     n = the length of the vector,
!     v = real array of length n containing the vector,
!     w = real array of length n containing weights,
!     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ).
!
!     DVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where
!     EWT is as set by subroutine DEWSET.
!
!     If the user supplies this function, it should return a nonnegative
!     value of DVNORM suitable for use in the error control in DLSODE.
!     None of the arguments should be altered by DVNORM.  For example, a
!     user-supplied DVNORM routine might:
!     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
!     - Ignore some components of v in the norm, with the effect of
!       suppressing the error control on those components of Y.
!  ---------------------------------------------------------------------
!***ROUTINES CALLED  DEWSET, DINTDY, DUMACH, DSTODE, DVNORM, XERRWD
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYYYMMDD)
! 19791129  DATE WRITTEN
! 19791213  Minor changes to declarations; DELP init. in STODE.
! 19800118  Treat NEQ as array; integer declarations added throughout;
!           minor changes to prologue.
! 19800306  Corrected TESCO(1,NQP1) setting in CFODE.
! 19800519  Corrected access of YH on forced order reduction;
!           numerous corrections to prologues and other comments.
! 19800617  In main driver, added loading of SQRT(UROUND) in RWORK;
!           minor corrections to main prologue.
! 19800923  Added zero initialization of HU and NQU.
! 19801218  Revised XERRWD routine; minor corrections to main prologue.
! 19810401  Minor changes to comments and an error message.
! 19810814  Numerous revisions: replaced EWT by 1/EWT; used flags
!           JCUR, ICF, IERPJ, IERSL between STODE and subordinates;
!           added tuning parameters CCMAX, MAXCOR, MSBP, MXNCF;
!           reorganized returns from STODE; reorganized type decls.;
!           fixed message length in XERRWD; changed default LUNIT to 6;
!           changed Common lengths; changed comments throughout.
! 19870330  Major update by ACH: corrected comments throughout;
!           removed TRET from Common; rewrote EWSET with 4 loops;
!           fixed t test in INTDY; added Cray directives in STODE;
!           in STODE, fixed DELP init. and logic around PJAC call;
!           combined routines to save/restore Common;
!           passed LEVEL = 0 in error message calls (except run abort).
! 19890426  Modified prologue to SLATEC/LDOC format.  (FNF)
! 19890501  Many improvements to prologue.  (FNF)
! 19890503  A few final corrections to prologue.  (FNF)
! 19890504  Minor cosmetic changes.  (FNF)
! 19890510  Corrected description of Y in Arguments section.  (FNF)
! 19890517  Minor corrections to prologue.  (FNF)
! 19920514  Updated with prologue edited 891025 by G. Shaw for manual.
! 19920515  Converted source lines to upper case.  (FNF)
! 19920603  Revised XERRWD calls using mixed upper-lower case.  (ACH)
! 19920616  Revised prologue comment regarding CFT.  (ACH)
! 19921116  Revised prologue comments regarding Common.  (ACH).
! 19930326  Added comment about non-reentrancy.  (FNF)
! 19930723  Changed D1MACH to DUMACH. (FNF)
! 19930801  Removed ILLIN and NTREP from Common (affects driver logic);
!           minor changes to prologue and internal comments;
!           changed Hollerith strings to quoted strings;
!           changed internal comments to mixed case;
!           replaced XERRWD with new version using character type;
!           changed dummy dimensions from 1 to *. (ACH)
! 19930809  Changed to generic intrinsic names; changed names of
!           subprograms and Common blocks to DLSODE etc. (ACH)
! 19930929  Eliminated use of REAL intrinsic; other minor changes. (ACH)
! 20010412  Removed all 'own' variables from Common block /DLS001/
!           (affects declarations in 6 routines). (ACH)
! 20010509  Minor corrections to prologue. (ACH)
! 20031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
! 20031112  Added SAVE statements for data-loaded constants.
!
!***END PROLOGUE  DLSODE
!
!*Internal Notes:
!
! Other Routines in the DLSODE Package.
!
! In addition to Subroutine DLSODE, the DLSODE package includes the
! following subroutines and function routines:
!  DINTDY   computes an interpolated value of the y vector at t = TOUT.
!  DSTODE   is the core integrator, which does one step of the
!           integration and the associated error control.
!  DCFODE   sets all method coefficients and test constants.
!  DPREPJ   computes and preprocesses the Jacobian matrix J = df/dy
!           and the Newton iteration matrix P = I - h*l0*J.
!  DSOLSY   manages solution of linear system in chord iteration.
!  DEWSET   sets the error weight vector EWT before each step.
!  DVNORM   computes the weighted R.M.S. norm of a vector.
!  DSRCOM   is a user-callable routine to save and restore
!           the contents of the internal Common block.
!  DGEFA and DGESL   are routines from LINPACK for solving full
!           systems of linear algebraic equations.
!  DGBFA and DGBSL   are routines from LINPACK for solving banded
!           linear systems.
!  DUMACH   computes the unit roundoff in a machine-independent manner.
!  XERRWD, XSETUN, XSETF, IXSAV, IUMACH   handle the printing of all
!           error messages and warnings.  XERRWD is machine-dependent.
! Note: DVNORM, DUMACH, IXSAV, and IUMACH are function routines.
! All the others are subroutines.
!
!**End
!
!  Declare externals.
!      EXTERNAL dprepj, dsolsy
!      DOUBLE PRECISION dumach, dvnorm
!
!  Declare all other variables.
      INTEGER init, mxstep, mxhnil, nhnil, nslast, nyh, iowns, icf, ierpj, &
        iersl, jcur, jstart, kflag, l, lyh, lewt, lacor, lsavf, lwm, liwm, &
        meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, i1, i2, iflag, imxer, kgo, lf0, leniw, lenrw, lenwm, ml, &
        mord, mu, mxhnl0, mxstp0
      DOUBLE PRECISION rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      DOUBLE PRECISION atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, tcrit, &
        tdist, tnext, tol, tolsf, tp, size, sum, w0
      DIMENSION mord(2)
      LOGICAL ihit
!      CHARACTER *80 msg
      CHARACTER(LEN=80) msg
      SAVE mord, mxstp0, mxhnl0
!-----------------------------------------------------------------------
! The following internal Common block contains
! (a) variables which are local to any subroutine but whose values must
!     be preserved between calls to the routine ("own" variables), and
! (b) variables which are communicated between subroutines.
! The block DLS001 is declared in subroutines DLSODE, DINTDY, DSTODE,
! DPREPJ, and DSOLSY.
! Groups of variables are replaced by dummy arrays in the Common
! declarations in routines where those variables are not used.
!-----------------------------------------------------------------------
      COMMON /dls001/rowns(209), ccmax, el0, h, hmin, hmxi, hu, rc, tn, &
        uround, init, mxstep, mxhnil, nhnil, nslast, nyh, iowns(6), icf, &
        ierpj, iersl, jcur, jstart, kflag, l, lyh, lewt, lacor, lsavf, lwm, &
        liwm, meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, &
        nqu
!
      DATA mord(1), mord(2)/12, 5/, mxstp0/500/, mxhnl0/10/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .GT. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DLSODE
      IF (istate<1 .OR. istate>3) GO TO 601
      IF (itask<1 .OR. itask>5) GO TO 602
      IF (istate==1) GO TO 10
      IF (init==0) GO TO 603
      IF (istate==2) GO TO 200
      GO TO 20
10    init = 0
      IF (tout==t) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all inputs and various initializations.
!
! First check legality of the non-optional inputs NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
20    IF (neq(1)<=0) GO TO 604
      IF (istate==1) GO TO 25
      IF (neq(1)>n) GO TO 605
25    n = neq(1)
      IF (itol<1 .OR. itol>4) GO TO 606
      IF (iopt<0 .OR. iopt>1) GO TO 607
      meth = mf/10
      miter = mf - 10*meth
      IF (meth<1 .OR. meth>2) GO TO 608
      IF (miter<0 .OR. miter>5) GO TO 608
      IF (miter<=3) GO TO 30
      ml = iwork(1)
      mu = iwork(2)
      IF (ml<0 .OR. ml>=n) GO TO 609
      IF (mu<0 .OR. mu>=n) GO TO 610
30    CONTINUE
! Next process and check the optional inputs. --------------------------
      IF (iopt==1) GO TO 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      IF (istate==1) h0 = 0.0D0
      hmxi = 0.0D0
      hmin = 0.0D0
      GO TO 60
40    maxord = iwork(5)
      IF (maxord<0) GO TO 611
      IF (maxord==0) maxord = 100
      maxord = min(maxord, mord(meth))
      mxstep = iwork(6)
      IF (mxstep<0) GO TO 612
      IF (mxstep==0) mxstep = mxstp0
      mxhnil = iwork(7)
      IF (mxhnil<0) GO TO 613
      IF (mxhnil==0) mxhnil = mxhnl0
      IF (istate/=1) GO TO 50
      h0 = rwork(5)
      IF ((tout-t)*h0<0.0D0) GO TO 614
50    hmax = rwork(6)
      IF (hmax<0.0D0) GO TO 615
      hmxi = 0.0D0
      IF (hmax>0.0D0) hmxi = 1.0D0/hmax
      hmin = rwork(7)
      IF (hmin<0.0D0) GO TO 616
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.
! Pointers to segments of RWORK and IWORK are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
!-----------------------------------------------------------------------
60    lyh = 21
      IF (istate==1) nyh = n
      lwm = lyh + (maxord+1)*nyh
      IF (miter==0) lenwm = 0
      IF (miter==1 .OR. miter==2) lenwm = n*n + 2
      IF (miter==3) lenwm = n + 2
      IF (miter>=4) lenwm = (2*ml+mu+1)*n + 2
      lewt = lwm + lenwm
      lsavf = lewt + n
      lacor = lsavf + n
      lenrw = lacor + n - 1
      iwork(17) = lenrw
      liwm = 1
      leniw = 20 + n
      IF (miter==0 .OR. miter==3) leniw = 20
      iwork(18) = leniw
      IF (lenrw>lrw) GO TO 617
      IF (leniw>liw) GO TO 618
! Check RTOL and ATOL for legality. ------------------------------------
      rtoli = rtol(1)
      atoli = atol(1)
      DO i = 1, n
        IF (itol>=3) rtoli = rtol(i)
        IF (itol==2 .OR. itol==4) atoli = atol(i)
        IF (rtoli<0.0D0) GO TO 619
        IF (atoli<0.0D0) GO TO 620
      ENDDO
      IF (istate==1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to DSTODE. -------
      jstart = -1
      IF (nq<=maxord) GO TO 90
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO i = 1, n
        rwork(i+lsavf-1) = rwork(i+lwm-1)
      ENDDO
! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
90    IF (miter>0) rwork(lwm) = sqrt(uround)
      IF (n==nyh) GO TO 200
! NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      i1 = lyh + l*nyh
      i2 = lyh + (maxord+1)*nyh - 1
      IF (i1>i2) GO TO 200
      DO i = i1, i2
        rwork(i) = 0.0D0
      ENDDO
      GO TO 200
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
100   uround = dumach()
      tn = t
      IF (itask/=4 .AND. itask/=5) GO TO 110
      tcrit = rwork(1)
      IF ((tcrit-tout)*(tout-t)<0.0D0) GO TO 625
      IF (h0/=0.0D0 .AND. (t+h0-tcrit)*h0>0.0D0) h0 = tcrit - t
110   jstart = 0
      IF (miter>0) rwork(lwm) = sqrt(uround)
      nhnil = 0
      nst = 0
      nje = 0
      nslast = 0
      hu = 0.0D0
      nqu = 0
      ccmax = 0.3D0
      maxcor = 3
      msbp = 20
      mxncf = 10
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      lf0 = lyh + nyh
      CALL f(neq, t, y, rwork(lf0))
      nfe = 1
! Load the initial value vector in YH. ---------------------------------
      DO i = 1, n
        rwork(i+lyh-1) = y(i)
      ENDDO
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      nq = 1
      h = 1.0D0
      CALL dewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      DO i = 1, n
        IF (rwork(i+lewt-1)<=0.0D0) GO TO 621
        rwork(i+lewt-1) = 1.0D0/rwork(i+lewt-1)
      ENDDO
!-----------------------------------------------------------------------
! The coding below computes the step size, H0, to be attempted on the
! first step, unless the user has supplied a value for this.
! First check that TOUT - T differs significantly from zero.
! A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
! if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
! so as to be between 100*UROUND and 1.0E-3.
! Then the computed value H0 is given by..
!                                      NEQ
!   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
!                                       1
! where   w0     = MAX ( ABS(T), ABS(TOUT) ),
!         f(i)   = i-th component of initial value of f,
!         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
! The sign of H0 is inferred from the initial values of TOUT and T.
!-----------------------------------------------------------------------
      IF (h0/=0.0D0) GO TO 180
      tdist = abs(tout-t)
      w0 = max(abs(t), abs(tout))
      IF (tdist<2.0D0*uround*w0) GO TO 622
      tol = rtol(1)
      IF (itol<=2) GO TO 140
      DO i = 1, n
        tol = max(tol, rtol(i))
      END DO
140   IF (tol>0.0D0) GO TO 160
      atoli = atol(1)
      DO i = 1, n
        IF (itol==2 .OR. itol==4) atoli = atol(i)
        ayi = abs(y(i))
        IF (ayi/=0.0D0) tol = max(tol, atoli/ayi)
      ENDDO
160   tol = max(tol, 100.0D0*uround)
      tol = min(tol, 0.001D0)
      sum = dvnorm(n, rwork(lf0), rwork(lewt))
      sum = 1.0D0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0D0/sqrt(sum)
      h0 = min(h0, tdist)
      h0 = sign(h0, tout-t)
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180   rh = abs(h0)*hmxi
      IF (rh>1.0D0) h0 = h0/rh
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      h = h0
      DO i = 1, n
        rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      ENDDO
      GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
200   nslast = nst
!      GO TO (210, 250, 220, 230, 240), itask
      SELECT CASE(itask)
        CASE(1) 
          GO TO 210
        CASE(2)
          GO TO 250
        CASE(3)
          GO TO 220
        CASE(4)
          GO TO 230
        CASE(5) 
          GO TO 240
      END SELECT 
210   IF ((tn-tout)*h<0.0D0) GO TO 250
      CALL dintdy(tout, 0, rwork(lyh), nyh, y, iflag)
      IF (iflag/=0) GO TO 627
      t = tout
      GO TO 420
220   tp = tn - hu*(1.0D0+100.0D0*uround)
      IF ((tp-tout)*h>0.0D0) GO TO 623
      IF ((tn-tout)*h<0.0D0) GO TO 250
      GO TO 400
230   tcrit = rwork(1)
      IF ((tn-tcrit)*h>0.0D0) GO TO 624
      IF ((tcrit-tout)*h<0.0D0) GO TO 625
      IF ((tn-tout)*h<0.0D0) GO TO 245
      CALL dintdy(tout, 0, rwork(lyh), nyh, y, iflag)
      IF (iflag/=0) GO TO 627
      t = tout
      GO TO 420
240   tcrit = rwork(1)
      IF ((tn-tcrit)*h>0.0D0) GO TO 624
245   hmx = abs(tn) + abs(h)
      ihit = abs(tn-tcrit) <= 100.0D0*uround*hmx
      IF (ihit) GO TO 400
      tnext = tn + h*(1.0D0+4.0D0*uround)
      IF ((tnext-tcrit)*h<=0.0D0) GO TO 250
      h = (tcrit-tn)*(1.0D0-4.0D0*uround)
      IF (istate==2) jstart = -2
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DSTODE.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
250   CONTINUE
      IF ((nst-nslast)>=mxstep) GO TO 500
      CALL dewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      DO i = 1, n
        IF (rwork(i+lewt-1)<=0.0D0) GO TO 510
        rwork(i+lewt-1) = 1.0D0/rwork(i+lewt-1)
      ENDDO
270   tolsf = uround*dvnorm(n, rwork(lyh), rwork(lewt))
      IF (tolsf<=1.0D0) GO TO 280
      tolsf = tolsf*2.0D0
      IF (nst==0) GO TO 626
      GO TO 520
280   IF ((tn+h)/=tn) GO TO 290
      nhnil = nhnil + 1
      IF (nhnil>mxhnil) GO TO 290
      msg = 'DLSODE-  Warning..internal T (=R1) and H (=R2) are'
      CALL xerrwd(msg, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      such that in the machine, T + H = T on the next step  '
      CALL xerrwd(msg, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      (H = step size). Solver will continue anyway'
      CALL xerrwd(msg, 50, 101, 0, 0, 0, 0, 2, tn, h)
      IF (nhnil<mxhnil) GO TO 290
      msg = 'DLSODE-  Above warning has been issued I1 times.  '
      CALL xerrwd(msg, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      It will not be issued again for this problem'
      CALL xerrwd(msg, 50, 102, 0, 1, mxhnil, 0, 0, 0.0D0, 0.0D0)
290   CONTINUE
!-----------------------------------------------------------------------
!  CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPREPJ,DSOLSY)
!-----------------------------------------------------------------------
      CALL dstode(neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), &
        rwork(lsavf), rwork(lacor), rwork(lwm), iwork(liwm), f, jac, dprepj, &
        dsolsy)
      kgo = 1 - kflag
!      GO TO (300, 530, 540), kgo
      SELECT CASE(kgo)
        CASE(1)
          GO TO 300
        CASE(2)
          GO TO 530
        CASE(3) 
          GO TO 540
      END SELECT
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
300   init = 1
!      GO TO (310, 400, 330, 340, 350), itask
      SELECT CASE(itask)
        CASE(1)
          GO TO 310
        CASE(2)
          GO TO 400
        CASE(3)
          GO TO 330
        CASE(4)
          GO TO 340
        CASE(5)
          GO TO 350
      END SELECT
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310   IF ((tn-tout)*h<0.0D0) GO TO 250
      CALL dintdy(tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330   IF ((tn-tout)*h>=0.0D0) GO TO 400
      GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340   IF ((tn-tout)*h<0.0D0) GO TO 345
      CALL dintdy(tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      GO TO 420
345   hmx = abs(tn) + abs(h)
      ihit = abs(tn-tcrit) <= 100.0D0*uround*hmx
      IF (ihit) GO TO 400
      tnext = tn + h*(1.0D0+4.0D0*uround)
      IF ((tnext-tcrit)*h<=0.0D0) GO TO 250
      h = (tcrit-tn)*(1.0D0-4.0D0*uround)
      jstart = -2
      GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350   hmx = abs(tn) + abs(h)
      ihit = abs(tn-tcrit) <= 100.0D0*uround*hmx
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from DLSODE.
! If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional outputs are loaded into the
! work arrays before returning.
!-----------------------------------------------------------------------
400   DO i = 1, n
        y(i) = rwork(i+lyh-1)
      ENDDO
      t = tn
      IF (itask/=4 .AND. itask/=5) GO TO 420
      IF (ihit) t = tcrit
420   istate = 2
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! If there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH and T is set to TN.  The optional outputs
! are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
500   msg = 'DLSODE-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL xerrwd(msg, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      taken on this call before reaching TOUT     '
      CALL xerrwd(msg, 50, 201, 0, 1, mxstep, 0, 1, tn, 0.0D0)
      istate = -1
      GO TO 580
! EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
510   ewti = rwork(lewt+i-1)
      msg = 'DLSODE-  At T (=R1), EWT(I1) has become R2 .LE. 0.'
      CALL xerrwd(msg, 50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      GO TO 580
! Too much accuracy requested for machine precision. -------------------
520   msg = 'DLSODE-  At T (=R1), too much accuracy requested  '
      CALL xerrwd(msg, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      for precision of machine..  see TOLSF (=R2) '
      CALL xerrwd(msg, 50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530   msg = 'DLSODE-  At T(=R1) and step size H(=R2), the error'
      CALL xerrwd(msg, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL xerrwd(msg, 50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540   msg = 'DLSODE-  At T (=R1) and step size H (=R2), the    '
      CALL xerrwd(msg, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      corrector convergence failed repeatedly     '
      CALL xerrwd(msg, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      or with ABS(H) = HMIN   '
      CALL xerrwd(msg, 30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
! Compute IMXER if relevant. -------------------------------------------
560   big = 0.0D0
      imxer = 1
      DO i = 1, n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        IF (big>=size) GO TO 570
        big = size
        imxer = i
570   ENDDO
      iwork(16) = imxer
! Set Y vector, T, and optional outputs. -------------------------------
580   DO i = 1, n
        y(i) = rwork(i+lyh-1)
      ENDDO
      t = tn
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.  If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
601   msg = 'DLSODE-  ISTATE (=I1) illegal '
      CALL xerrwd(msg, 30, 1, 0, 1, istate, 0, 0, 0.0D0, 0.0D0)
      IF (istate<0) GO TO 800
      GO TO 700
602   msg = 'DLSODE-  ITASK (=I1) illegal  '
      CALL xerrwd(msg, 30, 2, 0, 1, itask, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
603   msg = 'DLSODE-  ISTATE .GT. 1 but DLSODE not initialized '
      CALL xerrwd(msg, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
604   msg = 'DLSODE-  NEQ (=I1) .LT. 1     '
      CALL xerrwd(msg, 30, 4, 0, 1, neq(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
605   msg = 'DLSODE-  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL xerrwd(msg, 50, 5, 0, 2, n, neq(1), 0, 0.0D0, 0.0D0)
      GO TO 700
606   msg = 'DLSODE-  ITOL (=I1) illegal   '
      CALL xerrwd(msg, 30, 6, 0, 1, itol, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
607   msg = 'DLSODE-  IOPT (=I1) illegal   '
      CALL xerrwd(msg, 30, 7, 0, 1, iopt, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
608   msg = 'DLSODE-  MF (=I1) illegal     '
      CALL xerrwd(msg, 30, 8, 0, 1, mf, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
609   msg = 'DLSODE-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL xerrwd(msg, 50, 9, 0, 2, ml, neq(1), 0, 0.0D0, 0.0D0)
      GO TO 700
610   msg = 'DLSODE-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL xerrwd(msg, 50, 10, 0, 2, mu, neq(1), 0, 0.0D0, 0.0D0)
      GO TO 700
611   msg = 'DLSODE-  MAXORD (=I1) .LT. 0  '
      CALL xerrwd(msg, 30, 11, 0, 1, maxord, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
612   msg = 'DLSODE-  MXSTEP (=I1) .LT. 0  '
      CALL xerrwd(msg, 30, 12, 0, 1, mxstep, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
613   msg = 'DLSODE-  MXHNIL (=I1) .LT. 0  '
      CALL xerrwd(msg, 30, 13, 0, 1, mxhnil, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
614   msg = 'DLSODE-  TOUT (=R1) behind T (=R2)      '
      CALL xerrwd(msg, 40, 14, 0, 0, 0, 0, 2, tout, t)
      msg = '      Integration direction is given by H0 (=R1)  '
      CALL xerrwd(msg, 50, 14, 0, 0, 0, 0, 1, h0, 0.0D0)
      GO TO 700
615   msg = 'DLSODE-  HMAX (=R1) .LT. 0.0  '
      CALL xerrwd(msg, 30, 15, 0, 0, 0, 0, 1, hmax, 0.0D0)
      GO TO 700
616   msg = 'DLSODE-  HMIN (=R1) .LT. 0.0  '
      CALL xerrwd(msg, 30, 16, 0, 0, 0, 0, 1, hmin, 0.0D0)
      GO TO 700
617   CONTINUE
      msg = 'DLSODE-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL xerrwd(msg, 60, 17, 0, 2, lenrw, lrw, 0, 0.0D0, 0.0D0)
      GO TO 700
618   CONTINUE
      msg = 'DLSODE-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL xerrwd(msg, 60, 18, 0, 2, leniw, liw, 0, 0.0D0, 0.0D0)
      GO TO 700
619   msg = 'DLSODE-  RTOL(I1) is R1 .LT. 0.0        '
      CALL xerrwd(msg, 40, 19, 0, 1, i, 0, 1, rtoli, 0.0D0)
      GO TO 700
620   msg = 'DLSODE-  ATOL(I1) is R1 .LT. 0.0        '
      CALL xerrwd(msg, 40, 20, 0, 1, i, 0, 1, atoli, 0.0D0)
      GO TO 700
621   ewti = rwork(lewt+i-1)
      msg = 'DLSODE-  EWT(I1) is R1 .LE. 0.0         '
      CALL xerrwd(msg, 40, 21, 0, 1, i, 0, 1, ewti, 0.0D0)
      GO TO 700
622   CONTINUE
      msg = 'DLSODE-  TOUT (=R1) too close to T(=R2) to start integration'
      CALL xerrwd(msg, 60, 22, 0, 0, 0, 0, 2, tout, t)
      GO TO 700
623   CONTINUE
      msg = 'DLSODE-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL xerrwd(msg, 60, 23, 0, 1, itask, 0, 2, tout, tp)
      GO TO 700
624   CONTINUE
      msg = 'DLSODE-  ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL xerrwd(msg, 60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      GO TO 700
625   CONTINUE
      msg = 'DLSODE-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL xerrwd(msg, 60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      GO TO 700
626   msg = 'DLSODE-  At start of problem, too much accuracy   '
      CALL xerrwd(msg, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      msg = '      requested for precision of machine..  See TOLSF (=R1) '
      CALL xerrwd(msg, 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0D0)
      rwork(14) = tolsf
      GO TO 700
627   msg = 'DLSODE-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL xerrwd(msg, 50, 27, 0, 1, itask, 0, 1, tout, 0.0D0)
!
700   istate = -3
      RETURN
!
800   msg = 'DLSODE-  Run aborted.. apparent infinite loop     '
      CALL xerrwd(msg, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
!----------------------- END OF SUBROUTINE DLSODE ----------------------
    END SUBROUTINE dlsode

!DECK DUMACH
    DOUBLE PRECISION FUNCTION dumach()
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
      DOUBLE PRECISION u, comp
!***FIRST EXECUTABLE STATEMENT  DUMACH
      u = 1.0D0
10    u = u*0.5D0
      CALL dumsum(1.0D0, u, comp)
      IF (comp/=1.0D0) GO TO 10
      dumach = u*2.0D0
      RETURN
!----------------------- End of Function DUMACH ------------------------
    END FUNCTION dumach
    SUBROUTINE dumsum(a, b, c)
!     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION a, b, c

      c = a + b
      RETURN
    END SUBROUTINE dumsum
!DECK DCFODE
    SUBROUTINE dcfode(meth, elco, tesco)
!***BEGIN PROLOGUE  DCFODE
!***SUBSIDIARY
!***PURPOSE  Set ODE integrator coefficients.
!***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DCFODE is called by the integrator routine to set coefficients
!  needed there.  The coefficients for the current method, as
!  given by the value of METH, are set for all orders and saved.
!  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
!  (A smaller value of the maximum order is also allowed.)
!  DCFODE is called once at the beginning of the problem,
!  and is not called again unless and until METH is changed.
!
!  The ELCO array contains the basic method coefficients.
!  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
!  order nq are stored in ELCO(i,nq).  They are given by a genetrating
!  polynomial, i.e.,
!      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
!  For the implicit Adams methods, l(x) is given by
!      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
!  For the BDF methods, l(x) is given by
!      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
!  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
!  The TESCO array contains test constants used for the
!  local error test and the selection of step size and/or order.
!  At order nq, TESCO(k,nq) is used for the selection of step
!  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
!  nq + 1 if k = 3.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DCFODE
!**End
      INTEGER meth
      INTEGER i, ib, nq, nqm1, nqp1
      DOUBLE PRECISION elco, tesco
      DOUBLE PRECISION agamq, fnq, fnqm1, pc, pint, ragq, rqfac, rq1fac, &
        tsign, xpin
      DIMENSION elco(13, 12), tesco(3, 12)
      DIMENSION pc(12)
!
!***FIRST EXECUTABLE STATEMENT  DCFODE
!      GO TO (100, 200), meth
      SELECT CASE(meth)
        CASE(1)
          GO TO 100
        CASE(2)
          GO TO 200
      END SELECT
!
100   elco(1, 1) = 1.0D0
      elco(2, 1) = 1.0D0
      tesco(1, 1) = 0.0D0
      tesco(2, 1) = 2.0D0
      tesco(1, 2) = 1.0D0
      tesco(3, 12) = 0.0D0
      pc(1) = 1.0D0
      rqfac = 1.0D0
      DO nq = 2, 12
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq-1).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
        rq1fac = rqfac
        rqfac = rqfac/nq
        nqm1 = nq - 1
        fnqm1 = nqm1
        nqp1 = nq + 1
! Form coefficients of p(x)*(x+nq-1). ----------------------------------
        pc(nq) = 0.0D0
        DO ib = 1, nqm1
          i = nqp1 - ib
          pc(i) = pc(i-1) + fnqm1*pc(i)
        ENDDO
        pc(1) = fnqm1*pc(1)
! Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        pint = pc(1)
        xpin = pc(1)/2.0D0
        tsign = 1.0D0
        DO i = 2, nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/i
          xpin = xpin + tsign*pc(i)/(i+1)
        ENDDO
! Store coefficients in ELCO and TESCO. --------------------------------
        elco(1, nq) = pint*rq1fac
        elco(2, nq) = 1.0D0
        DO i = 2, nq
          elco(i+1, nq) = rq1fac*pc(i)/i
        ENDDO
        agamq = rqfac*xpin
        ragq = 1.0D0/agamq
        tesco(2, nq) = ragq
        IF (nq<12) tesco(1, nqp1) = ragq*rqfac/nqp1
        tesco(3, nqm1) = ragq
      ENDDO
      RETURN
!
200   pc(1) = 1.0D0
      rq1fac = 1.0D0
      DO nq = 1, 5
!-----------------------------------------------------------------------
! The PC array will contain the coefficients of the polynomial
!     p(x) = (x+1)*(x+2)*...*(x+nq).
! Initially, p(x) = 1.
!-----------------------------------------------------------------------
        fnq = nq
        nqp1 = nq + 1
! Form coefficients of p(x)*(x+nq). ------------------------------------
        pc(nqp1) = 0.0D0
        DO ib = 1, nq
          i = nq + 2 - ib
          pc(i) = pc(i-1) + fnq*pc(i)
        ENDDO
        pc(1) = fnq*pc(1)
! Store coefficients in ELCO and TESCO. --------------------------------
        DO i = 1, nqp1
          elco(i, nq) = pc(i)/pc(2)
        ENDDO
        elco(2, nq) = 1.0D0
        tesco(1, nq) = rq1fac
        tesco(2, nq) = nqp1/elco(1, nq)
        tesco(3, nq) = (nq+2)/elco(1, nq)
        rq1fac = rq1fac/fnq
      ENDDO
      RETURN
!----------------------- END OF SUBROUTINE DCFODE ----------------------
    END SUBROUTINE dcfode
!DECK DINTDY
    SUBROUTINE dintdy(t, k, yh, nyh, dky, iflag)
!***BEGIN PROLOGUE  DINTDY
!***SUBSIDIARY
!***PURPOSE  Interpolate solution derivatives.
!***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DINTDY computes interpolated values of the K-th derivative of the
!  dependent variable vector y, and stores it in DKY.  This routine
!  is called within the package with K = 0 and T = TOUT, but may
!  also be called by the user for any K up to the current order.
!  (See detailed instructions in the usage documentation.)
!
!  The computed values in DKY are gotten by interpolation using the
!  Nordsieck history array YH.  This array corresponds uniquely to a
!  vector-valued polynomial of degree NQCUR or less, and DKY is set
!  to the K-th derivative of this polynomial at T.
!  The formula for DKY is:
!               q
!   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
!              j=K
!  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
!  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
!  communicated by COMMON.  The above sum is done in reverse order.
!  IFLAG is returned negative if either K or T is out of bounds.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  XERRWD
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   050427  Corrected roundoff decrement in TP. (ACH)
!***END PROLOGUE  DINTDY
!**End
      INTEGER k, nyh, iflag
      DOUBLE PRECISION t, yh, dky
      DIMENSION yh(nyh, *), dky(*)
      INTEGER iownd, iowns, icf, ierpj, iersl, jcur, jstart, kflag, l, lyh, &
        lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      DOUBLE PRECISION rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      COMMON /dls001/rowns(209), ccmax, el0, h, hmin, hmxi, hu, rc, tn, &
        uround, iownd(6), iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, &
        lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, ic, j, jb, jb2, jj, jj1, jp1
      DOUBLE PRECISION c, r, s, tp
!      CHARACTER *80 msg
      CHARACTER(LEN=80) msg
!
!***FIRST EXECUTABLE STATEMENT  DINTDY
      iflag = 0
      IF (k<0 .OR. k>nq) GO TO 80
      tp = tn - hu - 100.0D0*uround*sign(abs(tn)+abs(hu), hu)
      IF ((t-tp)*(t-tn)>0.0D0) GO TO 90
!
      s = (t-tn)/h
      ic = 1
      IF (k==0) GO TO 15
      jj1 = l - k
      DO jj = jj1, nq
        ic = ic*jj
      ENDDO
15    c = ic
      DO i = 1, n
        dky(i) = c*yh(i, l)
      ENDDO
      IF (k==nq) GO TO 55
      jb2 = nq - k
      DO jb = 1, jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        IF (k==0) GO TO 35
        jj1 = jp1 - k
        DO jj = jj1, j
          ic = ic*jj
        ENDDO
35      c = ic
        DO i = 1, n
          dky(i) = c*yh(i, jp1) + s*dky(i)
        ENDDO
      ENDDO
      IF (k==0) RETURN
55    r = h**(-k)
      DO i = 1, n
        dky(i) = r*dky(i)
      ENDDO
      RETURN
!
80    msg = 'DINTDY-  K (=I1) illegal      '
      CALL xerrwd(msg, 30, 51, 0, 1, k, 0, 0, 0.0D0, 0.0D0)
      iflag = -1
      RETURN
90    msg = 'DINTDY-  T (=R1) illegal      '
      CALL xerrwd(msg, 30, 52, 0, 0, 0, 0, 1, t, 0.0D0)
      msg = '      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL xerrwd(msg, 60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      RETURN
!----------------------- END OF SUBROUTINE DINTDY ----------------------
    END SUBROUTINE dintdy
!DECK DPREPJ
    SUBROUTINE dprepj(neq, y, yh, nyh, ewt, ftem, savf, wm, iwm, f, jac)
!***BEGIN PROLOGUE  DPREPJ
!***SUBSIDIARY
!***PURPOSE  Compute and process Newton iteration matrix.
!***TYPE      DOUBLE PRECISION (SPREPJ-S, DPREPJ-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DPREPJ is called by DSTODE to compute and process the matrix
!  P = I - h*el(1)*J , where J is an approximation to the Jacobian.
!  Here J is computed by the user-supplied routine JAC if
!  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
!  If MITER = 3, a diagonal approximation to J is used.
!  J is stored in WM and replaced by P.  If MITER .ne. 3, P is then
!  subjected to LU decomposition in preparation for later solution
!  of linear systems with P as coefficient matrix.  This is done
!  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
!
!  In addition to variables described in DSTODE and DLSODE prologues,
!  communication with DPREPJ uses the following:
!  Y     = array containing predicted values on entry.
!  FTEM  = work array of length N (ACOR in DSTODE).
!  SAVF  = array containing f evaluated at predicted y.
!  WM    = real work space for matrices.  On output it contains the
!          inverse diagonal matrix if MITER = 3 and the LU decomposition
!          of P if MITER is 1, 2 , 4, or 5.
!          Storage of matrix elements starts at WM(3).
!          WM also contains the following matrix-related data:
!          WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
!          WM(2) = H*EL0, saved for later use if MITER = 3.
!  IWM   = integer work space containing pivot information, starting at
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!  EL0   = EL(1) (input).
!  IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
!          P matrix found to be singular.
!  JCUR  = output flag = 1 to indicate that the Jacobian matrix
!          (or approximation) is now current.
!  This routine also uses the COMMON variables EL0, H, TN, UROUND,
!  MITER, N, NFE, and NJE.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DGBFA, DGEFA, DVNORM
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890504  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DPREPJ
!**End
      EXTERNAL f, jac
      INTEGER neq, nyh, iwm
      DOUBLE PRECISION y, yh, ewt, ftem, savf, wm
      DIMENSION neq(*), y(*), yh(nyh, *), ewt(*), ftem(*), savf(*), wm(*), &
        iwm(*)
      INTEGER iownd, iowns, icf, ierpj, iersl, jcur, jstart, kflag, l, lyh, &
        lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      DOUBLE PRECISION rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      COMMON /dls001/rowns(209), ccmax, el0, h, hmin, hmxi, hu, rc, tn, &
        uround, iownd(6), iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, &
        lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, i1, i2, ier, ii, j, j1, jj, lenp, mba, mband, meb1, meband, &
        ml, ml3, mu, np1
      DOUBLE PRECISION con, di, fac, hl0, r, r0, srur, yi, yj, yjj !, dvnorm
!
!***FIRST EXECUTABLE STATEMENT  DPREPJ
      nje = nje + 1
      ierpj = 0
      jcur = 1
      hl0 = h*el0
!      GO TO (100, 200, 300, 400, 500), miter
      SELECT CASE(miter)
        CASE(1)
          GO TO 100
        CASE(2)
          GO TO 200
        CASE(3)
          GO TO 300
        CASE(4)
          GO TO 400
        CASE(5)
          GO TO 500
      END SELECT
! If MITER = 1, call JAC and multiply by scalar. -----------------------
100   lenp = n*n
      DO i = 1, lenp
        wm(i+2) = 0.0D0
      ENDDO
      CALL jac(neq, tn, y, 0, 0, wm(3), n)
      con = -hl0
      DO i = 1, lenp
        wm(i+2) = wm(i+2)*con
      ENDDO
      GO TO 240
! If MITER = 2, make N calls to F to approximate J. --------------------
200   fac = dvnorm(n, savf, ewt)
      r0 = 1000.0D0*abs(h)*uround*n*fac
      IF (r0==0.0D0) r0 = 1.0D0
      srur = wm(1)
      j1 = 2
      DO j = 1, n
        yj = y(j)
        r = max(srur*abs(yj), r0/ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        CALL f(neq, tn, y, ftem)
        DO i = 1, n
          wm(i+j1) = (ftem(i)-savf(i))*fac
        ENDDO
        y(j) = yj
        j1 = j1 + n
      ENDDO
      nfe = nfe + n
! Add identity matrix. -------------------------------------------------
240   j = 3
      np1 = n + 1
      DO i = 1, n
        wm(j) = wm(j) + 1.0D0
        j = j + np1
      ENDDO
! Do LU decomposition on P. --------------------------------------------
      CALL dgefa(wm(3), n, n, iwm(21), ier)
      IF (ier/=0) ierpj = 1
      RETURN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
300   wm(2) = hl0
      r = el0*0.1D0
      DO i = 1, n
        y(i) = y(i) + r*(h*savf(i)-yh(i,2))
      ENDDO
      CALL f(neq, tn, y, wm(3))
      nfe = nfe + 1
      DO i = 1, n
        r0 = h*savf(i) - yh(i, 2)
        di = 0.1D0*r0 - h*(wm(i+2)-savf(i))
        wm(i+2) = 1.0D0
        IF (abs(r0)<uround/ewt(i)) GO TO 320
        IF (abs(di)==0.0D0) GO TO 330
        wm(i+2) = 0.1D0*r0/di
320   ENDDO
      RETURN
330   ierpj = 1
      RETURN
! If MITER = 4, call JAC and multiply by scalar. -----------------------
400   ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      DO i = 1, lenp
        wm(i+2) = 0.0D0
      ENDDO
      CALL jac(neq, tn, y, ml, mu, wm(ml3), meband)
      con = -hl0
      DO i = 1, lenp
        wm(i+2) = wm(i+2)*con
      ENDDO
      GO TO 570
! If MITER = 5, make MBAND calls to F to approximate J. ----------------
500   ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min(mband, n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = dvnorm(n, savf, ewt)
      r0 = 1000.0D0*abs(h)*uround*n*fac
      IF (r0==0.0D0) r0 = 1.0D0
      DO j = 1, mba
        DO i = j, n, mband
          yi = y(i)
          r = max(srur*abs(yi), r0/ewt(i))
          y(i) = y(i) + r
        ENDDO
        CALL f(neq, tn, y, ftem)
        DO jj = j, n, mband
          y(jj) = yh(jj, 1)
          yjj = y(jj)
          r = max(srur*abs(yjj), r0/ewt(jj))
          fac = -hl0/r
          i1 = max(jj-mu, 1)
          i2 = min(jj+ml, n)
          ii = jj*meb1 - ml + 2
          DO i = i1, i2
            wm(ii+i) = (ftem(i)-savf(i))*fac
          ENDDO
        ENDDO
      ENDDO
      nfe = nfe + mba
! Add identity matrix. -------------------------------------------------
570   ii = mband + 2
      DO i = 1, n
        wm(ii) = wm(ii) + 1.0D0
        ii = ii + meband
      ENDDO
! Do LU decomposition of P. --------------------------------------------
      CALL dgbfa(wm(3), meband, n, ml, mu, iwm(21), ier)
      IF (ier/=0) ierpj = 1
      RETURN
!----------------------- END OF SUBROUTINE DPREPJ ----------------------
    END SUBROUTINE dprepj
!DECK DSOLSY
    SUBROUTINE dsolsy(wm, iwm, x, tem)
!***BEGIN PROLOGUE  DSOLSY
!***SUBSIDIARY
!***PURPOSE  ODEPACK linear system solver.
!***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This routine manages the solution of the linear system arising from
!  a chord iteration.  It is called if MITER .ne. 0.
!  If MITER is 1 or 2, it calls DGESL to accomplish this.
!  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
!  matrix, and then computes the solution.
!  If MITER is 4 or 5, it calls DGBSL.
!  Communication with DSOLSY uses the following variables:
!  WM    = real work space containing the inverse diagonal matrix if
!          MITER = 3 and the LU decomposition of the matrix otherwise.
!          Storage of matrix elements starts at WM(3).
!          WM also contains the following matrix-related data:
!          WM(1) = SQRT(UROUND) (not used here),
!          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
!  IWM   = integer work space containing pivot information, starting at
!          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
!          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
!  X     = the right-hand side vector on input, and the solution vector
!          on output, of length N.
!  TEM   = vector of work space of length N, not used in this version.
!  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
!          IERSL = 1 if a singular matrix arose with MITER = 3.
!  This routine also uses the COMMON variables EL0, H, MITER, and N.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DGBSL, DGESL
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSOLSY
!**End
      INTEGER iwm
      DOUBLE PRECISION wm, x, tem
      DIMENSION wm(*), iwm(*), x(*), tem(*)
      INTEGER iownd, iowns, icf, ierpj, iersl, jcur, jstart, kflag, l, lyh, &
        lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      DOUBLE PRECISION rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      COMMON /dls001/rowns(209), ccmax, el0, h, hmin, hmxi, hu, rc, tn, &
        uround, iownd(6), iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, &
        lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, maxcor, msbp, &
        mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, meband, ml, mu
      DOUBLE PRECISION di, hl0, phl0, r
!
!***FIRST EXECUTABLE STATEMENT  DSOLSY
      iersl = 0
!      GO TO (100, 100, 300, 400, 400), miter
      SELECT CASE(miter)
        CASE(1)
          GO TO 100
        CASE(2) 
          GO TO 100
        CASE(3)
          GO TO 300
        CASE(4)
          GO TO 400
        CASE(5)
          GO TO 400
      END SELECT
100   CALL dgesl(wm(3), n, n, iwm(21), x, 0)
      RETURN
!
300   phl0 = wm(2)
      hl0 = h*el0
      wm(2) = hl0
      IF (hl0==phl0) GO TO 330
      r = hl0/phl0
      DO i = 1, n
        di = 1.0D0 - r*(1.0D0-1.0D0/wm(i+2))
        IF (abs(di)==0.0D0) GO TO 390
        wm(i+2) = 1.0D0/di
      ENDDO
330   DO i = 1, n
        x(i) = wm(i+2)*x(i)
      ENDDO
      RETURN
390   iersl = 1
      RETURN
!
400   ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      CALL dgbsl(wm(3), meband, n, ml, mu, iwm(21), x, 0)
      RETURN
!----------------------- END OF SUBROUTINE DSOLSY ----------------------
    END SUBROUTINE dsolsy
!DECK DSRCOM
    SUBROUTINE dsrcom(rsav, isav, job)
!***BEGIN PROLOGUE  DSRCOM
!***SUBSIDIARY
!***PURPOSE  Save/restore ODEPACK COMMON blocks.
!***TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This routine saves or restores (depending on JOB) the contents of
!  the COMMON block DLS001, which is used internally
!  by one or more ODEPACK solvers.
!
!  RSAV = real array of length 218 or more.
!  ISAV = integer array of length 37 or more.
!  JOB  = flag indicating to save or restore the COMMON blocks:
!         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
!         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
!         A call with JOB = 2 presumes a prior call with JOB = 1.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   921116  Deleted treatment of block /EH0001/.  (ACH)
!   930801  Reduced Common block length by 2.  (ACH)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced Common block length by 209+12. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!   031112  Added SAVE statement for data-loaded constants.
!***END PROLOGUE  DSRCOM
!**End
      INTEGER isav, job
      INTEGER ils
      INTEGER i, lenils, lenrls
      DOUBLE PRECISION rsav, rls
      DIMENSION rsav(*), isav(*)
      SAVE lenrls, lenils
      COMMON /dls001/rls(218), ils(37)
      DATA lenrls/218/, lenils/37/
!
!***FIRST EXECUTABLE STATEMENT  DSRCOM
      IF (job==2) GO TO 100
!
      DO i = 1, lenrls
        rsav(i) = rls(i)
      ENDDO
      DO i = 1, lenils
        isav(i) = ils(i)
      ENDDO
      RETURN
!
100   CONTINUE
      DO i = 1, lenrls
        rls(i) = rsav(i)
      ENDDO
      DO i = 1, lenils
        ils(i) = isav(i)
      ENDDO
      RETURN
!----------------------- END OF SUBROUTINE DSRCOM ----------------------
    END SUBROUTINE dsrcom
!DECK DSTODE
    SUBROUTINE dstode(neq, y, yh, nyh, yh1, ewt, savf, acor, wm, iwm, f, jac, &
      pjac, slvs)
!***BEGIN PROLOGUE  DSTODE
!***SUBSIDIARY
!***PURPOSE  Performs one step of an ODEPACK integration.
!***TYPE      DOUBLE PRECISION (SSTODE-S, DSTODE-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  DSTODE performs one step of the integration of an initial value
!  problem for a system of ordinary differential equations.
!  Note:  DSTODE is independent of the value of the iteration method
!  indicator MITER, when this is .ne. 0, and hence is independent
!  of the type of chord method used, or the Jacobian structure.
!  Communication with DSTODE is done with the following variables:
!
!  NEQ    = integer array containing problem size in NEQ(1), and
!           passed as the NEQ argument in all calls to F and JAC.
!  Y      = an array of length .ge. N used as the Y argument in
!           all calls to F and JAC.
!  YH     = an NYH by LMAX array containing the dependent variables
!           and their approximate scaled derivatives, where
!           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!           j-th derivative of y(i), scaled by h**j/factorial(j)
!           (j = 0,1,...,NQ).  on entry for the first step, the first
!           two columns of YH must be set from the initial values.
!  NYH    = a constant integer .ge. N, the first dimension of YH.
!  YH1    = a one-dimensional array occupying the same space as YH.
!  EWT    = an array of length N containing multiplicative weights
!           for local error measurements.  Local errors in Y(i) are
!           compared to 1.0/EWT(i) in various error tests.
!  SAVF   = an array of working storage, of length N.
!           Also used for input of YH(*,MAXORD+2) when JSTART = -1
!           and MAXORD .lt. the current order NQ.
!  ACOR   = a work array of length N, used for the accumulated
!           corrections.  On a successful return, ACOR(i) contains
!           the estimated one-step local error in Y(i).
!  WM,IWM = real and integer work arrays associated with matrix
!           operations in chord iteration (MITER .ne. 0).
!  PJAC   = name of routine to evaluate and preprocess Jacobian matrix
!           and P = I - h*el0*JAC, if a chord method is being used.
!  SLVS   = name of routine to solve linear system in chord iteration.
!  CCMAX  = maximum relative change in h*el0 before PJAC is called.
!  H      = the step size to be attempted on the next step.
!           H is altered by the error control algorithm during the
!           problem.  H can be either positive or negative, but its
!           sign must remain constant throughout the problem.
!  HMIN   = the minimum absolute value of the step size h to be used.
!  HMXI   = inverse of the maximum absolute value of h to be used.
!           HMXI = 0.0 is allowed and corresponds to an infinite hmax.
!           HMIN and HMXI may be changed at any time, but will not
!           take effect until the next change of h is considered.
!  TN     = the independent variable. TN is updated on each step taken.
!  JSTART = an integer used for input only, with the following
!           values and meanings:
!                0  perform the first step.
!            .gt.0  take a new step continuing from the last.
!               -1  take the next step with a new value of H, MAXORD,
!                     N, METH, MITER, and/or matrix parameters.
!               -2  take the next step with a new value of H,
!                     but with other inputs unchanged.
!           On return, JSTART is set to 1 to facilitate continuation.
!  KFLAG  = a completion code with the following meanings:
!                0  the step was succesful.
!               -1  the requested error could not be achieved.
!               -2  corrector convergence could not be achieved.
!               -3  fatal error in PJAC or SLVS.
!           A return with KFLAG = -1 or -2 means either
!           abs(H) = HMIN or 10 consecutive failures occurred.
!           On a return with KFLAG negative, the values of TN and
!           the YH array are as of the beginning of the last
!           step, and H is the last step size attempted.
!  MAXORD = the maximum order of integration method to be allowed.
!  MAXCOR = the maximum number of corrector iterations allowed.
!  MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
!  MXNCF  = maximum number of convergence failures allowed.
!  METH/MITER = the method flags.  See description in driver.
!  N      = the number of first-order differential equations.
!  The values of CCMAX, H, HMIN, HMXI, TN, JSTART, KFLAG, MAXORD,
!  MAXCOR, MSBP, MXNCF, METH, MITER, and N are communicated via COMMON.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  DCFODE, DVNORM
!***COMMON BLOCKS    DLS001
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!   010418  Reduced size of Common block /DLS001/. (ACH)
!   031105  Restored 'own' variables to Common block /DLS001/, to
!           enable interrupt/restart feature. (ACH)
!***END PROLOGUE  DSTODE
!**End
      EXTERNAL f, jac, pjac, slvs
      INTEGER neq, nyh, iwm
      DOUBLE PRECISION y, yh, yh1, ewt, savf, acor, wm
      DIMENSION neq(*), y(*), yh(nyh, *), yh1(*), ewt(*), savf(*), acor(*), &
        wm(*), iwm(*)
      INTEGER iownd, ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj, iersl, &
        jcur, jstart, kflag, l, lyh, lewt, lacor, lsavf, lwm, liwm, meth, &
        miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      INTEGER i, i1, iredo, iret, j, jb, m, ncf, newq
      DOUBLE PRECISION conit, crate, el, elco, hold, rmax, tesco, ccmax, el0, &
        h, hmin, hmxi, hu, rc, tn, uround
      DOUBLE PRECISION dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, &
        rh, rhdn, rhsm, rhup, told! , dvnorm
      COMMON /dls001/conit, crate, el(13), elco(13, 12), hold, rmax, &
        tesco(3, 12), ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(6), &
        ialth, ipup, lmax, meo, nqnyh, nslp, icf, ierpj, iersl, jcur, jstart, &
        kflag, l, lyh, lewt, lacor, lsavf, lwm, liwm, meth, miter, maxord, &
        maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
!
!***FIRST EXECUTABLE STATEMENT  DSTODE
      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0D0
      IF (jstart>0) GO TO 200
      IF (jstart==-1) GO TO 100
      IF (jstart==-2) GO TO 160
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  RMAX is the maximum ratio by which H can be increased
! in a single step.  It is initially 1.E4 to compensate for the small
! initial H, but then is normally equal to 10.  If a failure
! occurs (in corrector convergence or error test), RMAX is set to 2
! for the next increase.
!-----------------------------------------------------------------------
      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0D0
      rc = 0.0D0
      el0 = 1.0D0
      crate = 0.7D0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      GO TO 140
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! IPUP is set to MITER to force a matrix update.
! If an order increase is about to be considered (IALTH = 1),
! IALTH is reset to 2 to postpone consideration one more step.
! If the caller has changed METH, DCFODE is called to reset
! the coefficients of the method.
! If the caller has changed MAXORD to a value less than the current
! order NQ, NQ is reduced to MAXORD, and a new H chosen accordingly.
! If H is to be changed, YH must be rescaled.
! If H or METH is being changed, IALTH is reset to L = NQ + 1
! to prevent further changes in H for that many steps.
!-----------------------------------------------------------------------
100   ipup = miter
      lmax = maxord + 1
      IF (ialth==1) ialth = 2
      IF (meth==meo) GO TO 110
      CALL dcfode(meth, elco, tesco)
      meo = meth
      IF (nq>maxord) GO TO 120
      ialth = l
      iret = 1
      GO TO 150
110   IF (nq<=maxord) GO TO 160
120   nq = maxord
      l = lmax
      DO i = 1, l
        el(i) = elco(i, nq)
      ENDDO
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5D0/(nq+2)
      ddn = dvnorm(n, savf, ewt)/tesco(1, l)
      exdn = 1.0D0/l
      rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
      rh = min(rhdn, 1.0D0)
      iredo = 3
      IF (h==hold) GO TO 170
      rh = min(rh, abs(h/hold))
      h = hold
      GO TO 175
!-----------------------------------------------------------------------
! DCFODE is called to get all the integration coefficients for the
! current METH.  Then the EL vector and related constants are reset
! whenever the order NQ is changed, or at the start of the problem.
!-----------------------------------------------------------------------
140   CALL dcfode(meth, elco, tesco)
150   DO i = 1, l
        el(i) = elco(i, nq)
      ENDDO
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5D0/(nq+2)
!      GO TO (160, 170, 200), iret
      SELECT CASE(iret)
        CASE(1)
          GO TO 160
        CASE(2)
          GO TO 170
        CASE(3)
          GO TO 200
      END SELECT
!-----------------------------------------------------------------------
! If H is being changed, the H ratio RH is checked against
! RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
! L = NQ + 1 to prevent a change of H for that many steps, unless
! forced by a convergence or error test failure.
!-----------------------------------------------------------------------
160   IF (h==hold) GO TO 200
      rh = h/hold
      h = hold
      iredo = 3
      GO TO 175
170   rh = max(rh, hmin/abs(h))
175   rh = min(rh, rmax)
      rh = rh/max(1.0D0, abs(h)*hmxi*rh)
      r = 1.0D0
      DO j = 2, l
        r = r*rh
        DO i = 1, n
          yh(i, j) = yh(i, j)*r
        ENDDO
      ENDDO
      h = h*rh
      rc = rc*rh
      ialth = l
      IF (iredo==0) GO TO 690
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal Triangle matrix.
! RC is the ratio of new to old values of the coefficient  H*EL(1).
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force PJAC to be called, if a Jacobian is involved.
! In any case, PJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
200   IF (abs(rc-1.0D0)>ccmax) ipup = miter
      IF (nst>=nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      DO jb = 1, nq
        i1 = i1 - nyh
!dir$ ivdep
        DO i = i1, nqnyh
          yh1(i) = yh1(i) + yh1(i+nyh)
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the R.M.S. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220   m = 0
      DO i = 1, n
        y(i) = yh(i, 1)
      ENDDO
      CALL f(neq, tn, y, savf)
      nfe = nfe + 1
      IF (ipup<=0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*el(1)*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      CALL pjac(neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac)
      ipup = 0
      rc = 1.0D0
      nslp = nst
      crate = 0.7D0
      IF (ierpj/=0) GO TO 430
250   DO i = 1, n
        acor(i) = 0.0D0
      ENDDO
270   IF (miter/=0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO i = 1, n
        savf(i) = h*savf(i) - yh(i, 2)
        y(i) = savf(i) - acor(i)
      ENDDO
      del = dvnorm(n, y, ewt)
      DO i = 1, n
        y(i) = yh(i, 1) + el(1)*savf(i)
        acor(i) = savf(i)
      ENDDO
      GO TO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.
!-----------------------------------------------------------------------
350   DO i = 1, n
        y(i) = h*savf(i) - (yh(i,2)+acor(i))
      ENDDO
      CALL slvs(wm, iwm, y, savf)
      IF (iersl<0) GO TO 430
      IF (iersl>0) GO TO 410
      del = dvnorm(n, y, ewt)
      DO i = 1, n
        acor(i) = acor(i) + y(i)
        y(i) = yh(i, 1) + el(1)*acor(i)
      ENDDO
!-----------------------------------------------------------------------
! Test for convergence.  If M.gt.0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
400   IF (m/=0) crate = max(0.2D0*crate, del/delp)
      dcon = del*min(1.0D0, 1.5D0*crate)/(tesco(2,nq)*conit)
      IF (dcon<=1.0D0) GO TO 450
      m = m + 1
      IF (m==maxcor) GO TO 410
      IF (m>=2 .AND. del>2.0D0*delp) GO TO 410
      delp = del
      CALL f(neq, tn, y, savf)
      nfe = nfe + 1
      GO TO 270
!-----------------------------------------------------------------------
! The corrector iteration failed to converge.
! If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
! the next try.  Otherwise the YH array is retracted to its values
! before prediction, and H is reduced, if possible.  If H cannot be
! reduced or MXNCF failures have occurred, exit with KFLAG = -2.
!-----------------------------------------------------------------------
410   IF (miter==0 .OR. jcur==1) GO TO 430
      icf = 1
      ipup = miter
      GO TO 220
430   icf = 2
      ncf = ncf + 1
      rmax = 2.0D0
      tn = told
      i1 = nqnyh + 1
      DO jb = 1, nq
        i1 = i1 - nyh
!dir$ ivdep
        DO i = i1, nqnyh
          yh1(i) = yh1(i) - yh1(i+nyh)
        ENDDO
      ENDDO
      IF (ierpj<0 .OR. iersl<0) GO TO 680
      IF (abs(h)<=hmin*1.00001D0) GO TO 670
      IF (ncf==mxncf) GO TO 670
      rh = 0.25D0
      ipup = miter
      iredo = 1
      GO TO 170
!-----------------------------------------------------------------------
! The corrector has converged.  JCUR is set to 0
! to signal that the Jacobian involved may need updating later.
! The local error test is made and control passes to statement 500
! if it fails.
!-----------------------------------------------------------------------
450   jcur = 0
      IF (m==0) dsm = del/tesco(2, nq)
      IF (m>0) dsm = dvnorm(n, acor, ewt)/tesco(2, nq)
      IF (dsm>1.0D0) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH array.
! Consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
! If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
! use in a possible order increase on the next step.
! If a change in H is considered, an increase or decrease in order
! by one is considered also.  A change in H is made only if it is by a
! factor of at least 1.1.  If not, IALTH is set to 3 to prevent
! testing for that many steps.
!-----------------------------------------------------------------------
      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      DO j = 1, l
        DO i = 1, n
          yh(i, j) = yh(i, j) + el(j)*acor(i)
        ENDDO
      ENDDO
      ialth = ialth - 1
      IF (ialth==0) GO TO 520
      IF (ialth>1) GO TO 700
      IF (l==lmax) GO TO 700
      DO i = 1, n
        yh(i, lmax) = acor(i)
      ENDDO
      GO TO 700
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for this or
! one lower order.  After 2 or more failures, H is forced to decrease
! by a factor of 0.2 or less.
!-----------------------------------------------------------------------
500   kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      DO jb = 1, nq
        i1 = i1 - nyh
!dir$ ivdep
        DO i = i1, nqnyh
          yh1(i) = yh1(i) - yh1(i+nyh)
        ENDDO
      ENDDO
      rmax = 2.0D0
      IF (abs(h)<=hmin*1.00001D0) GO TO 660
      IF (kflag<=-3) GO TO 640
      iredo = 2
      rhup = 0.0D0
      GO TO 540
!-----------------------------------------------------------------------
! Regardless of the success or failure of the step, factors
! RHDN, RHSM, and RHUP are computed, by which H could be multiplied
! at order NQ - 1, order NQ, or order NQ + 1, respectively.
! In the case of failure, RHUP = 0.0 to avoid an order increase.
! The largest of these is determined and the new order chosen
! accordingly.  If the order is to be increased, we compute one
! additional scaled derivative.
!-----------------------------------------------------------------------
520   rhup = 0.0D0
      IF (l==lmax) GO TO 540
      DO i = 1, n
        savf(i) = acor(i) - yh(i, lmax)
      ENDDO
      dup = dvnorm(n, savf, ewt)/tesco(3, nq)
      exup = 1.0D0/(l+1)
      rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
540   exsm = 1.0D0/l
      rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
      rhdn = 0.0D0
      IF (nq==1) GO TO 560
      ddn = dvnorm(n, yh(1,l), ewt)/tesco(1, nq)
      exdn = 1.0D0/nq
      rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
560   IF (rhsm>=rhup) GO TO 570
      IF (rhup>rhdn) GO TO 590
      GO TO 580
570   IF (rhsm<rhdn) GO TO 580
      newq = nq
      rh = rhsm
      GO TO 620
580   newq = nq - 1
      rh = rhdn
      IF (kflag<0 .AND. rh>1.0D0) rh = 1.0D0
      GO TO 620
590   newq = l
      rh = rhup
      IF (rh<1.1D0) GO TO 610
      r = el(l)/l
      DO i = 1, n
        yh(i, newq+1) = acor(i)*r
      ENDDO
      GO TO 630
610   ialth = 3
      GO TO 700
620   IF ((kflag==0) .AND. (rh<1.1D0)) GO TO 610
      IF (kflag<=-2) rh = min(rh, 0.2D0)
!-----------------------------------------------------------------------
! If there is a change of order, reset NQ, l, and the coefficients.
! In any case H is reset according to RH and the YH array is rescaled.
! Then exit from 690 if the step was OK, or redo the step otherwise.
!-----------------------------------------------------------------------
      IF (newq==nq) GO TO 170
630   nq = newq
      l = nq + 1
      iret = 2
      GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more failures have occured.
! If 10 failures have occurred, exit with KFLAG = -1.
! It is assumed that the derivatives that have accumulated in the
! YH array have errors of the wrong order.  Hence the first
! derivative is recomputed, and the order is set to 1.  Then
! H is reduced by a factor of 10, and the step is retried,
! until it succeeds or H reaches HMIN.
!-----------------------------------------------------------------------
640   IF (kflag==-10) GO TO 660
      rh = 0.1D0
      rh = max(hmin/abs(h), rh)
      h = h*rh
      DO i = 1, n
        y(i) = yh(i, 1)
      ENDDO
      CALL f(neq, tn, y, savf)
      nfe = nfe + 1
      DO i = 1, n
        yh(i, 2) = h*savf(i)
      ENDDO
      ipup = miter
      ialth = 5
      IF (nq==1) GO TO 200
      nq = 1
      l = 2
      iret = 3
      GO TO 150
!-----------------------------------------------------------------------
! All returns are made through this section.  H is saved in HOLD
! to allow the caller to change H on the next step.
!-----------------------------------------------------------------------
660   kflag = -1
      GO TO 720
670   kflag = -2
      GO TO 720
680   kflag = -3
      GO TO 720
690   rmax = 10.0D0
700   r = 1.0D0/tesco(2, nqu)
      DO i = 1, n
        acor(i) = acor(i)*r
      ENDDO
720   hold = h
      jstart = 1
      RETURN
!----------------------- END OF SUBROUTINE DSTODE ----------------------
    END SUBROUTINE dstode
!DECK DEWSET
    SUBROUTINE dewset(n, itol, rtol, atol, ycur, ewt)
!***BEGIN PROLOGUE  DEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DEWSET
!**End
      INTEGER n, itol
      INTEGER i
      DOUBLE PRECISION rtol, atol, ycur, ewt
      DIMENSION rtol(*), atol(*), ycur(n), ewt(n)
!
!***FIRST EXECUTABLE STATEMENT  DEWSET
!      GO TO (10, 20, 30, 40), itol
      SELECT CASE(itol)
        CASE(1)
          GO TO 10
        CASE(2)
          GO TO 20
        CASE(3)
          GO TO 30
        CASE(4)
          GO TO 40
      END SELECT
10    CONTINUE
      DO i = 1, n
        ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      ENDDO
      RETURN
20    CONTINUE
      DO i = 1, n
        ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      ENDDO
      RETURN
30    CONTINUE
      DO i = 1, n
        ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      ENDDO
      RETURN
40    CONTINUE
      DO i = 1, n
        ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      ENDDO
      RETURN
!----------------------- END OF SUBROUTINE DEWSET ----------------------
    END SUBROUTINE dewset
!DECK DVNORM
    DOUBLE PRECISION FUNCTION dvnorm(n, v, w)
!***BEGIN PROLOGUE  DVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the array V, with weights
!  contained in the array W of length N:
!    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DVNORM
!**End
      INTEGER n, i
      DOUBLE PRECISION v, w, sum
      DIMENSION v(n), w(n)
!
!***FIRST EXECUTABLE STATEMENT  DVNORM
      sum = 0.0D0
      DO i = 1, n
        sum = sum + (v(i)*w(i))**2
      ENDDO
      dvnorm = sqrt(sum/n)
      RETURN
!----------------------- END OF FUNCTION DVNORM ------------------------
    END FUNCTION dvnorm

!DECK DGEFA
    SUBROUTINE dgefa(a, lda, n, ipvt, info)
!***BEGIN PROLOGUE  DGEFA
!***PURPOSE  Factor a matrix using Gaussian elimination.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
!***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGEFA factors a double precision matrix by Gaussian elimination.
!
!     DGEFA is usually called by DGECO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be factored.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL or DGEDI will divide by zero
!                     if called.  Use  RCOND  in DGECO for a reliable
!                     indication of singularity.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGEFA
      INTEGER lda, n, ipvt(*), info
      DOUBLE PRECISION a(lda, *)
!
      DOUBLE PRECISION t
      INTEGER idamax, j, k, kp1, l, nm1
!
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
!***FIRST EXECUTABLE STATEMENT  DGEFA
      info = 0
      nm1 = n - 1
      IF (nm1<1) GO TO 70
      DO k = 1, nm1
        kp1 = k + 1
!
!        FIND L = PIVOT INDEX
!
        l = idamax(n-k+1, a(k,k), 1) + k - 1
        ipvt(k) = l
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
        IF (a(l,k)==0.0D0) GO TO 40
!
!           INTERCHANGE IF NECESSARY
!
        IF (l==k) GO TO 10
        t = a(l, k)
        a(l, k) = a(k, k)
        a(k, k) = t
10      CONTINUE
!
!           COMPUTE MULTIPLIERS
!
        t = -1.0D0/a(k, k)
        CALL dscal(n-k, t, a(k+1,k), 1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
        DO j = kp1, n
          t = a(l, j)
          IF (l==k) GO TO 20
          a(l, j) = a(k, j)
          a(k, j) = t
20        CONTINUE
          CALL daxpy(n-k, t, a(k+1,k), 1, a(k+1,j), 1)
        ENDDO
        GO TO 50
40      CONTINUE
        info = k
50      CONTINUE
      ENDDO
70    CONTINUE
      ipvt(n) = n
      IF (a(n,n)==0.0D0) info = n
      RETURN
    END SUBROUTINE dgefa
!DECK DGESL
    SUBROUTINE dgesl(a, lda, n, ipvt, b, job)
!***BEGIN PROLOGUE  DGESL
!***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
!            factors computed by DGECO or DGEFA.
!***CATEGORY  D2A1
!***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGESL solves the double precision system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B  where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND .GT. 0.0
!        or DGEFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGESL
      INTEGER lda, n, ipvt(*), job
      DOUBLE PRECISION a(lda, *), b(*)
!
      DOUBLE PRECISION ddot, t
      INTEGER k, kb, l, nm1
!***FIRST EXECUTABLE STATEMENT  DGESL
      nm1 = n - 1
      IF (job/=0) GO TO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE  L*Y = B
!
      IF (nm1<1) GO TO 30
      DO k = 1, nm1
        l = ipvt(k)
        t = b(l)
        IF (l==k) GO TO 10
        b(l) = b(k)
        b(k) = t
10      CONTINUE
        CALL daxpy(n-k, t, a(k+1,k), 1, b(k+1), 1)
      ENDDO
30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
      DO kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/a(k, k)
        t = -b(k)
        CALL daxpy(k-1, t, a(1,k), 1, b(1), 1)
      ENDDO
      GO TO 100
50    CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
      DO k = 1, n
        t = ddot(k-1, a(1,k), 1, b(1), 1)
        b(k) = (b(k)-t)/a(k, k)
      ENDDO
!
!        NOW SOLVE TRANS(L)*X = Y
!
      IF (nm1<1) GO TO 90
      DO kb = 1, nm1
        k = n - kb
        b(k) = b(k) + ddot(n-k, a(k+1,k), 1, b(k+1), 1)
        l = ipvt(k)
        IF (l==k) GO TO 70
        t = b(l)
        b(l) = b(k)
        b(k) = t
70      CONTINUE
      ENDDO
90    CONTINUE
100   CONTINUE
      RETURN
    END SUBROUTINE dgesl
!DECK DGBFA
    SUBROUTINE dgbfa(abd, lda, n, ml, mu, ipvt, info)
!***BEGIN PROLOGUE  DGBFA
!***PURPOSE  Factor a band matrix using Gaussian elimination.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBFA factors a double precision band matrix by elimination.
!
!     DGBFA is usually called by DGBCO, but it can be called
!     directly with a saving in time if  RCOND  is not needed.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                contains the matrix in band storage.  The columns
!                of the matrix are stored in the columns of  ABD  and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of  ABD .
!                See the comments below for details.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!                LDA must be .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 .LE. ML .LT.  N .
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 .LE. MU .LT.  N .
!                More efficient if  ML .LE. MU .
!     On Return
!
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written  A = L*U  where
!                L  is a product of permutation and unit lower
!                triangular matrices and  U  is upper triangular.
!
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL will divide by zero if
!                     called.  Use  RCOND  in DGBCO for a reliable
!                     indication of singularity.
!
!     Band Storage
!
!           If  A  is a band matrix, the following program segment
!           will set up the input.
!
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
!           In addition, the first  ML  rows in  ABD  are used for
!           elements generated during the triangularization.
!           The total number of rows needed in  ABD  is  2*ML+MU+1 .
!           The  ML+MU by ML+MU  upper left triangle and the
!           ML by ML  lower right triangle are not referenced.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBFA
      INTEGER lda, n, ml, mu, ipvt(*), info
      DOUBLE PRECISION abd(lda, *)
!
      DOUBLE PRECISION t
      INTEGER i, idamax, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1
!
!***FIRST EXECUTABLE STATEMENT  DGBFA
      m = ml + mu + 1
      info = 0
!
!     ZERO INITIAL FILL-IN COLUMNS
!
      j0 = mu + 2
      j1 = min(n, m) - 1
      IF (j1<j0) GO TO 30
      DO jz = j0, j1
        i0 = m + 1 - jz
        DO i = i0, ml
          abd(i, jz) = 0.0D0
        ENDDO
      ENDDO
30    CONTINUE
      jz = j1
      ju = 0
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      nm1 = n - 1
      IF (nm1<1) GO TO 130
      DO k = 1, nm1
        kp1 = k + 1
!
!        ZERO NEXT FILL-IN COLUMN
!
        jz = jz + 1
        IF (jz>n) GO TO 50
        IF (ml<1) GO TO 50
        DO i = 1, ml
          abd(i, jz) = 0.0D0
        ENDDO
50      CONTINUE
!
!        FIND L = PIVOT INDEX
!
        lm = min(ml, n-k)
        l = idamax(lm+1, abd(m,k), 1) + m - 1
        ipvt(k) = l + k - m
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
        IF (abd(l,k)==0.0D0) GO TO 100
!
!           INTERCHANGE IF NECESSARY
!
        IF (l==m) GO TO 60
        t = abd(l, k)
        abd(l, k) = abd(m, k)
        abd(m, k) = t
60      CONTINUE
!
!           COMPUTE MULTIPLIERS
!
        t = -1.0D0/abd(m, k)
        CALL dscal(lm, t, abd(m+1,k), 1)
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
        ju = min(max(ju,mu+ipvt(k)), n)
        mm = m
        IF (ju<kp1) GO TO 90
        DO j = kp1, ju
          l = l - 1
          mm = mm - 1
          t = abd(l, j)
          IF (l==mm) GO TO 70
          abd(l, j) = abd(mm, j)
          abd(mm, j) = t
70        CONTINUE
          CALL daxpy(lm, t, abd(m+1,k), 1, abd(mm+1,j), 1)
        ENDDO
90      CONTINUE
        GO TO 110
100     CONTINUE
        info = k
110     CONTINUE
      ENDDO
130   CONTINUE
      ipvt(n) = n
      IF (abd(m,n)==0.0D0) info = n
      RETURN
    END SUBROUTINE dgbfa
!DECK DGBSL
    SUBROUTINE dgbsl(abd, lda, n, ml, mu, ipvt, b, job)
!***BEGIN PROLOGUE  DGBSL
!***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
!            the factors computed by DGBCO or DGBFA.
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGBSL solves the double precision band system
!     A * X = B  or  TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DGBCO or DGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!        JOB     INTEGER
!                = 0         to solve  A*X = B ,
!                = nonzero   to solve  TRANS(A)*X = B , where
!                            TRANS(A)  is the transpose.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal.  Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA .  It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND .GT. 0.0
!        or DGBFA has set INFO .EQ. 0 .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGBSL
      INTEGER lda, n, ml, mu, ipvt(*), job
      DOUBLE PRECISION abd(lda, *), b(*)
!
      DOUBLE PRECISION ddot, t
      INTEGER k, kb, l, la, lb, lm, m, nm1
!***FIRST EXECUTABLE STATEMENT  DGBSL
      m = mu + ml + 1
      nm1 = n - 1
      IF (job/=0) GO TO 50
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
      IF (ml==0) GO TO 30
      IF (nm1<1) GO TO 30
      DO k = 1, nm1
        lm = min(ml, n-k)
        l = ipvt(k)
        t = b(l)
        IF (l==k) GO TO 10
        b(l) = b(k)
        b(k) = t
10      CONTINUE
        CALL daxpy(lm, t, abd(m+1,k), 1, b(k+1), 1)
      ENDDO
30    CONTINUE
!
!        NOW SOLVE  U*X = Y
!
      DO kb = 1, n
        k = n + 1 - kb
        b(k) = b(k)/abd(m, k)
        lm = min(k, m) - 1
        la = m - lm
        lb = k - lm
        t = -b(k)
        CALL daxpy(lm, t, abd(la,k), 1, b(lb), 1)
      ENDDO
      GO TO 100
50    CONTINUE
!
!        JOB = NONZERO, SOLVE  TRANS(A) * X = B
!        FIRST SOLVE  TRANS(U)*Y = B
!
      DO k = 1, n
        lm = min(k, m) - 1
        la = m - lm
        lb = k - lm
        t = ddot(lm, abd(la,k), 1, b(lb), 1)
        b(k) = (b(k)-t)/abd(m, k)
      ENDDO
!
!        NOW SOLVE TRANS(L)*X = Y
!
      IF (ml==0) GO TO 90
      IF (nm1<1) GO TO 90
      DO kb = 1, nm1
        k = n - kb
        lm = min(ml, n-k)
        b(k) = b(k) + ddot(lm, abd(m+1,k), 1, b(k+1), 1)
        l = ipvt(k)
        IF (l==k) GO TO 70
        t = b(l)
        b(l) = b(k)
        b(k) = t
70      CONTINUE
      ENDDO
90    CONTINUE
100   CONTINUE
      RETURN
    END SUBROUTINE dgbsl
!DECK DAXPY
    SUBROUTINE daxpy(n, da, dx, incx, dy, incy)
!***BEGIN PROLOGUE  DAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)
!
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DAXPY
      DOUBLE PRECISION dx(*), dy(*), da
      INTEGER n, incx, incy, ix,iy,i,m,mp1,ns
!***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (n<=0 .OR. da==0.0D0) RETURN
!      IF (incx==incy) IF (incx-1) 5, 20, 60
      IF (incx==incy) THEN
        IF (incx-1 < 0) GO TO 5
        IF (incx-1 == 0) GO TO 20
        IF (incx-1 > 0 ) GO TO 60
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
5     ix = 1
      iy = 1
      IF (incx<0) ix = (-n+1)*incx + 1
      IF (incy<0) iy = (-n+1)*incy + 1
      DO i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      ENDDO
      RETURN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 4.
!
20    m = mod(n, 4)
      IF (m==0) GO TO 40
      DO i = 1, m
        dy(i) = dy(i) + da*dx(i)
      END DO
      IF (n<4) RETURN
40    mp1 = m + 1
      DO i = mp1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i+1) = dy(i+1) + da*dx(i+1)
        dy(i+2) = dy(i+2) + da*dx(i+2)
        dy(i+3) = dy(i+3) + da*dx(i+3)
      ENDDO
      RETURN
!
!     Code for equal, positive, non-unit increments.
!
60    ns = n*incx
      DO i = 1, ns, incx
        dy(i) = da*dx(i) + dy(i)
      ENDDO
      RETURN
    END SUBROUTINE daxpy
!DECK DDOT
    DOUBLE PRECISION FUNCTION ddot(n, dx, incx, dy, incy)
!***BEGIN PROLOGUE  DDOT
!***PURPOSE  Compute the inner product of two vectors.
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!     DDOT  double precision dot product (zero if N .LE. 0)
!
!     Returns the dot product of double precision DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DDOT
      DOUBLE PRECISION dx(*), dy(*)
      INTEGER n,incx,incy,ix,iy,i,m,mp1,ns
!***FIRST EXECUTABLE STATEMENT  DDOT
      ddot = 0.0D0
      IF (n<=0) RETURN
!      IF (incx==incy) IF (incx-1) 5, 20, 60
      IF (incx==incy) THEN
        IF (incx-1 < 0) GO TO 5
        IF (incx-1 == 0) GO TO 20
        IF (incx-1 > 0) GO TO 60
      ENDIF
!
!     Code for unequal or nonpositive increments.
!
5     ix = 1
      iy = 1
      IF (incx<0) ix = (-n+1)*incx + 1
      IF (incy<0) iy = (-n+1)*incy + 1
      DO i = 1, n
        ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      ENDDO
      RETURN
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
20    m = mod(n, 5)
      IF (m==0) GO TO 40
      DO i = 1, m
        ddot = ddot + dx(i)*dy(i)
      ENDDO
      IF (n<5) RETURN
40    mp1 = m + 1
      DO i = mp1, n, 5
        ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + &
          dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
      ENDDO
      RETURN
!
!     Code for equal, positive, non-unit increments.
!
60    ns = n*incx
      DO i = 1, ns, incx
        ddot = ddot + dx(i)*dy(i)
      ENDDO
      RETURN
    END FUNCTION ddot
!DECK DSCAL
    SUBROUTINE dscal(n, da, dx, incx)
!***BEGIN PROLOGUE  DSCAL
!***PURPOSE  Multiply a vector by a constant.
!***CATEGORY  D1A6
!***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scale factor
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!       DX  double precision result (unchanged if N.LE.0)
!
!     Replace double precision DX by double precision DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSCAL
      DOUBLE PRECISION da, dx(*)
      INTEGER i, incx, ix, m, mp1, n
!***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (n<=0) RETURN
      IF (incx==1) GO TO 20
!
!     Code for increment not equal to 1.
!
      ix = 1
      IF (incx<0) ix = (-n+1)*incx + 1
      DO i = 1, n
        dx(ix) = da*dx(ix)
        ix = ix + incx
      ENDDO
      RETURN
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
20    m = mod(n, 5)
      IF (m==0) GO TO 40
      DO i = 1, m
        dx(i) = da*dx(i)
      ENDDO
      IF (n<5) RETURN
40    mp1 = m + 1
      DO i = mp1, n, 5
        dx(i) = da*dx(i)
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
      ENDDO
      RETURN
    END SUBROUTINE dscal
!DECK IDAMAX
    INTEGER FUNCTION idamax(n, dx, incx)
!***BEGIN PROLOGUE  IDAMAX
!***PURPOSE  Find the smallest index of that component of a vector
!            having the maximum magnitude.
!***CATEGORY  D1A2
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N .LE. 0)
!
!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IDAMAX
      DOUBLE PRECISION dx(*), dmax, xmag
      INTEGER i, incx, ix, n
!***FIRST EXECUTABLE STATEMENT  IDAMAX
      idamax = 0
      IF (n<=0) RETURN
      idamax = 1
      IF (n==1) RETURN
!
      IF (incx==1) GO TO 20
!
!     Code for increments not equal to 1.
!
      ix = 1
      IF (incx<0) ix = (-n+1)*incx + 1
      dmax = abs(dx(ix))
      ix = ix + incx
      DO i = 2, n
        xmag = abs(dx(ix))
        IF (xmag>dmax) THEN
          idamax = i
          dmax = xmag
        END IF
        ix = ix + incx
      ENDDO
      RETURN
!
!     Code for increments equal to 1.
!
20    dmax = abs(dx(1))
      DO i = 2, n
        xmag = abs(dx(i))
        IF (xmag>dmax) THEN
          idamax = i
          dmax = xmag
        END IF
      ENDDO
      RETURN
    END FUNCTION idamax
!DECK XERRWD
    SUBROUTINE xerrwd(msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.
!
      DOUBLE PRECISION r1, r2
      INTEGER nmes, nerr, level, ni, i1, i2, nr
!      CHARACTER *(*) msg
      CHARACTER(LEN=*) msg
!
!  Declare local variables.
!
      INTEGER lunit,mesflg!, ixsav, mesflg
!
!  Get logical unit number and message print flag.
!
!***FIRST EXECUTABLE STATEMENT  XERRWD
      lunit = ixsav(1, 0, .FALSE.)
      mesflg = ixsav(2, 0, .FALSE.)
      IF (mesflg==0) GO TO 100
!
!  Write the message.
!
      WRITE (lunit, 10) msg
10    FORMAT (1X, A)
      IF (ni==1) WRITE (lunit, 20) i1
20    FORMAT (6X, 'In above message,  I1 =', I10)
      IF (ni==2) WRITE (lunit, 30) i1, i2
30    FORMAT (6X, 'In above message,  I1 =', I10, 3X, 'I2 =', I10)
      IF (nr==1) WRITE (lunit, 40) r1
40    FORMAT (6X, 'In above message,  R1 =', D21.13)
      IF (nr==2) WRITE (lunit, 50) r1, r2
50    FORMAT (6X, 'In above,  R1 =', D21.13, 3X, 'R2 =', D21.13)
!
!  Abort the run if LEVEL = 2.
!
100   IF (level/=2) RETURN
      STOP
!----------------------- End of Subroutine XERRWD ----------------------
    END SUBROUTINE xerrwd
!DECK IXSAV
    INTEGER FUNCTION ixsav(ipar, ivalue, iset)
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!    LUNIT, the logical unit number to which messages are printed, and
!    MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!   LUNIT  = Logical unit number for messages.  The default is obtained
!            by a call to IUMACH (may be machine-dependent).
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  On input..
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!    ISET   = Logical flag to indicate whether to read or write.
!             If ISET = .TRUE., the parameter will be given
!             the value IVALUE.  If ISET = .FALSE., the parameter
!             will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!    IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Modified prologue to SLATEC format. (FNF)
!   930915  Added IUMACH call to get default output unit.  (ACH)
!   930922  Minor cosmetic changes. (FNF)
!   010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
      LOGICAL iset
      INTEGER ipar, ivalue
!-----------------------------------------------------------------------
      INTEGER lunit,mesflg!iumach, lunit, mesflg
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      SAVE lunit, mesflg
      DATA lunit/ -1/, mesflg/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (ipar==1) THEN
        IF (lunit==-1) lunit = iumach()
        ixsav = lunit
        IF (iset) lunit = ivalue
      END IF
!
      IF (ipar==2) THEN
        ixsav = mesflg
        IF (iset) mesflg = ivalue
      END IF
!
      RETURN
!----------------------- End of Function IXSAV -------------------------
    END FUNCTION ixsav
!DECK IUMACH
    INTEGER FUNCTION iumach()
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
      iumach = 6
!
      RETURN
!----------------------- End of Function IUMACH ------------------------
    END FUNCTION iumach
END MODULE mo_fplume_dlsode
