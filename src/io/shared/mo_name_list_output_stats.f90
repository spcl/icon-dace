! Module generating a PostScript file visualizing the asynchronous I/O.
!
! The output write intervals are visualized as bars on a timeline,
! each processor on its separate row. The time line is given in
! seconds of wall-clock time. The interval bars are annotated by the
! corresponding filenames. Where there is no sufficient space left
! (on the paper), the annotation is ommitted.
!
! The start and end of each interval is triggered by user
! calls. These calls must happen consecutively (START1, END1, START2,
! END2, ...) but are otherwise not checked for consistency.
!
! Note that this Fortran module is restricted to high level calls
! while the actual drawing and scaling is automatically handled by
! the generated PostScript code.
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

MODULE mo_name_list_output_stats

  USE mo_kind,                      ONLY: wp
  USE mo_exception,                 ONLY: finish
  USE mo_impl_constants,            ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_io_units,                  ONLY: find_next_free_unit, &
    &                                     FILENAME_MAX
  USE mo_mpi,                       ONLY: p_gather, p_send, p_recv, &
    &                                     p_comm_rank, p_comm_size

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_reference_time
  PUBLIC :: interval_start
  PUBLIC :: interval_end
  PUBLIC :: interval_write_psfile

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_stats'

  !> No. of additional entries to allocate when resizing interval list
  INTEGER, PARAMETER :: ALLOCATE_AMOUNT = 1000

  !> Type definition: A single interval
  TYPE t_interval
    CHARACTER(LEN=FILENAME_MAX) :: annotation          !< annotation string for this interval
    REAL(wp)                    :: sec_start, sec_end  !< start, end of this interval
  END TYPE t_interval

  !> Type definition: List of intervals
  TYPE t_list
    INTEGER                       :: nentries          !< no. of intervals in list (may be smaller than "SIZE(intvl)")
    INTEGER                       :: current_idx       !< index of open interval
    TYPE(t_interval), ALLOCATABLE :: intvl(:)          !< list of intervals
  END TYPE t_list

  !> Type definition for a "list of lists": several lists collected
  !> from different PEs
  TYPE t_global_list
    TYPE(t_list)                                 :: list         !< several lists collected from different PEs
    CHARACTER(LEN=MAX_CHAR_LENGTH), ALLOCATABLE  :: pe_names(:)  !< name string for each PE
    INTEGER,                        ALLOCATABLE  :: start_idx(:) !< start index for each PEs intervals in total list
  END TYPE t_global_list

  !> Type definition: time stamp
  TYPE t_timestamp
    ! Julian day ; note that we ignore the UTC shift when querying
    ! time stamps since were are eventually calculating differences.
    REAL(wp) :: julian_day
  END TYPE t_timestamp

  !> List of intervals
  TYPE (t_list) :: list

  !> Initial time stamp used as reference (0 seconds wall-time)
  TYPE (t_timestamp) :: start_time = t_timestamp(-1._wp)

CONTAINS

  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript file header.
  !
  SUBROUTINE ps_define_header(psfile)
    INTEGER, INTENT(IN) :: psfile   !< PostScript file handle

    WRITE (psfile, '(a)') "%!PS"
    WRITE (psfile, '(a)') "% ICON Output Timings - Auto-Generated PostScript Plot"
    WRITE (psfile, '(a)') ""
    WRITE (psfile, '(a)') "%%BeginPageSetup"
    WRITE (psfile, '(a)') "  90 rotate 0 -595 translate"
    WRITE (psfile, '(a)') "%%EndPageSetup"
  END SUBROUTINE ps_define_header


  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript procedure block.
  !
  SUBROUTINE ps_define_routines(psfile)
    INTEGER, INTENT(IN) :: psfile   !< PostScript file handle

    WRITE (psfile, '(a)') "% choose font, font height is on stack                                      "
    WRITE (psfile, '(a)') "/chgfont { /hgt exch def  /Helvetica findfont hgt scalefont setfont } def   "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% draw a box, where width,x,y,height are on stack                           "
    WRITE (psfile, '(a)') "/rect {                                                                     "
    WRITE (psfile, '(a)') "  /hgt exch def                                                             "
    WRITE (psfile, '(a)') "  newpath moveto                                                            "
    WRITE (psfile, '(a)') "   0 hgt rlineto 0 rlineto                                                  "
    WRITE (psfile, '(a)') "   0 hgt neg rlineto                                                        "
    WRITE (psfile, '(a)') "  closepath                                                                 "
    WRITE (psfile, '(a)') "  0.5 setlinewidth                                                          "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% annotate if current position < xlimit                                     "
    WRITE (psfile, '(a)') "/annotate {                                                                 "
    WRITE (psfile, '(a)') "  fontsize chgfont                                                          "
    WRITE (psfile, '(a)') "  % copy top of stack to variable                                           "
    WRITE (psfile, '(a)') "  /textstr exch def                                                         "
    WRITE (psfile, '(a)') "  currentpoint exch xlimit                                                  "
    WRITE (psfile, '(a)') "  gt { /xlimit currentpoint pop textstr stringwidth pop unscale mul add def "
    WRITE (psfile, '(a)') "       gsave unscale 1. scale textstr show grestore } if                    "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% annotated, filled box, where x1,x2 are on stack                           "
    WRITE (psfile, '(a)') "/drawintvl {                                                                "
    WRITE (psfile, '(a)') "  % copy top of stack to variables                                          "
    WRITE (psfile, '(a)') "  /thisx2 exch frstcol add def                                              "
    WRITE (psfile, '(a)') "  /thisx1 exch frstcol add def                                              "
    WRITE (psfile, '(a)') "  thisx2 thisx1 sub thisx1 line linehgt mul boxhgt rect                     "
    WRITE (psfile, '(a)') "  gsave 0.8 setgray fill grestore 0 setgray stroke                          "
    WRITE (psfile, '(a)') "  thisx1 line linehgt mul boxhgt add textsep add moveto annotate            "
    WRITE (psfile, '(a)') "  % xlimit = max(x2,xlimit+string width)                                    "
    WRITE (psfile, '(a)') "  xlimit thisx2 lt { /xlimit thisx2 def } if                                "
    WRITE (psfile, '(a)') "  % store rightmost and topmost position                                    "
    WRITE (psfile, '(a)') "  imgwdth thisx2 lt { /imgwdth thisx2 def } if                              "
    WRITE (psfile, '(a)') "  line 1 add linehgt mul dup imghgt gt { /imghgt exch def } if              "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% print grid lines                                                          "
    WRITE (psfile, '(a)') "/drawgrid {                                                                 "
    WRITE (psfile, '(a)') "  0.1 setlinewidth [0.5 1.] 0 setdash                                       "
    WRITE (psfile, '(a)') "  frstcol gridintvl frstcol imgwdth add {                                   "
    WRITE (psfile, '(a)') "   newpath -15 moveto                                                       "
    WRITE (psfile, '(a)') "   0 imghgt 10 add rlineto stroke                                           "
    WRITE (psfile, '(a)') "  } for                                                                     "
    WRITE (psfile, '(a)') "  % print xtick values                                                      "
    WRITE (psfile, '(a)') "  /nstr 7 string def                                                        "
    WRITE (psfile, '(a)') "  0 gridintvl imgwdth {                                                     "
    WRITE (psfile, '(a)') "    /xval exch def                                                          "
    WRITE (psfile, '(a)') "    newpath xval frstcol add 2 add -15 moveto                               "
    WRITE (psfile, '(a)') "    gsave unscale 1. scale scalefct dup scale                               "
    WRITE (psfile, '(a)') "    xval nstr cvs show                                                      "
    WRITE (psfile, '(a)') "    grestore                                                                "
    WRITE (psfile, '(a)') "  } for                                                                     "
    WRITE (psfile, '(a)') "  [] 0 setdash                                                              "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% draw page title                                                           "
    WRITE (psfile, '(a)') "/drawtitle {                                                                "
    WRITE (psfile, '(a)') "  0.1 setlinewidth                                                          "
    WRITE (psfile, '(a)') "  titlehgt chgfont                                                          "
    WRITE (psfile, '(a)') "  0 imghgt titlehgt add moveto                                              "
    WRITE (psfile, '(a)') "  true charpath gsave fill grestore stroke                                  "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% write processor number, where proc string and line no. are on string      "
    WRITE (psfile, '(a)') "/drawprocnum {                                                              "
    WRITE (psfile, '(a)') "  % copy top of stack to variables                                          "
    WRITE (psfile, '(a)') "  /line    exch def                                                         "
    WRITE (psfile, '(a)') "  /procnum exch def                                                         "
    WRITE (psfile, '(a)') "  linehgt 0.5 mul chgfont                                                   "
    WRITE (psfile, '(a)') "  procnum stringwidth pop unscale mul frstcol exch sub 10 unscale mul sub 0 "
    WRITE (psfile, '(a)') "  line linehgt mul add moveto procnum                                       "
    WRITE (psfile, '(a)') "  gsave unscale 1. scale show grestore                                      "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% increment line counter, processor ID is on stack                          "
    WRITE (psfile, '(a)') "/newline {                                                                  "
    WRITE (psfile, '(a)') "  /line line 1 add def                                                      "
    WRITE (psfile, '(a)') "  line drawprocnum                                                          "
    WRITE (psfile, '(a)') "  /xlimit -1 def                                                            "
    WRITE (psfile, '(a)') "} def                                                                       "
    WRITE (psfile, '(a)') "                                                                            "
    WRITE (psfile, '(a)') "% Assemble whole page                                                       "
    WRITE (psfile, '(a)') "/drawpage {                                                                 "
    WRITE (psfile, '(a)') "  /imgwdth  0  def                                                          "
    WRITE (psfile, '(a)') "  /imghgt   0  def                                                          "
    WRITE (psfile, '(a)') "  /unscale  1. def                                                          "
    WRITE (psfile, '(a)') "  /scalefct 1. def                                                          "
    WRITE (psfile, '(a)') "  % draw to determine size                                                  "
    WRITE (psfile, '(a)') "  drawimage                                                                 "
    WRITE (psfile, '(a)') "  20 550 imghgt sub translate                                               "
    WRITE (psfile, '(a)') "  750 imgwdth div dup 1. lt { /scalefct exch def  } if                      "
    WRITE (psfile, '(a)') "  % scale to fit page width                                                 "
    WRITE (psfile, '(a)') "  1. scalefct div /unscale exch def                                         "
    WRITE (psfile, '(a)') "  scalefct 1. scale                                                         "
    WRITE (psfile, '(a)') "  /frstcol frstcol unscale mul def                                          "
    WRITE (psfile, '(a)') "  % redraw scaled image                                                     "
    WRITE (psfile, '(a)') "  erasepage drawgrid drawimage                                              "
    WRITE (psfile, '(a)') "  unscale 1. scale                                                          "
    WRITE (psfile, '(a)') "  thistitle drawtitle                                                       "
    WRITE (psfile, '(a)') "  showpage                                                                  "
    WRITE (psfile, '(a)') "} def                                                                       "
  END SUBROUTINE ps_define_routines


  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript constant definition block.
  !
  !  Change font sizes, line heights and other layout settings here.
  !
  SUBROUTINE ps_define_constants(psfile, title_str)
    INTEGER,          INTENT(IN) :: psfile     !< PostScript file handle
    CHARACTER(LEN=*), INTENT(IN) :: title_str  !< title string

    WRITE (psfile, '(a)') "% constants                                        " 
    WRITE (psfile, '(a)') "/boxhgt     3 def % interval box height            "
    WRITE (psfile, '(a)') "/linehgt   13 def % total row height               "
    WRITE (psfile, '(a)') "/fontsize   5 def % annotation font size           "
    WRITE (psfile, '(a)') "/textsep    2 def % gap between box and annotation "
    WRITE (psfile, '(a)') "/titlehgt   8 def % height of title string         "
    WRITE (psfile, '(a)') "/frstcol   35 def % width of first column          "
    WRITE (psfile, '(a)') "/gridintvl 30 def % grid interval                  "
    WRITE (psfile, '(a,a,a)') "/thistitle (", TRIM(title_str), ") def"
  END SUBROUTINE ps_define_constants


  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript prologue before interval definitions.
  !
  SUBROUTINE ps_define_prologue(psfile)
    INTEGER, INTENT(IN) :: psfile   !< PostScript file handle

    WRITE (psfile, '(a)') "% -------------------------------------------------"
    WRITE (psfile, '(a)') ""
    WRITE (psfile, '(a)') "/drawimage {"
    WRITE (psfile, '(a)') "  /line -1 def    "
  END SUBROUTINE ps_define_prologue


  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript block for a new row in interval diagram.
  !
  SUBROUTINE ps_define_newrow(psfile, row_title)
    INTEGER,          INTENT(IN) :: psfile     !< PostScript file handle
    CHARACTER(LEN=*), INTENT(IN) :: row_title  !< row title string

    WRITE (psfile, '(a,a,a)') "  (", TRIM(row_title), ") newline"
  END SUBROUTINE ps_define_newrow


  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript block defining a single interval.
  !
  SUBROUTINE ps_define_interval(psfile, annotation, sec_start, sec_end)
    INTEGER,          INTENT(IN) :: psfile              !< PostScript file handle
    CHARACTER(LEN=*), INTENT(IN) :: annotation          !< annotation string for this interval
    REAL(wp),         INTENT(IN) :: sec_start, sec_end  !< start, end of this interval

    WRITE (psfile, '(a,a,a,f15.2,a,f15.2,a)') "(", TRIM(annotation), ")  ", &
      &                                       sec_start, "  ", sec_end, "  drawintvl"
  END SUBROUTINE ps_define_interval

  !-------------------------------------------------------------------------------------------------
  !> Generate PostScript epilogue after interval definitions.
  !
  SUBROUTINE ps_define_epilogue(psfile)
    INTEGER, INTENT(IN) :: psfile   !< PostScript file handle

    WRITE (psfile, '(a)') "} def"
    WRITE (psfile, '(a)') "drawpage"
  END SUBROUTINE ps_define_epilogue


  !-------------------------------------------------------------------------------------------------
  !> Enlarge local list of intervals (if necessary).
  !
  SUBROUTINE resize_list(list, new_size)
    TYPE(t_list), INTENT(INOUT)     :: list      !< list of intervals
    INTEGER,      INTENT(IN)        :: new_size  !< new list size
    ! local parameters
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::resize_list'
    TYPE (t_interval), ALLOCATABLE :: tmp(:)
    INTEGER                        :: max_entries, new_max_entries, ierrstat

    ! triangle copy
    IF (ALLOCATED(list%intvl)) THEN
      max_entries = SIZE(list%intvl)
    ELSE
      max_entries      =  0
      list%current_idx = -1 
    ENDIF
    IF (max_entries > 0) THEN
      ALLOCATE(tmp(max_entries), STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE of tmp failed!")
      tmp(1:max_entries) = list%intvl(1:max_entries)
      DEALLOCATE(list%intvl, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE of old array failed!")
    END IF
    new_max_entries = new_size
    ALLOCATE(list%intvl(new_max_entries), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE of new array failed!")
    IF (max_entries > 0) THEN
      list%intvl(1:MIN(new_max_entries, max_entries)) = tmp(1:MIN(new_max_entries, max_entries))
      DEALLOCATE(tmp, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish(routine, "DEALLOCATE of tmp failed!")
    END IF
  END SUBROUTINE resize_list


  !-------------------------------------------------------------------------------------------------
  !> Free local list of intervals.
  !
  SUBROUTINE clear_list(list)
    TYPE(t_list), INTENT(INOUT) :: list !< list of intervals
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::clear_list'
    INTEGER :: ierrstat

    IF (ALLOCATED(list%intvl)) THEN
      DEALLOCATE(list%intvl, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
    list%nentries = 0
  END SUBROUTINE clear_list


  !-------------------------------------------------------------------------------------------------
  !> @return current time in the form of a "t_timestamp" data structure.
  !
  FUNCTION get_timestamp()
    TYPE (t_timestamp) :: get_timestamp !< current time stamp
    ! local variables
    INTEGER :: a, y, m, ivalues(8)

    ! The Fortran routine "date_and_time" returns the following 8 values:
    !
    ! year ; month of the year ; day of the month ; time diff. wrt. UTC ; hour of day ; 
    !  [1]           [2]                [3]               [4]                [5]
    !
    ! minute of hour ; seconds of minute ; milliseconds of second
    !  [6]                  [7]                [8]
    CALL date_and_time(values=ivalues)
    ! compute Julian day out of this time stamp; we ignore the UTC
    ! shift here since were are eventually calculating differences.
    a = FLOOR((14._wp - ivalues(2))/12._wp)
    y = ivalues(1) + 4800 - a
    m = ivalues(2) + 12*a - 3
    
    get_timestamp%julian_day =  &
      &       (ivalues(3) + FLOOR((53._wp*m + 2._wp)/5._wp) + 365._wp*y + FLOOR(y/4.)      &
      &         - FLOOR(y/100._wp) + FLOOR(y/400._wp) - 32045._wp) * 24._wp*3600._wp       &
      &     + ivalues(5) * 3600._wp + ivalues(6)*60._wp + ivalues(7) + ivalues(8)/1000._wp
  END FUNCTION get_timestamp


  !-------------------------------------------------------------------------------------------------
  !> @return current time difference (t1-t2) in seconds.
  !
  FUNCTION get_time_diff(t1, t2)
    REAL(wp) :: get_time_diff                  !< time difference in seconds.
    TYPE (t_timestamp), INTENT(IN) :: t1, t2   !< time stamps

    get_time_diff = t1%julian_day - t2%julian_day;
  END FUNCTION get_time_diff


  !-------------------------------------------------------------------------------------------------
  !> Define initial time stamp used as reference (0 seconds wall-time).
  !
  SUBROUTINE set_reference_time()
    start_time = get_timestamp()
  END SUBROUTINE set_reference_time


  !-------------------------------------------------------------------------------------------------
  !> Append the begin of a new interval to list.
  !
  SUBROUTINE interval_start(name)
    CHARACTER(LEN=*), INTENT(IN) :: name      !< name string for this interval
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::interval_start'
    INTEGER :: nlen

    ! append a new interval to our list of intervals:
    IF (.NOT. ALLOCATED(list%intvl)) THEN
      list%nentries = 1
    ELSE
      list%nentries = list%nentries + 1
    END IF
    CALL resize_list(list, list%nentries)
    ! check if we have currently no open interval:
    IF (list%current_idx /= -1) &
      &  CALL finish(routine, "Internal error: Open interval already exists!")
    list%current_idx = list%nentries
    ! store string for annotation:
    list%intvl(list%current_idx)%annotation = ""
    nlen = MIN(FILENAME_MAX, LEN(name))
    list%intvl(list%current_idx)%annotation(1:nlen) = name(1:nlen)
    ! save current time stamp:
    IF (start_time%julian_day < 0._wp) THEN
      CALL finish(routine, "Internal error: No start time set!")
    ELSE
      list%intvl(list%current_idx)%sec_start = get_time_diff(get_timestamp(), start_time)
    END IF
  END SUBROUTINE interval_start


  !-------------------------------------------------------------------------------------------------
  !> Close the currently opened interval in the list.
  !
  SUBROUTINE interval_end(name)
    CHARACTER(LEN=*), INTENT(IN) :: name      !< name string for this interval
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::interval_end'
    INTEGER :: nlen

    ! check if there is an interval currently open:
    IF (list%current_idx == -1) &
      &  CALL finish(routine, "Internal error: No interval currently open!")
    ! check if the name string matches the name of the currently
    ! opened interval:
    nlen = MIN(FILENAME_MAX, LEN(name))
    IF (name(1:nlen) /= list%intvl(list%current_idx)%annotation(1:nlen)) THEN
      CALL finish(routine, "Internal error: String does not match the currently opened interval!")
    END IF
    ! save current time stamp:
    IF (start_time%julian_day == -1) THEN
      CALL finish(routine, "Internal error: No start time set!")
    ELSE
      list%intvl(list%current_idx)%sec_end = get_time_diff(get_timestamp(), start_time)
    END IF
    list%current_idx = -1 ! current interval is now complete
  END SUBROUTINE interval_end


  !-------------------------------------------------------------------------------------------------
  !> Collect intervals from multiple MPI tasks on a single MPI task.
  !
  SUBROUTINE interval_gather_procs(iroot, mpi_comm, this_pe_name, tot_list)
    INTEGER,          INTENT(IN)    :: iroot         !< MPI root process rank
    INTEGER,          INTENT(IN)    :: mpi_comm      !< MPI process communicator
    CHARACTER(LEN=MAX_CHAR_LENGTH), INTENT(IN)    :: this_pe_name  !< local PE name string
    TYPE(t_global_list),            INTENT(INOUT) :: tot_list      !< gathered result
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::interval_gather_procs'
    INTEGER :: np, ierrstat, this_pe, istart, iend, max_size(1), pe, asize
    INTEGER,                     ALLOCATABLE :: nentries(:)
    REAL(wp),                    ALLOCATABLE :: real_buf(:)
    CHARACTER(LEN=FILENAME_MAX), ALLOCATABLE :: char_buf(:)

    ! get local MPI rank
    this_pe = p_comm_rank(mpi_comm)
    ! get the total no. of list entries
    np = p_comm_size(mpi_comm)
    ALLOCATE(nentries(np), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! handle the case that this PE did not have any output work:
    IF (.NOT. ALLOCATED(list%intvl))  list%nentries = 0

    CALL p_gather(list%nentries, nentries, iroot, mpi_comm)

    max_size(1) = list%nentries
    IF (this_pe == iroot) THEN
      max_size = MAXVAL(nentries)
      tot_list%list%nentries    = SUM(nentries(:))
      tot_list%list%current_idx = -1
      CALL resize_list(tot_list%list, tot_list%list%nentries)
      ! allocate array is "t_global_list" data structure
      asize = np
    ELSE
      ! MoHa: tot_list%pe_names needs to be allocated, otherwise it must not be
      !       passed to p_gather
      asize = 1
    END IF
    ! allocate arrays:
    ALLOCATE(tot_list%start_idx(asize), tot_list%pe_names(asize), &
      &      real_buf(max_size(1)), char_buf(max_size(1)), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    ! gather PE names:
    CALL p_gather(this_pe_name, tot_list%pe_names, iroot, mpi_comm)

    IF (this_pe /= iroot) THEN
      DEALLOCATE(tot_list%pe_names, tot_list%start_idx, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF

    IF (this_pe == iroot) THEN
      ! loop over all processes
      iend = 0
      DO pe=0,(np-1)
        istart = iend + 1
        iend   = istart+nentries(pe+1)-1
        tot_list%start_idx(pe+1) = istart
        IF (nentries(pe+1) == 0) CYCLE
        IF (pe == this_pe) THEN
          ! local copy 
          tot_list%list%intvl(istart:iend) = list%intvl(1:list%nentries)
        ELSE
          ! WRITE (0,*) "Receive ", nentries(pe+1), " entries from PE ", pe
          ! receive data
          CALL p_recv(real_buf, pe, 0, nentries(pe+1), mpi_comm)
          tot_list%list%intvl(istart:iend)%sec_start = real_buf(1:nentries(pe+1))
          CALL p_recv(real_buf, pe, 1, nentries(pe+1), mpi_comm)
          tot_list%list%intvl(istart:iend)%sec_end = real_buf(1:nentries(pe+1))
          CALL p_recv(char_buf, pe, 2, nentries(pe+1), mpi_comm)
          tot_list%list%intvl(istart:iend)%annotation = char_buf(1:nentries(pe+1))
        END IF
      END DO
    ELSE IF (list%nentries > 0) THEN
      ! WRITE (0,*) "Send ", list%nentries, " entries to iroot"
      ! send data
      real_buf(1:list%nentries) = list%intvl(1:list%nentries)%sec_start
      CALL p_send(real_buf, iroot, 0, list%nentries, mpi_comm)
      real_buf(1:list%nentries) = list%intvl(1:list%nentries)%sec_end
      CALL p_send(real_buf, iroot, 1, list%nentries, mpi_comm)
      char_buf(1:list%nentries) = list%intvl(1:list%nentries)%annotation
      CALL p_send(char_buf, iroot, 2, list%nentries, mpi_comm)
    END IF
    ! clean up
    DEALLOCATE(nentries, real_buf, char_buf, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE interval_gather_procs


  !-------------------------------------------------------------------------------------------------
  !> Collect intervals from other MPI tasks on rank #0 and write
  !  PostScript file.
  !
  SUBROUTINE interval_write_psfile(filename, title_str, this_pe_name, mpi_comm)
    CHARACTER(LEN=*),               INTENT(IN) :: filename      !< file name string (only used on PE#0)
    CHARACTER(LEN=*),               INTENT(IN) :: title_str     !< title string     (only used on PE#0)
    CHARACTER(LEN=*),               INTENT(IN) :: this_pe_name  !< local PE name (row title)
    INTEGER,                        INTENT(IN) :: mpi_comm      !< MPI process communicator
    !> PE rank that collects and writes the output:
    INTEGER, PARAMETER :: iroot = 0 
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::interval_write_psfile'
    INTEGER                        :: psfile, i, j, start_idx, end_idx, ierrstat
    TYPE (t_global_list)           :: global_list
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: pe_name

    ! --- collect all process data on rank #0:
    pe_name = this_pe_name
    CALL interval_gather_procs(iroot, mpi_comm, pe_name, global_list)

    ! --- generate PostScript file:
    IF (ALLOCATED(global_list%pe_names)) THEN
      psfile = find_next_free_unit(10,100)
      OPEN(psfile, file=TRIM(filename))
      
      CALL ps_define_header(psfile)
      CALL ps_define_routines(psfile)
      CALL ps_define_constants(psfile, title_str)
      CALL ps_define_prologue(psfile)
      ! loop over all PEs
      DO j=1,SIZE(global_list%pe_names)
        CALL ps_define_newrow(psfile, global_list%pe_names(j))
        ! loop over all intervals for PE "i":
        start_idx = global_list%start_idx(j)
        end_idx   = global_list%list%nentries
        IF (j<SIZE(global_list%pe_names))  end_idx = global_list%start_idx(j+1)-1
        DO i=start_idx,end_idx
          CALL ps_define_interval(psfile, TRIM(global_list%list%intvl(i)%annotation), &
            &                     global_list%list%intvl(i)%sec_start,           &
            &                     global_list%list%intvl(i)%sec_end)
        END DO
      END DO
      CALL ps_define_epilogue(psfile)
      
      CLOSE(psfile)
    END IF

    ! --- clean up
    CALL clear_list(list)
    IF (ALLOCATED(global_list%pe_names))  THEN
      DEALLOCATE(global_list%pe_names, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
    IF (ALLOCATED(global_list%start_idx))  THEN
      DEALLOCATE(global_list%start_idx, STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
  END SUBROUTINE interval_write_psfile

END MODULE mo_name_list_output_stats

