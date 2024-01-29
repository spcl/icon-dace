!
! mo_art_fplume_utilities
! This module contains all the utility subroutines for FPLUME.
! Original code by A. Folch, G. Macedonio, A. Costa (2016)
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

MODULE mo_art_fplume_utilities
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish, message, message_text
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: get_input_npar, get_input_rea, get_input_cha,        &
         &  get_granulometry_nclass, get_granulometry_value

  INTEGER, PARAMETER  :: s_mess = 256
  INTEGER, PARAMETER  :: s_long = 512
  INTEGER, PARAMETER  :: nwormax = 128
  INTEGER, PARAMETER  :: nparmax = 128

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_granulometry_nclass(fname,nc,istat,mess)
  !**************************************************************************
  !*
  !*    Gets the number of granulometric classes
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*
  !*    OUTPUT:
  !*    integer         nc       Number of granulometric classes
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !**************************************************************************
  IMPLICIT NONE
  !
  CHARACTER(LEN=*)    :: fname, mess
  INTEGER             :: nc,istat
  !
  CHARACTER(LEN=s_mess) :: mymessage
  INTEGER               :: msglen
  !
  !***  Initializations
  !
  msglen = LEN(mess)
  mess(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  !
  OPEN(90,FILE=fname(1:LEN_TRIM(fname)),STATUS='old',ERR=100)
  !
  !***  Reads
  !
  READ(90,*,ERR=101) nc
  CLOSE(90)
  !
  !***  Successful end
  RETURN
  !
  !***  List of errors
100 istat = -1
  mymessage = 'get_granulometry_nclass: error opening the '//'granulometry file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
101 istat = -1
  mymessage = 'get_granulometry_nclass: error reading the '//'granulometry file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
  !
END SUBROUTINE get_granulometry_nclass
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_granulometry_value(fname,word,val,istat,mess)
  !**************************************************************************
  !*
  !*    Gets the a granulometric property
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the granulometry file
  !*    character*(*)   word     Property to extract. Possible values are
  !*                             'DIAMETER' (in mm)
  !*                             'DENSITY'
  !*                             'FRACTION'
  !*                             'SPHERICITY'
  !*
  !*    OUTPUT:
  !*    real            val(nc)  Values of the property
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character*(*)   message  Exit message
  !*
  !*************************************************************************
  IMPLICIT NONE
  !
  CHARACTER(LEN=*)  :: fname,mess,word
  INTEGER           :: nc,istat
  REAL   (wp)       :: val(*)
  !
  CHARACTER(LEN=s_mess) :: mymessage
  INTEGER               :: msglen,ipos,ic
  REAL   (wp)           :: rvoid(4)
  !
  !***  Initializations
  !
  msglen = LEN(mess)
  mess(:)  = ' '
  istat = 0
  !
  !***  Checks that word exists in the dictionary
  !
  IF (word(1:LEN_TRIM(word))=='DIAMETER') THEN
    ipos = 1
  ELSE IF (word(1:LEN_TRIM(word))=='DENSITY') THEN
    ipos = 2
  ELSE IF (word(1:LEN_TRIM(word))=='SPHERICITY') THEN
    ipos = 3
  ELSE IF (word(1:LEN_TRIM(word))=='FRACTION') THEN
    ipos = 4
  ELSE
    istat = -1
    mymessage = 'get_granulometry_value: word not found in the dictionary'
    mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
    RETURN
  END IF
  !
  !***  Opens the file
  OPEN(90,FILE=TRIM(fname),STATUS='old',ERR=100)
  !
  !***  Reads
  READ(90,*,ERR=101) nc
  DO ic = 1,nc
    READ(90,*,ERR=101) rvoid(1),rvoid(2),rvoid(3),rvoid(4)
    val(ic) = rvoid(ipos)
  END DO
  CLOSE(90)
  !
  !***  Successful end
  RETURN
  !
  !***  List of errors
100 istat = -1
  mymessage = 'get_granulometry_value: error opening the '//'granulometry file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
101 istat = -1
  mymessage = 'get_granulometry_value: error reading the '//'granulometry file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
  !
END SUBROUTINE get_granulometry_value
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_input_int(fname,sblock,line,rvalue,nval,istat,mess)
  !**************************************************************************
  !*
  !*    Gets nval integer inputs from the file fname
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   sblock    Block to search
  !*    character*(*)   line     Line block to search
  !*    integer         nval     Number of integers to read
  !*
  !*    OUTPUT:
  !*    integer          istat    -1 ERROR  0 OK  1 WARNING
  !*    integer          rvalue    Values of the nval integers read
  !*    character        message  Output message with the error description
  !*
  !**************************************************************************
  IMPLICIT NONE
  CHARACTER(LEN=*)      :: mess,fname,sblock,line
  INTEGER               :: nval,istat
  INTEGER               :: rvalue(nval)
  !
  CHARACTER(LEN=s_mess) :: mymessage
  CHARACTER(LEN=s_long) :: card
  CHARACTER(LEN=s_long),DIMENSION(128) :: words(nwormax)
  REAL(wp)              :: words_all(nparmax)
  LOGICAL               :: linefound,blockfound
  INTEGER               :: msglen,nword,npar,ival,ipar,nword_old
  REAL(wp)              :: param(nparmax),x0,xf,dx
  !
  !***  Initializations
  msglen = LEN(mess)
  mess(:)  = ' '
  words(:)(:) = ' '
  istat = 0
  !
  !***  Opens the file
  OPEN(90,FILE=TRIM(fname),STATUS='old',ERR=101)
  !
  !***  Search the line
  blockfound = .false.
  linefound  = .false.
  DO WHILE (.NOT.linefound)
    DO WHILE (.NOT.blockfound)
      READ(90,'(a256)',END=102) card
      CALL sdecode(card,words,param,nword,npar)
      IF (words(1)(1:LEN_TRIM(sblock))==sblock(1:LEN_TRIM(sblock))) blockfound=.true.
    END DO
    READ(90,'(a256)',END=103) card
    CALL sdecode(card,words,param,nword,npar)
!    IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) linefound = .true.
      IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) THEN
        linefound = .true.
        words_all(1:npar) = param
        nword_old = npar +1             ! if last symbol is '/' --> overwrite in next step, first is 'PHASES'
        DO WHILE (TRIM(words(nword))=='/') ! neu
          READ(90,'(a256)',END=102) card
          CALL sdecode(card,words,param,nword,npar)
          words_all(nword_old:(nword_old+npar)) = param(1:npar)
!          IF (words(nword)=='/') THEN
!            words_all(nword_old:(nword_old+nword-2)) = words(1:nword-1)
!          ELSE
!            words_all(nword_old:(nword_old+nword-1)) = words(1:nword)
!          ENDIF
          nword_old= nword_old+npar
        ENDDO
      ENDIF

  END DO
  !
  !***  Line format  FROM x0 TO xf INCREMENT dx
  !***  Calculate npar
  !
  IF (TRIM(words(2))=='FROM'.OR.TRIM(words(2))=='from') THEN
    x0 = param(1)
    xf = param(2)
    IF (x0>xf) GOTO 104
    dx = MIN(xf-x0,param(3))
    npar = INT((xf-x0)/dx)+1
    IF (npar>nparmax) THEN
      npar = nparmax
      istat = 1
      mymessage = 'get_input_int: warning. Too big number of parameters'
      mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
    END IF
    DO ipar = 1,npar
      param(ipar) = x0 + (ipar-1)*dx
    END DO
  END IF
  !
  npar = nword_old - 1
  IF (npar<nval) GOTO 105
  !
  DO ival = 1,nval
    rvalue(ival) = INT(words_all(ival))
  END DO
  !
  !***  Successful end
  CLOSE(90)
  RETURN
  !
  !***  List of errors
101 istat = -1
  CLOSE(90)
  mymessage = 'get_input_int: error opening the input file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
102 istat = -1
  CLOSE(90)
  mymessage = 'get_input_int: sblock '//TRIM(sblock)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
103 istat = -1
  CLOSE(90)
  mymessage = 'get_input_int: line '//TRIM(line)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
104 istat = -1
  CLOSE(90)
  mymessage = 'get_input_int: error in line '//line(1:LEN_TRIM(line))
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
105 istat = 1
  CLOSE(90)
  mymessage = 'get_input_int: too few parameters in line '//TRIM(line)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
  !
END SUBROUTINE get_input_int
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_input_rea(fname,sblock,line,rvalue,nval,istat,mess)
  !**************************************************************************
  !*
  !*    Gets nval real inputs from the file fname
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   sblock    Block to search
  !*    character*(*)   line     Line block to search
  !*    integer         nval     Number of integers to read
  !*
  !*    OUTPUT:
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    real            rvalue    Values of the nval integers read
  !*    character       message  Output message with the error description
  !*
  !**************************************************************************
  IMPLICIT NONE
  CHARACTER(LEN=*)      :: mess,fname,sblock,line
  INTEGER               :: nval,istat
  REAL(wp)              :: rvalue(nval)
  !
  CHARACTER(LEN=s_mess) :: mymessage
  CHARACTER(LEN=s_long) :: card
  CHARACTER(LEN=s_long),DIMENSION(128) :: words(nwormax)
  REAL(wp)              :: words_all(nparmax)
  LOGICAL               :: linefound,blockfound
  INTEGER               :: msglen,nword,npar,ival,ipar,nword_old
  REAL(wp)              :: param(nparmax),x0,xf,dx
  !
  !***  Initializations
  msglen = LEN(mess)
  mess(:)  = ' '
  words(:)(:) = ' '
  istat = 0
  !
  !***  Opens the file
  OPEN(90,FILE=TRIM(fname),STATUS='old',ERR=101)
  !
  !***  Search the line
  blockfound = .FALSE.
  linefound  = .FALSE.
  DO WHILE (.NOT.linefound)
    DO WHILE (.NOT.blockfound)
      READ(90,'(a256)',END=102) card
      CALL sdecode(card,words,param,nword,npar)
      IF (words(1)(1:LEN_TRIM(sblock))==sblock(1:LEN_TRIM(sblock))) blockfound=.TRUE.
    END DO
    READ(90,'(a256)',END=103) card
    CALL sdecode(card,words,param,nword,npar)
!      IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) linefound = .true.
      IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) THEN
        linefound = .TRUE.
        words_all(1:npar) = param(1:npar)
        nword_old = npar+1             ! if last symbol is '/' --> overwrite in next step, first is 'PHASES'
        DO WHILE (TRIM(words(nword))=='/') ! neu
          READ(90,'(a256)',END=102) card
          CALL sdecode(card,words,param,nword,npar)
          words_all(nword_old:(nword_old+npar)) = param
!          IF (words(nword)=='/') THEN
!            words_all(nword_old:(nword_old+nword-2)) = words(1:nword-1)
!          ELSE
!            words_all(nword_old:(nword_old+nword-1)) = words(1:nword)
!          ENDIF
          nword_old= nword_old+npar
        ENDDO
      ENDIF
  END DO
  !
  !***  Line format  FROM x0 TO xf INCREMENT dx
  !***  Calculate npar
  IF (TRIM(words(2))=='FROM'.OR.TRIM(words(2))=='from') THEN
    x0 = param(1)
    xf = param(2)
    IF (x0>xf) GOTO 104
    dx = MIN(xf-x0,param(3))
    npar = INT((xf-x0)/dx)+1
    IF (npar>nparmax) THEN
      npar = nparmax
      istat = 1
      mymessage = 'get_input_rea: warning. Too big number of parameters'
      mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
    END IF
    DO ipar = 1,npar
      param(ipar) = x0 + (ipar-1)*dx
    END DO
  END IF
  !
  npar = nword_old - 1
  IF (npar<nval) GOTO 105
  !
  DO ival = 1,nval
    rvalue(ival) = words_all(ival)
  END DO
  !
  !***  Successful end
  CLOSE(90)
  RETURN
  !
101 istat = -1
  CLOSE(90)
  mymessage = 'get_input_rea: error opening the input file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
102 istat = -1
  CLOSE(90)
  mymessage = 'get_input_rea: sblock '//TRIM(sblock)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
103 istat = -1
  CLOSE(90)
  mymessage = 'get_input_rea: line '//TRIM(line)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
104 istat = -1
  CLOSE(90)
  mymessage = 'get_input_rea: error in line '//line(1:LEN_TRIM(line))
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
105 istat = -1
  CLOSE(90)
  mymessage = 'get_input_rea: too few parameters in line '//TRIM(line)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
END SUBROUTINE get_input_rea
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_input_cha(fname,sblock,line,rvalue,nval,istat,mess)
  !**************************************************************************
  !*
  !*    Gets nval character inputs from the file fname
  !*
  !*    NOTE: words are converted to UPPER case
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   sblock   Block to search
  !*    character*(*)   line     Line block to search
  !*    integer         nval     Number of integers to read
  !*
  !*    OUTPUT:
  !*    integer          istat    -1 ERROR  0 OK  1 WARNING
  !*    character        rvalue    Values of the nval integers read
  !*    character        message  Output message with the error description
  !*
  !**************************************************************************
  IMPLICIT NONE
  INTEGER               :: nval,istat
  CHARACTER(LEN=*)      :: mess,fname,sblock,line
  CHARACTER(LEN=*),DIMENSION(100)      :: rvalue(nval)
  !
  CHARACTER(LEN=s_mess) :: mymessage
  CHARACTER(LEN=s_long) :: card
  CHARACTER(LEN=s_long),DIMENSION(128) :: words(nwormax)
  CHARACTER(LEN=s_long),DIMENSION(128) :: words_all(nwormax)
  LOGICAL               :: linefound,blockfound
  INTEGER               :: msglen,nword,npar,ival,j,nword_old
  REAL(wp)              :: param(nparmax)
  !
  !***  Initializations
  msglen = LEN(mess)
  mess(:)  = ' '
  istat = 0
  !
  !***  Opens the file
  OPEN(90,FILE=TRIM(fname),STATUS='old',ERR=101)
  !
  !***  Search the line
  blockfound = .FALSE.
  linefound  = .FALSE.
  DO WHILE (.NOT.linefound)
    DO WHILE (.NOT.blockfound)
      READ(90,'(a256)',END=102) card
      CALL sdecode(card,words,param,nword,npar)! param,npar leer im Fall von get_input_cha
      IF (words(1)(1:LEN_TRIM(sblock))==sblock(1:LEN_TRIM(sblock))) blockfound=.TRUE.
    END DO
    READ(90,'(a256)',END=103) card
    CALL sdecode(card,words,param,nword,npar)     
    IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) THEN 
      linefound = .TRUE.
      words_all(1:nword) = words(1:nword) 
      nword_old = nword ! if last symbol is '/' --> overwrite in next step, first is 'PHASES'
      DO WHILE (words(nword)=='/') ! neu ueberschreiben      
        READ(90,'(a256)',END=102) card 
        CALL sdecode(card,words,param,nword,npar)    
        IF (words(nword)=='/') THEN
          words_all(nword_old:(nword_old+nword-2)) = words(1:nword-1)(1:19) 
        ELSE
          words_all(nword_old:(nword_old+nword-1)) = words(1:nword)(1:19)  
        ENDIF
        nword_old= nword_old+nword-1   
      ENDDO  
    ENDIF
  END DO
  !
  IF (line/='PHASES') THEN
    IF ((nword-1)<nval) GOTO 104
  ENDIF
  !
  DO ival = 1,nval
    rvalue(ival)(1:LEN_TRIM(words_all(ival+1))) = words_all(ival+1)(1:LEN_TRIM(words_all(ival+1)))
    !
    !***     Fill the rest with ' '
    DO j=LEN_TRIM(words_all(ival+1))+1,LEN(rvalue(ival))
      rvalue(ival)(j:j)=' '
    END DO
  END DO
  !
  !***  Convert to upper case
  DO ival = 1,nval
    CALL upcase(rvalue(ival))
  END DO
  !
  !***  Successful end
  CLOSE(90)
  RETURN
  !
  !***  List of errors
101 istat = -1
  CLOSE(90)
  mymessage = 'get_input_cha: error opening the input file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
102 istat = -1
  CLOSE(90)
  mymessage = 'get_input_cha: sblock '//TRIM(sblock)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
103 istat = -2
  CLOSE(90)
  mymessage = 'get_input_cha: line '//TRIM(line)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
104 istat = -1
  CLOSE(90)
  mymessage = 'get_input_cha: too few parameters in line '//TRIM(line)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
  !
END SUBROUTINE get_input_cha
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_input_npar(fname,sblock,line,nword_old,istat,mess) !nword_old ist npar
  !**************************************************************************
  !*
  !*    Gets the number of parameters associated to a line
  !*
  !*    INPUT:
  !*    character*(*)   fname    Name of the file
  !*    character*(*)   sblock   Block to search
  !*    character*(*)   line     Line block to search
  !*
  !*    OUTPUT:
  !*    integer         npar     Number of parameters
  !*    integer         istat    -1 ERROR  0 OK  1 WARNING
  !*    character       message  Output message with the error description
  !*
  !**************************************************************************
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: mess,fname,sblock,line
  INTEGER          :: npar,istat,ival
  !
  CHARACTER(LEN=s_mess) :: mymessage
  CHARACTER(LEN=s_long) :: card
  CHARACTER(LEN=s_long),DIMENSION(128) :: words(nwormax)
  REAL(wp)              :: words_all(nparmax)
  LOGICAL               :: linefound,blockfound
  INTEGER               :: msglen,nword,nword_old
  REAL(wp)              :: param(nparmax),x0,xf,dx
  !
  !***  Initializations
  msglen = LEN(mess)
  mess(:)  = ' '
  words(:)(:) = ' '
  istat = 0
  !
  !***  Opens the file
  OPEN(90,FILE=TRIM(fname),STATUS='old',ERR=101)
  !
  !***  Search the line
  blockfound = .FALSE.
  linefound  = .FALSE.
  DO WHILE (.NOT.linefound)
    DO WHILE (.NOT.blockfound)
      READ(90,'(a256)',END=102) card
      CALL sdecode(card,words,param,nword,npar)
      IF (words(1)(1:LEN_TRIM(sblock)) == sblock(1:LEN_TRIM(sblock))) blockfound=.TRUE.
    END DO
    READ(90,'(a256)',END=103) card
    CALL sdecode(card,words,param,nword,npar)
!    IF (words(1)(1:LEN_TRIM(line)) == line(1:LEN_TRIM(line))) linefound = .true.
      IF (words(1)(1:LEN_TRIM(line))==line(1:LEN_TRIM(line))) THEN
        linefound = .TRUE.
        words_all(1:npar) = param(1:npar)
        nword_old = npar             ! if last symbol is '/' --> overwrite in next step, first is 'PHASES'
        DO WHILE (words(nword)=='/') ! neu
          READ(90,'(a256)',END=102) card
          CALL sdecode(card,words,param,nword,npar)
            words_all(nword_old:(nword_old+npar)) = param(1:npar)
!          IF (words(nword)=='/') THEN
!            words_all(nword_old:npar) = words(1:nword-1)
!          ELSE
!            words_all(nword_old:(nword_old+nword-1)) = words(1:nword)
!          ENDIF
          nword_old= nword_old+npar
        ENDDO
      ENDIF
  nword_old=nword_old
  END DO
  !
  !***  Line format  FROM x0 TO xf INCREMENT dx
  !***  Calculate npar
  IF (TRIM(words(2)) == 'FROM' .OR. TRIM(words(2)) == 'from') THEN
    x0 = param(1)
    xf = param(2)
    IF (x0>xf) GOTO 104
    dx = MIN(xf-x0,param(3))
    npar = INT((xf-x0)/dx)+1
    IF (npar>nparmax) THEN
      npar = nparmax
      istat = 1
      mymessage = 'get_input_int4: warning. Too big number of parameters'
      mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
    END IF
  END IF
  !
  npar = nword_old - 1
  !
  !***  Successful end
  CLOSE(90)
  RETURN
  !
  !***  List of errors
101 istat = -1
  CLOSE(90)
  mymessage = 'get_input_npar : error opening the input file '//TRIM(fname)
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
102 istat = -1
  CLOSE(90)
  mymessage = 'get_input_npar : sblock '//TRIM(sblock)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
103 istat = -1
  CLOSE(90)
  mymessage = 'get_input_npar : line '//TRIM(line)//' not found in the input file'
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
104 istat = -1
  CLOSE(90)
  mymessage = 'get_input_npar: error in line '//line(1:LEN_TRIM(line))
  mess(1:MIN(msglen,s_mess)) = mymessage(1:MIN(msglen,s_mess))
  RETURN
  !
END SUBROUTINE get_input_npar
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE sdecode(card,words,param,nword,npar)
  !********************************************************************
  !*
  !*    This routine decodes a string card(s_long) into words and parameters
  !*
  !********************************************************************
  IMPLICIT NONE
  INTEGER      ::    nword,npar
  CHARACTER(LEN=s_long) ::  card
  CHARACTER(LEN=s_long) ::  words(nwormax)
  CHARACTER(LEN=1)      ::  sstring(s_long)
  REAL(wp)              ::  param(nparmax)
  !
  INTEGER     ::  ipos,first,last,nstr,lflag,i,ii
  REAL(wp)    ::  digit
  !
  !***  Initializations
  nword = 0
  npar  = 0
  ipos  = 0

  DO
    ipos = ipos + 1
    IF (ipos>s_long) RETURN
11  IF (card(ipos:ipos)==' '.OR.card(ipos:ipos)=='=') THEN  ! kein Wert an Position
      ipos = ipos + 1
      IF (ipos>s_long) RETURN
      GOTO 11
    END IF
    first = ipos
    !
    ipos = ipos + 1
    IF (ipos>s_long) RETURN
21  IF (card(ipos:ipos)/=' '.AND.card(ipos:ipos)/='=') THEN ! lese den character an Position ipos
      ipos = ipos + 1
      IF (ipos>s_long) RETURN
      GOTO 21
    END IF
    last = ipos-1
    !
    nstr = last-first+1 !Laenge des Parameters
    !
    ii = 0
    DO i=first,last
      ii = ii + 1
      sstring(ii) = card(i:i)
    END DO
    CALL decod1(sstring,nstr,lflag,digit)
    IF (lflag==0) THEN !number
      npar = npar + 1
      param(npar)= digit
    ELSE IF (lflag==1) THEN !string
      nword = nword + 1
      words(nword)(:) = ' '
      words(nword)(1:nstr) = card(first:last) !Speichern des parameters
    END IF
    !
  END DO
  RETURN
END SUBROUTINE sdecode
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE upcase(word)
  !***********************************************************************
  !*
  !*    This routine converts word to upper case
  !*
  !***********************************************************************
  IMPLICIT NONE
  CHARACTER(LEN=*) :: word
  INTEGER          :: iposi,ioctv,item1,item2,item3
  INTEGER          :: msglen
  !
  item1 = INT(o'141')
  item2 = INT(o'172')
  item3 = INT(o'40')
  !
  msglen=LEN_TRIM(word)
  !
  DO iposi=1,msglen                                 ! process all positions
    ioctv=ICHAR(word(iposi:iposi))                 ! octal value
    IF (item1<=ioctv.AND.item2>=ioctv) THEN        ! it is a lower case
      ioctv=ioctv-item3                           ! equivalent upper case
      word(iposi:iposi)=CHAR(ioctv)               ! convert it to upcase
    END IF
  END DO ! iposi=1,msglen
  RETURN
END SUBROUTINE upcase
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE decod1(string,nstr,lflag,digit)
  !*******************************************************************
  !*
  !*    This subroutine decodes a single string(1:nstr)
  !*
  !*    If string(1:nstr) is a string returns lflag = 1
  !*    If string(1:nstr) is a number returns lflag = 0 and the digit
  !*
  !*******************************************************************
  IMPLICIT NONE
  INTEGER     :: nstr,lflag
  INTEGER     :: istr,decim
  REAL(wp)    :: digit
  CHARACTER(LEN=1) :: string(s_long)
  !
  lflag = 0                                             ! Number by default
  istr  = 1
  !
  DO WHILE (istr<=nstr)
    decim = ICHAR(string(istr))                         ! decimal value
    IF (decim<48.OR.decim>57) THEN                      ! It is not a num.
      IF (decim/=43.AND.decim/=45.AND.  &              ! discard + -
        & decim/=68.AND.decim/=69.AND.  &              ! discard D E
        & decim/=100.AND.decim/=101.AND.  &            ! discard d e
        & decim/=46) THEN                              ! discard .
        lflag = 1
        istr  = nstr
       END IF
     END IF
       istr = istr+1
  END DO
  !
  IF (lflag==0) digit = stof(string,nstr)               ! It's a number
  RETURN
END SUBROUTINE decod1
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION stof(string,nstr)
  !**************************************************************
  !*
  !*    This routine converts a real/integer number stored in a
  !*    string(1:nstr) into a real(wp) digit format
  !*
  !**************************************************************
  IMPLICIT NONE
  INTEGER          :: nstr
  CHARACTER(LEN=1) :: string(*)
  !
  INTEGER          :: i,ipos,nsign,esign,nvalu
  INTEGER          :: expo,valu(s_long)
  LOGICAL          :: next
  !
  stof = 0.0_wp
  !
  !***  Sing decoding
  ipos = 1
  IF (ichar(string(ipos))==43) THEN         !  + sign
    nsign = 1
    ipos  = ipos + 1
  ELSE IF (ichar(string(ipos))==45) THEN    !  - sign
    nsign = -1
    ipos  = ipos + 1
  ELSE                                      !  no sing (+)
    nsign = 1
    ipos  = ipos
  END IF
  !
  !***  Base decoding
  nvalu = 0
  next  = .TRUE.
  DO while(next)
    IF ((ICHAR(string(ipos))==68 ).OR. &         ! D
      & (ICHAR(string(ipos))==69 ).OR. &         ! E
      & (ICHAR(string(ipos))==100).OR. &         ! d
      & (ICHAR(string(ipos))==101).OR. &         ! e
      & (ICHAR(string(ipos))==46 )) THEN         ! .
      next = .FALSE.
    ELSE
      nvalu = nvalu + 1
      valu(nvalu) = stof1(string(ipos))
      ipos = ipos + 1
      IF (ipos==(nstr+1)) THEN
        next = .FALSE.
        ipos = ipos - 1
      END IF
    END IF
  END DO
  DO i = 1,nvalu
    stof = stof + valu(i)*1d1**(nvalu-i)
  END DO
  !
  !***  Decimal decoding
  IF ((ichar(string(ipos))==46   ).AND.  &
    &  ipos  /=nstr) THEN
    ipos = ipos + 1
    nvalu = 0
    next  = .TRUE.
    DO WHILE (next)
      IF ((ichar(string(ipos))==68 ).OR. &       ! D
        & (ichar(string(ipos))==69 ).OR. &       ! E
        & (ichar(string(ipos))==100).OR. &       ! d
        & (ichar(string(ipos))==101)) THEN       ! e
        next = .FALSE.
      ELSE
        nvalu = nvalu + 1
        valu(nvalu) = stof1(string(ipos))
        ipos = ipos + 1
        IF (ipos==(nstr+1)) THEN
          next = .FALSE.
          ipos = ipos - 1
        END IF
      END IF
    END DO
    DO i = 1,nvalu
      stof = stof + valu(i)*1d1**(-i)
    END DO
  END IF
  !
  !***  Exponent
  IF (((ichar(string(ipos))==68 ).OR. &         ! D
    &  (ichar(string(ipos))==69 ).OR. &         ! E
    &  (ichar(string(ipos))==100).OR. &         ! d
    &  (ichar(string(ipos))==101)).AND. &       ! e
    &   ipos  /=nstr) THEN
    ipos = ipos + 1
    IF (ichar(string(ipos))==43) THEN          !  + sign
      esign = 1
      ipos  = ipos + 1
    ELSE IF (ichar(string(ipos))==45) THEN     !  - sign
      esign = -1
      ipos  = ipos + 1
    ELSE                                       !  no sing (+)
      esign = 1
      ipos  = ipos
    END IF
    !
    nvalu = 0
    next  = .TRUE.
    DO WHILE (next)
      nvalu = nvalu + 1
      valu(nvalu) = stof1(string(ipos))
      ipos = ipos + 1
      IF (ipos==(nstr+1)) THEN
        next = .FALSE.
        ipos = ipos - 1
      END IF
    END DO
    expo = 0
    DO i = 1,nvalu
      expo = expo + valu(i)*10**(nvalu-i)
    END DO
    !
    IF (esign==1) THEN
      stof = stof*(10.0_wp**expo)
    ELSE IF (esign==-1) THEN
      stof = stof/(10.0_wp**expo)
    END IF
    !
  END IF
  !
  stof = nsign*stof
END FUNCTION stof
!!
!!-------------------------------------------------------------------------
!!
INTEGER FUNCTION stof1(string1)
  !**************************************************************
  !*
  !*    Decodes a character*1 string
  !*
  !**************************************************************
  IMPLICIT NONE
  CHARACTER(LEN=1) :: string1
  !
  IF (string1=='0') THEN
    stof1 = 0
  ELSE IF (string1=='1') THEN
    stof1 = 1
  ELSE IF (string1=='2') THEN
    stof1 = 2
  ELSE IF (string1=='3') THEN
    stof1 = 3
  ELSE IF (string1=='4') THEN
    stof1 = 4
  ELSE IF (string1=='5') THEN
    stof1 = 5
  ELSE IF (string1=='6') THEN
    stof1 = 6
  ELSE IF (string1=='7') THEN
    stof1 = 7
  ELSE IF (string1=='8') THEN
    stof1 = 8
  ELSE IF (string1=='9') THEN
    stof1 = 9
  END IF
  RETURN
END FUNCTION stof1
!!
!!-------------------------------------------------------------------------
!!
END  MODULE mo_art_fplume_utilities
