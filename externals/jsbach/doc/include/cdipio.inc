! This file was automatically generated, don't edit!
!
! Fortran interface for CDI library version 1.8.1
!
! Author:
! -------
! Uwe Schulzweida, MPI-MET, Hamburg,   April 2017
!

!
!  parallel IO IOMode
!
      INTEGER    PIO_NONE
      PARAMETER (PIO_NONE               =  0)
      INTEGER    PIO_MPI
      PARAMETER (PIO_MPI                =  1)
      INTEGER    PIO_WRITER
      PARAMETER (PIO_WRITER             =  2)
      INTEGER    PIO_ASYNCH
      PARAMETER (PIO_ASYNCH             =  3)
      INTEGER    PIO_FPGUARD
      PARAMETER (PIO_FPGUARD            =  4)
      INTEGER    PIO_MPI_FW_ORDERED
      PARAMETER (PIO_MPI_FW_ORDERED     =  5)
      INTEGER    PIO_MPI_FW_AT_ALL
      PARAMETER (PIO_MPI_FW_AT_ALL      =  6)
      INTEGER    PIO_MPI_FW_AT_REBLOCK
      PARAMETER (PIO_MPI_FW_AT_REBLOCK  =  7)
      INTEGER    PIO_MINIOMODE
      PARAMETER (PIO_MINIOMODE          = PIO_NONE)
      INTEGER    PIO_MAXIOMODE
      PARAMETER (PIO_MAXIOMODE          = PIO_MPI_FW_AT_REBLOCK)
      INTEGER    PIO_ROLE_CLIENT
      PARAMETER (PIO_ROLE_CLIENT        =  0)
      INTEGER    PIO_ROLE_COLLECTOR
      PARAMETER (PIO_ROLE_COLLECTOR     =  1)
      INTEGER    PIO_ROLE_WRITER
      PARAMETER (PIO_ROLE_WRITER        =  2)
      INTEGER    PIO_ROLE_WRITER_COLLECTOR
      PARAMETER (PIO_ROLE_WRITER_COLLECTOR =  3)
      INTEGER    PIO_ROLE_FPGUARD
      PARAMETER (PIO_ROLE_FPGUARD       =  4)
!
!  parallel IO routines
!
!                     pioEndDef
      EXTERNAL        pioEndDef

!                     pioEndTimestepping
      EXTERNAL        pioEndTimestepping

!                     pioFinalize
      EXTERNAL        pioFinalize

!                     cdiPioNoPostCommSetup
      EXTERNAL        cdiPioNoPostCommSetup

      INTEGER         pioInit
!                                    (INTEGER         commSuper,
!                                     INTEGER         nProcsIO,
!                                     INTEGER         IOMode,
!                                     INTEGER         pioNamespace,
!                                     REAL            partInflate,
!                                     PROCEDURE       postCommSetupActions)
      EXTERNAL        pioInit

      INTEGER         cdiPioInit
!                                    (INTEGER         commSuper,
!                                     INTEGER         confResH,
!                                     INTEGER         pioNamespace)
      EXTERNAL        cdiPioInit

!                     pioWriteTimestep
      EXTERNAL        pioWriteTimestep

!                     cdiPioRDMAProgress
      EXTERNAL        cdiPioRDMAProgress

!                     streamWriteVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     CHOICE          data,
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteVarPart

!                     streamWriteScatteredVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     CHOICE          data,
!                                     INTEGER         numBlocks,
!                                     INTEGER         blocklengths(*),
!                                     INTEGER         displacements(*),
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteScatteredVarPart

      INTEGER         cdiPioCSRLastN
!                                    (INTEGER         commSuper,
!                                     INTEGER         IOMode,
!                                     INTEGER         nProcsIO)
      EXTERNAL        cdiPioCSRLastN

      INTEGER         cdiPioCSRFirstN
!                                    (INTEGER         commSuper,
!                                     INTEGER         IOMode,
!                                     INTEGER         nProcsIO)
      EXTERNAL        cdiPioCSRFirstN

      INTEGER         cdiPioCSRBalanced
!                                    (INTEGER         commSuper,
!                                     INTEGER         IOMode,
!                                     INTEGER         nProcsIO)
      EXTERNAL        cdiPioCSRBalanced

      INTEGER         cdiPioStr2IOMode
!                                    (CHARACTER*(*)   modeStr)
      EXTERNAL        cdiPioStr2IOMode

      INTEGER         cdiPioConfCreate
      EXTERNAL        cdiPioConfCreate

!                     cdiPioConfDestroy
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfDestroy

!                     cdiPioConfSetPartInflate
!                                    (INTEGER         confResH,
!                                     REAL            partInflate)
      EXTERNAL        cdiPioConfSetPartInflate

      REAL            cdiPioConfGetPartInflate
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetPartInflate

!                     cdiPioConfSetIOMode
!                                    (INTEGER         confResH,
!                                     INTEGER         IOMode)
      EXTERNAL        cdiPioConfSetIOMode

      INTEGER         cdiPioConfGetIOMode
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetIOMode

!                     cdiPioConfSetCSRole
!                                    (INTEGER         confResH,
!                                     INTEGER         CSRole)
      EXTERNAL        cdiPioConfSetCSRole

      INTEGER         cdiPioConfGetCSRole
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetCSRole

!                     cdiPioConfSetPostCommSetupActions
!                                    (INTEGER         confResH,
!                                     PROCEDURE       postCommSetupActions)
      EXTERNAL        cdiPioConfSetPostCommSetupActions

!
!  cdiPioConfSetLargePageAlign should block buffer be aligned to
!  large pages instead of normal pages?
!
!                     cdiPioConfSetLargePageAlign
!                                    (INTEGER         confResH,
!                                     INTEGER         largePageAlign)
      EXTERNAL        cdiPioConfSetLargePageAlign

      INTEGER         cdiPioConfGetLargePageAlign
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetLargePageAlign

!                     cdiPioConfSetRedistCache
!                                    (INTEGER         confResH,
!                                     INTEGER         doCache)
      EXTERNAL        cdiPioConfSetRedistCache

      INTEGER         cdiPioConfGetRedistCache
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetRedistCache

!
!  convert index lists to stripes prior to intersection computation,
!  defaults to true
!
!                     cdiPioConfSetStripeConversion
!                                    (INTEGER         confResH,
!                                     INTEGER         doStripify)
      EXTERNAL        cdiPioConfSetStripeConversion

!
!  are index lists of parts converted stripes before being passed
!  to the xmap constructor?
!
      INTEGER         cdiPioConfGetStripeConversion
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetStripeConversion

