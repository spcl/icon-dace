! This file was automatically generated, don't edit!
!
! Fortran interface for CDI library version 2.2.4
!
! Author:
! -------
! Uwe Schulzweida, MPI-MET, Hamburg,   August 2023
!

!
!  Start of fortran interface for the following routines (make_fint.c)
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

!                     cdiPioStreamDefDecomposedVlist
!                                    (INTEGER         streamID,
!                                     INTEGER         vlistID,
!                                     TYPE(XT_IDXLIST)partDesc(*),
!                                     INTEGER         conversion(*))
      EXTERNAL        cdiPioStreamDefDecomposedVlist

!                     streamWriteVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION data,
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteVarPart

!                     streamWriteVarPartF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     REAL            data,
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteVarPartF

!                     streamWriteScatteredVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION data,
!                                     INTEGER         numBlocks,
!                                     INTEGER         blocklengths(*),
!                                     INTEGER         displacements(*),
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteScatteredVarPart

!                     streamWriteScatteredVarPartF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     REAL            data,
!                                     INTEGER         numBlocks,
!                                     INTEGER         blocklengths(*),
!                                     INTEGER         displacements(*),
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteScatteredVarPartF

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

      CHARACTER(80)   cdiPioIOMode2Str
!                                    (INTEGER         IOMode)
      EXTERNAL        cdiPioIOMode2Str

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

      INTEGER    CDIPIO_CALLBACK_POSTCOMMSETUP
      PARAMETER (CDIPIO_CALLBACK_POSTCOMMSETUP =  0)
      INTEGER    CDIPIO_CALLBACK_POSTSTREAMCLOSE
      PARAMETER (CDIPIO_CALLBACK_POSTSTREAMCLOSE =  1)
      INTEGER    CDIPIO_CALLBACK_POSTWRITEBATCH
      PARAMETER (CDIPIO_CALLBACK_POSTWRITEBATCH =  2)
!                     cdiPioConfSetCallBackActions
!                                    (INTEGER         confResH,
!                                     INTEGER         trigger,
!                                     PROCEDURE       action)
      EXTERNAL        cdiPioConfSetCallBackActions

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

!                     cdiPioConfSetRecordAggBufLim
!                                    (INTEGER         confResH,
!                                     INTEGER         lim_mb)
      EXTERNAL        cdiPioConfSetRecordAggBufLim

      INTEGER         cdiPioConfGetRecordAggBufLim
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetRecordAggBufLim

!                     cdiPioConfSetWriteAggBufLim
!                                    (INTEGER         confResH,
!                                     INTEGER         lim_mb)
      EXTERNAL        cdiPioConfSetWriteAggBufLim

      INTEGER         cdiPioConfGetWriteAggBufLim
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetWriteAggBufLim

!                     cdiPioConfSetAioQueueDepth
!                                    (INTEGER         confResH,
!                                     INTEGER         queue_depth)
      EXTERNAL        cdiPioConfSetAioQueueDepth

      INTEGER         cdiPioConfGetAioQueueDepth
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetAioQueueDepth

!                     cdiPioConfSetMaxPathLen
!                                    (INTEGER         confResH,
!                                     INTEGER         max_path_len)
      EXTERNAL        cdiPioConfSetMaxPathLen

      INTEGER         cdiPioConfGetMaxPathLen
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetMaxPathLen

!                     cdiPioConfSetRedistCache
!                                    (INTEGER         confResH,
!                                     INTEGER         doCache)
      EXTERNAL        cdiPioConfSetRedistCache

      INTEGER         cdiPioConfGetRedistCache
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetRedistCache

!                     cdiPioConfSetXmapCache
!                                    (INTEGER         confResH,
!                                     INTEGER         doCache)
      EXTERNAL        cdiPioConfSetXmapCache

      INTEGER         cdiPioConfGetXmapCache
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetXmapCache

!                     cdiPioConfSetStripeConversion
!                                    (INTEGER         confResH,
!                                     INTEGER         doStripify)
      EXTERNAL        cdiPioConfSetStripeConversion

      INTEGER         cdiPioConfGetStripeConversion
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetStripeConversion

!                     cdiPioConfSetBatchedRMA
!                                    (INTEGER         confResH,
!                                     INTEGER         doBatchedRMA)
      EXTERNAL        cdiPioConfSetBatchedRMA

      INTEGER         cdiPioConfGetBatchedRMA
!                                    (INTEGER         confResH)
      EXTERNAL        cdiPioConfGetBatchedRMA

      INTEGER         cdiPioDistGridCreate
!                                    (INTEGER         gridtype,
!                                     INTEGER         size,
!                                     INTEGER         xsize,
!                                     INTEGER         ysize,
!                                     INTEGER         nvertex,
!                                     INTEGER         xy_decomposition_optional(2,*),
!                                     TYPE(XT_IDXLIST)partDesc2D,
!                                     TYPE(XT_IDXLIST)partDescX,
!                                     TYPE(XT_IDXLIST)partDescY)
      EXTERNAL        cdiPioDistGridCreate

!                     cdiPioDistGridEnableIndividualQueries
!                                    (INTEGER         gridID)
      EXTERNAL        cdiPioDistGridEnableIndividualQueries

!                     cdiPioDistGridDisableIndividualQueries
!                                    (INTEGER         gridID)
      EXTERNAL        cdiPioDistGridDisableIndividualQueries

      LOGICAL         cdiPioDistGridIndividualQueriesEnabled
!                                    (INTEGER         gridID)
      EXTERNAL        cdiPioDistGridIndividualQueriesEnabled

      INTEGER         cdiPioInqInterComm
      EXTERNAL        cdiPioInqInterComm

!
!  End of fortran interface
!
