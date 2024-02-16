/*
  CDI PIO C header file

  Include this file in applications to make use of the parallel I/O interfaces
  of CDI.
*/

#ifndef CDIPIO_H_
#define CDIPIO_H_

// clang-format off

//FINT_ON  <--- don't change or remove this line!!!
// Start of fortran interface for the following routines (make_fint.c)

#include <mpi.h>

/* parallel IO IOMode */

#define PIO_NONE 0
#define PIO_MPI 1
#define PIO_WRITER 2
#define PIO_ASYNCH 3
#define PIO_FPGUARD 4
#define PIO_MPI_FW_ORDERED 5
#define PIO_MPI_FW_AT_ALL 6
#define PIO_MPI_FW_AT_REBLOCK 7

#define PIO_MINIOMODE PIO_NONE
#define PIO_MAXIOMODE PIO_MPI_FW_AT_REBLOCK

#define PIO_ROLE_CLIENT 0
#define PIO_ROLE_COLLECTOR 1
#define PIO_ROLE_WRITER 2
#define PIO_ROLE_WRITER_COLLECTOR 3
#define PIO_ROLE_FPGUARD 4

/* parallel IO routines */
#include <yaxt.h>

void pioEndDef(void);
void pioEndTimestepping(void);
void pioFinalize(void);
/* cdiPioNoPostCommSetup: Dummy default function to use as argument to
 * cdiPioConfSetCallBackActions or pioInit
 * if no actions are necessary after I/O servers initialize communication */
void cdiPioNoPostCommSetup(void);
/* pioInit: initialize I/O server processes and communication
 * Deprecated, use cdiPioInit instead! */
MPI_Comm pioInit(MPI_Comm commSuper, int nProcsIO, int IOMode, int *pioNamespace, float partInflate,
                 void (*postCommSetupActions)(void));
/*      cdiPioInit: initialize I/O server processes and communication */
MPI_Comm cdiPioInit(MPI_Comm commSuper, int confResH, int *pioNamespace);
/* pioWriteTimestep: flush data from all client RMA buffers to server */
void pioWriteTimestep(void);
/* cdiPioRDMAProgress: devote some resources to make RMA progress This
 * call is meant for systems where the hardware and/or MPI make
 * insufficient progress when only calling
 * MPI_Win_post/MPI_Win_wait+MPI_Win_start/MPI_Win_get/MPI_Win_complete */
void cdiPioRDMAProgress(void);

/* cdiPioStreamDefDecomposedVlist: collectively define the vlist assigned to a
 * stream together with a fixed decomposition where for each variable
 * an entry of partDesc specifies the local part at the caller, and
 * an entry of conversion specifies if data will be written with
 * streamWriteVar (i.e. as double) with
 * conversion[varID] == CDI_DATATYPE_FLT64 or with streamWriteVarF
 * (i.e. as float) with conversion[varID] == CDI_DATATYPE_FLT32
 */
void cdiPioStreamDefDecomposedVlist(int streamID, int vlistID, const Xt_idxlist partDesc[], const int conversion[]);

/* streamWriteVarPart: Write part of the data making up variable varID
 * of stream streamID.
 *
 * The processes in the communicator returned from cdiPioInit or
 * pioInit must call this routine collectively and data must point to
 * M items, where partDesc is a YAXT index list describing with M
 * indices in 0 to N-1 the data items stored as doubles. N is the
 * number of values per time step in the variable, i.e. the size of
 * the corresponding array passed to streamWriteVar in the serial version.
 * The group of processes collectively calling streamWriteVarPart
 * must provide data for all indices or the behaviour is undefined. */
void streamWriteVarPart(int streamID, int varID, const double *data, int nmiss, Xt_idxlist partDesc);

/* streamWriteVarPartF: Write part of the data making up variable
 * varID of stream streamID.
 *
 * Single-precision version of streamWriteVarPart.
 */
void streamWriteVarPartF(int streamID, int varID, const float *data, int nmiss, Xt_idxlist partDesc);

/* streamWriteScatteredVarPart: Write part of the data making up
 * variable varID of stream streamID.
 *
 * In contrast to streamWriteVarPart, the data is not read from data as one
 * contiguous sequence but instead the numBlocks chunks of length
 * blocklengths[i] and starting displacements[i] each for i in [0,numBlocks)
 */
void streamWriteScatteredVarPart(int streamID, int varID, const double *data, int numBlocks, const int blocklengths[],
                                 const int displacements[], int nmiss, Xt_idxlist partDesc);

/* streamWriteScatteredVarPartF: Write part of the data making up
 * variable varID of stream streamID.
 *
 * Single-precision version of streamWriteScatteredVarPart.
 */
void streamWriteScatteredVarPartF(int streamID, int varID, const float *data, int numBlocks, const int blocklengths[],
                                  const int displacements[], int nmiss, Xt_idxlist partDesc);
/* cdiPioCSRLastN: return role codes appropriate to use the last
   \textit{nProcsIO} tasks as I/O servers */
int cdiPioCSRLastN(MPI_Comm commSuper, int IOMode, int nProcsIO);

/* cdiPioCSRFirstN: return role codes appropriate to use the first
   \textit{nProcsIO} tasks as I/O servers */
int cdiPioCSRFirstN(MPI_Comm commSuper, int IOMode, int nProcsIO);

/* cdiPioCSRBalanced: return role codes appropriate to use \textit{nProcsIO}
 * tasks distributed on evenly spaced ranks as I/O servers */
int cdiPioCSRBalanced(MPI_Comm commSuper, int IOMode, int nProcsIO);

/* cdiPioStr2IOMode: return integer code corresponding to string
 * representation of mode or -1 if no match was found */
int cdiPioStr2IOMode(const char *modeStr);

/* cdiPioStr2IOMode: return string corresponding to integer
 * code of mode or empty string if code is not valid */
const char *cdiPioIOMode2Str(int IOMode);

/* cdiPioConfCreate: create new configuration object and return its handle */
int cdiPioConfCreate(void);

/* cdiPioConfDestroy: delete configuration object */
void cdiPioConfDestroy(int confResH);

/* cdiPioConfSetPartInflate: set partition imbalance attribute of
 * configuration object */
void cdiPioConfSetPartInflate(int confResH, float partInflate);

/* cdiPioConfGetPartInflate: query partition imbalance attribute of
 * configuration object */
float cdiPioConfGetPartInflate(int confResH);

/* cdiPioConfSetIOMode: set IOMode attribute of configuration object */
void cdiPioConfSetIOMode(int confResH, int IOMode);

/* cdiPioConfGetIOMode: query IOMode attribute of configuration object */
int cdiPioConfGetIOMode(int confResH);

/* cdiPioConfSetCSRole: set role attribute of configuration object */
void cdiPioConfSetCSRole(int confResH, int CSRole);

/* cdiPioConfGetCSRole: query role attribute of configuration object */
int cdiPioConfGetCSRole(int confResH);

/* cdiPioConfSetPostCommSetupActions: set function to be called after
 * setup of client/server communications of configuration object.
 * Deprecated: use cdiPioConfSetCallBackActions with
 * trigger == CDIPIO_CALLBACK_POSTCOMMSETUP in new programs! */
void cdiPioConfSetPostCommSetupActions(int confResH, void (*postCommSetupActions)(void));

/* cdiPioConfGetPostCommSetupActions: get function to be called after
 * setup of client/server communications from configuration object.
 * Deprecated: use cdiPioConfGetCallBackActions with
 * trigger == CDIPIO_CALLBACK_POSTCOMMSETUP in new programs. */
void (*cdiPioConfGetPostCommSetupActions(int confResH))(void);

/* CDIPIO_CALLBACK_POSTCOMMSETUP: trigger number of the hook called
 * after communication has been established. This is the same hook
 * previously setup with cdiPioConfSetPostCommSetupActions, takes no
 * argument */
#define CDIPIO_CALLBACK_POSTCOMMSETUP 0
/* CDIPIO_CALLBACK_POSTSTREAMCLOSE: trigger number for callback
 * invoked after each streamClose on the collector side.
 * Accepts the streamID as int parameter, i.e. use INTEGER, VALUE and
 * BIND(C) on Fortran side
 */
#define CDIPIO_CALLBACK_POSTSTREAMCLOSE 1
/* CDIPIO_CALLBACK_POSTWRITEBATCH: trigger number for callback called
 * on server side after the processing for all operations initiated by
 * a client-side pioWriteTimestep have completed */
#define CDIPIO_CALLBACK_POSTWRITEBATCH 2

/* cdiPioConfSetCallBack: set function to be called at
 * indicated trigger of configuration object, action will be cast to
 * the appropriate type as indicated for the respective trigger */
void cdiPioConfSetCallBackActions(int confResH, int trigger, void (*action)(void));

/* cdiPioConfGetCallBack: query function to be called at
 * indicated trigger of configuration object */
void (*cdiPioConfGetCallBackActions(int confResH, int trigger))(void);

/* cdiPioConfSetLargePageAlign should block buffer be aligned to
 * large pages instead of normal pages? */
void cdiPioConfSetLargePageAlign(int confResH, int largePageAlign);

/* cdiPioConfSetLargePageAlign: should block buffer be aligned to
 * large pages instead of normal pages? */
int cdiPioConfGetLargePageAlign(int confResH);

/* cdiPioConfSetRecordAggBufLim: Set limit to pre-encoding step
 * aggregation of data in Mebibyte. Increasing this value trades fewer
 * communication operations between I/O servers for higher packet
 * sizes by transposing more data at once but increases total memory
 * consumed by I/O server processes.  Default size is 128MiB. */
void cdiPioConfSetRecordAggBufLim(int confResH, int lim_mb);

/* cdiPioConfGetRecordAggBufLim: Query size of pre-encoding
 * aggregation buffer in MiB */
int cdiPioConfGetRecordAggBufLim(int confResH);

/* cdiPioConfSetWriteAggBufLim: Set limit for write buffer aggregation
 * (default: 16MiB or value of BUFSIZE environment variable (whichever
 * is larger)). Before writing encoded data records to disk, data
 * up to this size is concatenated. For this reason this must be at
 * least equal in size to the largest GRIB record written.  Increasing
 * this size increases memory consumed by I/O server ranks
 * proportionally and can reduce the number of write operations. This
 * value is rounded to the next (large) page size in many implementations. */
void cdiPioConfSetWriteAggBufLim(int confResH, int lim_mb);

/* cdiPioConfGetWriteAggBufLim: Query the size of write
 * aggregation buffers. */
int cdiPioConfGetWriteAggBufLim(int confResH);

/* cdiPioConfSetAioQueueDepth: Set number of concurrent async I/O
 * requests to create. Depending on implementation, this might
 * increase throughput by increasing the number of concurrent
 * operations but also increases buffer size requirements.
 * (Default value: 4) */
void cdiPioConfSetAioQueueDepth(int confResH, int queue_depth);

/* cdiPioConfGetAioQueueDepth: Query depth of AIO queue. */
int cdiPioConfGetAioQueueDepth(int confResH);

/* cdiPioConfSetMaxPathLen: Set maximal path length allowed in RPC
 * operations. This defaults to 2*PATH_MAX and therefore should be
 * safe in almost any environment. In case of deeply nested directory
 * structures it might be necessary to adjust this value. */
void cdiPioConfSetMaxPathLen(int confResH, int max_path_len);

/* cdiPioConfGetMaxPathLen: Query maximal path length supported in
 * streamOpen operations by some paths of parallel I/O. */
int cdiPioConfGetMaxPathLen(int confResH);

/* cdiPioConfSetRedistCache: set doCache to anything non-zero if data
 * for internal data exchanges is to be cached. This makes sense when
 * the data passed via streamWriteVarPart or streamWriteScatteredVarPart
 * is always decomposed statically using the same partitioning
 * description objects and the sequence of calls to streamWriteVarPart
 * or streamWriteScatteredVarPart for a stream matches the sequence
 * of the previous sequence (divided by pioWriteTimestep) */
void cdiPioConfSetRedistCache(int confResH, int doCache);

/* cdiPioConfSetRedistCache: will data for internal data exchanges
 * be cached? */
int cdiPioConfGetRedistCache(int confResH);

/* cdiPioConfSetXmapCache: set doCache to any non-zero value if data
 * for internal data exchange maps is to be cached. This is helpful
 * when the different variables of a file are decomposed identically
 * and is active by default. Set doCache to 0 if decompositions are
 * dynamic in both variables and time */
void cdiPioConfSetXmapCache(int confResH, int doCache);

/* cdiPioConfGetXmapCache: query xmap caching status */
int cdiPioConfGetXmapCache(int confResH);

/* cdiPioConfSetXmapNew: set method to compute part intersections,
 * defaults to xt_xmap_dist_dir_new */
void cdiPioConfSetXmapNew(int confResH, Xt_xmap (*xmap_new)(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm));

/* cdiPioConfSetXmapNew: get method to compute part intersections */
Xt_xmap (*cdiPioConfGetXmapNew(int confResH))(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist, MPI_Comm comm);
/* cdiPioConfSetStripeConversion: Convert index lists to stripes prior
 * to intersection computation, defaults to true. Only if parts are
 * strictly sections, it can be beneficial to set this to 0, i.e. false */
void cdiPioConfSetStripeConversion(int confResH, int doStripify);

/* cdiPioConfGetStripeConversion: Are index lists of parts converted
 * into stripes before being passed to the xmap constructor? */
int cdiPioConfGetStripeConversion(int confResH);

/* cdiPioConfSetBatchedRMA: (de-)activate batched transfer of data
 * for all streams before any data is written to disk, i.e. if
 * doBatchedRMA is set to 0, data is written as soon as it can be
 * retrieved from client ranks, any other value results in all RMA
 * transfers occurring immediately. doBatchedRMA == 0 implies less
 * memory used on server ranks at the cost of distributing
 * disturbances over a longer time */
void cdiPioConfSetBatchedRMA(int confResH, int doBatchedRMA);

/* cdiPioConfGetBatchedRMA: query if batched RMA is active, see
 * cdiPioConfSetBatchedRMA */
int cdiPioConfGetBatchedRMA(int confResH);

/* cdiPioDistGridCreate: create a grid data structure where the
 * per-coordinate data is distributed over the client tasks
 * chunk_decomposition specifies how to distribute the data of each of
 * the following arrays:
 * x-values, y-values, x-bounds, y-bounds, area, mask and gme mask
 *
 * for 2D grids (all but gridtype == GRID_UNSTRUCTURED), four values specifiy the start and size of each dimension
 * where e.g. chunk_decomposition[1][0] would
 * specify the part size of the y-dimension and chunk_decomposition[0][0] the x dimension start
 * index. For unstructured grids only the x-dimension applies.
 */
int cdiPioDistGridCreate(int gridtype, int size, int xsize, int ysize, int nvertex, const int xy_decomposition_optional[][2],
                         Xt_idxlist partDesc2D, Xt_idxlist partDescX, Xt_idxlist partDescY);

/* cdiPioDistGridEnableIndividualQueries: for the provided gridID, the
 * queries
 * gridInqXval, gridInqYval, gridInqXinc, and gridInqYinc
 *
 * must not be called for a distributed grid before that has been
 * enabled with this routine, also
 *
 * gridInqXvals, gridInqYvals, gridInqXbounds, gridInqYbounds,
 * gridInqArea, gridInqMaskGME and gridInqMask
 *
 * can only be called collectively before this routine was
 * called */
void cdiPioDistGridEnableIndividualQueries(int gridID);

/* cdiPioDistGridDisableIndividualQueries: inverse of
 * cdiPioDistGridEnableIndividualQueries, i.e. after calling this
 * routine, inidividual grid queries must not be made for gridID */
void cdiPioDistGridDisableIndividualQueries(int gridID);

/* cdiPioDistGridIndividualQueriesEnabled: determine if the inquiries
 * listed for cdiPioDistGridEnableIndividualQueries can be issued for
 * gridID */
bool cdiPioDistGridIndividualQueriesEnabled(int gridID);

/* cdiPioInqInterComm: query the intercommunicator of active CDI-PIO
 * namespace */
MPI_Comm cdiPioInqInterComm(void);

// End of fortran interface
//FINT_OFF  <--- don't change or remove this line!!!

// clang-format on

#endif
