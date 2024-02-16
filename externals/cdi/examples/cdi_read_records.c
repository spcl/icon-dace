#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/cdi.h"

#define PARIO 1

#ifdef PARIO
#define HAVE_LIBPTHREAD 1
#endif

#if defined(HAVE_LIBPTHREAD)
#include <pthread.h>
#endif

typedef struct
{
  int streamID;
  int *varID, *levelID;
  size_t *nmiss;
  double *array;
} read_arg_t;

typedef struct
{
  int varID, levelID;
  size_t nmiss;
  double *array;
  int array_size;
  int recID, nrecs;
  read_arg_t read_arg;
  void *iothread;
} par_io_t;

typedef struct work_st
{
  void (*routine)(void *);
  void *arg;
} work_t;

typedef struct
{
  int varID, levelID;
  size_t nmiss;
  double *array;
  int array_size;
  int recID, nrecs;
  int shutdown;
  int used;
  work_t *work;
#if defined(HAVE_LIBPTHREAD)
  pthread_t thrID;
  pthread_attr_t attr;
  pthread_mutex_t lock;
  pthread_cond_t not_empty;
  pthread_cond_t empty;
#endif
} iothread_t;

void
readRecord(void *arg)
{
  read_arg_t *read_arg = (read_arg_t *) arg;
  int streamID;
  int *varID, *levelID;
  size_t *nmiss;
  double *array;

  streamID = read_arg->streamID;
  varID = read_arg->varID;
  levelID = read_arg->levelID;
  nmiss = read_arg->nmiss;
  array = read_arg->array;

  streamInqRecord(streamID, varID, levelID);
  streamReadRecord(streamID, array, nmiss);
}

#if defined(HAVE_LIBPTHREAD)
/* This function is the work function of the thread */
void *
do_work(void *p)
{
  iothread_t *iothread = (iothread_t *) p;

  while (1)
    {
      pthread_mutex_lock(&(iothread->lock));

      while (iothread->used == 0)
        {
          if (iothread->shutdown)
            {
              pthread_mutex_unlock(&(iothread->lock));
              //    pthread_exit(NULL);
              return (NULL);
            }
          // wait until the condition says its no emtpy and give up the lock.
          pthread_mutex_unlock(&(iothread->lock));
          pthread_cond_wait(&(iothread->not_empty), &(iothread->lock));

          if (iothread->shutdown)
            {
              pthread_mutex_unlock(&(iothread->lock));
              //    pthread_exit(NULL);
              return (NULL);
            }
        }

      (iothread->work->routine)(iothread->work->arg);  // actually do work.
      free(iothread->work);

      iothread->used = 0;

      // now signal that its empty.
      pthread_cond_signal(&(iothread->empty));

      pthread_mutex_unlock(&(iothread->lock));
    }
}
#endif

#if defined(HAVE_LIBPTHREAD)
iothread_t *
create_iothread()
{
  iothread_t *iothread = NULL;

  iothread = (iothread_t *) malloc(sizeof(iothread_t));
  if (iothread == NULL)
    {
      fprintf(stderr, "Out of memory creating a new iothread!\n");
      return (NULL);
    }

  // initialize mutex and condition variables.
  if (pthread_mutex_init(&iothread->lock, NULL))
    {
      fprintf(stderr, "Mutex initiation error!\n");
      return (NULL);
    }

  if (pthread_cond_init(&(iothread->empty), NULL))
    {
      fprintf(stderr, "CV initiation error!\n");
      return (NULL);
    }

  if (pthread_cond_init(&(iothread->not_empty), NULL))
    {
      fprintf(stderr, "CV initiation error!\n");
      return (NULL);
    }

  pthread_attr_init(&(iothread->attr));
  pthread_attr_setdetachstate(&(iothread->attr), PTHREAD_CREATE_JOINABLE);

  iothread->shutdown = 0;
  iothread->used = 0;

  // make thread
  if (pthread_create(&(iothread->thrID), &(iothread->attr), do_work, (void *) iothread))
    {
      fprintf(stderr, "Thread initiation error!\n");
      return NULL;
    }

  return (iothread);
}
#endif

#if defined(HAVE_LIBPTHREAD)
void
check_iothread(void *p)
{
  iothread_t *iothread = (iothread_t *) p;

  pthread_mutex_lock(&(iothread->lock));

  while (iothread->used == 1)
    {
      // wait until the condition says its emtpy and give up the lock.
      pthread_mutex_unlock(&(iothread->lock));  // get the lock.
      pthread_cond_wait(&(iothread->empty), &(iothread->lock));
    }

  pthread_mutex_unlock(&(iothread->lock));
}
#endif

#if defined(HAVE_LIBPTHREAD)
typedef void (*dispatch_fn)(void *);

void
dispatch(void *p, dispatch_fn dispatch_to_here, void *arg)
{
  iothread_t *iothread = (iothread_t *) p;
  work_t *work;

  // make a work element.
  work = (work_t *) malloc(sizeof(work_t));
  if (work == NULL)
    {
      fprintf(stderr, "Out of memory creating a work struct!\n");
      return;
    }

  work->routine = dispatch_to_here;
  work->arg = arg;

  pthread_mutex_lock(&(iothread->lock));

  if (iothread->used == 1) fprintf(stderr, "dispatch: Internal synchronization problem!\n");

  iothread->work = work;
  iothread->used = 1;
  pthread_cond_signal(&(iothread->not_empty));

  pthread_mutex_unlock(&(iothread->lock));
}
#endif

#if defined(HAVE_LIBPTHREAD)
void
destroy_iothread(void *p)
{
  iothread_t *iothread = (iothread_t *) p;
  int status;

  pthread_mutex_lock(&(iothread->lock));

  if (iothread->used == 1) fprintf(stderr, "destroy_iothread: Internal synchronization problem!\n");

  iothread->shutdown = 1;
  pthread_cond_signal(&(iothread->not_empty));

  pthread_mutex_unlock(&(iothread->lock));

  status = pthread_join(iothread->thrID, NULL);
  if (status > 0)
    {
      fprintf(stderr, "pthread_join error!\n");
      return;
    }

  pthread_mutex_destroy(&(iothread->lock));
  pthread_cond_destroy(&(iothread->empty));
  pthread_cond_destroy(&(iothread->not_empty));
  pthread_attr_destroy(&(iothread->attr));

  free(iothread);
}
#endif

void
stream_read_record_par(int streamID, int *varID, int *levelID, double *array, size_t *nmiss, par_io_t *parIO)
{
  int lpario = 0;
  int recID = 0, nrecs = 0;

#if defined(HAVE_LIBPTHREAD)
  if (parIO)
    {
      lpario = 1;
      recID = parIO->recID;
      nrecs = parIO->nrecs;
    }
#endif

  if (recID == 0 || lpario == 0)
    {
      read_arg_t read_arg;
      read_arg.streamID = streamID;
      read_arg.varID = varID;
      read_arg.levelID = levelID;
      read_arg.nmiss = nmiss;
      read_arg.array = array;

      readRecord(&read_arg);
    }
#if defined(HAVE_LIBPTHREAD)
  else
    {
      check_iothread(parIO->iothread);

      *varID = parIO->varID;
      *levelID = parIO->levelID;
      *nmiss = parIO->nmiss;

      memcpy(array, parIO->array, parIO->array_size * sizeof(double));
    }

  if (lpario && nrecs > 1)
    {
      if ((recID + 1) < nrecs)
        {
          read_arg_t *read_arg = &(parIO->read_arg);

          if (recID == 0) parIO->iothread = create_iothread();

          read_arg->streamID = streamID;
          read_arg->varID = &parIO->varID;
          read_arg->levelID = &parIO->levelID;
          read_arg->nmiss = &parIO->nmiss;
          read_arg->array = parIO->array;

          dispatch(parIO->iothread, readRecord, read_arg);
        }
      else
        {
          destroy_iothread(parIO->iothread);
        }
    }
#endif
}

void
stream_read_record(int streamID, int *varID, int *levelID, double *data, size_t *nmiss)
{
  streamInqRecord(streamID, varID, levelID);
  streamReadRecord(streamID, data, nmiss);
}

int
main(int argc, char *argv[])
{
  int taxisID, vlistID, varID, levelID, streamID, tsID;
  size_t nmiss;
  int vdate, vtime;
  int nrecs, recID, code;
  int gridsize, i;
  double *data;
  double fmin, fmax, fmean;
  par_io_t parIO;
  char *fname = NULL;

  if (argc != 2)
    {
      fprintf(stderr, "usage: %s filename\n", argv[0]);
      return (-1);
    }

  fname = argv[1];

  /* Open the dataset */
  streamID = streamOpenRead(fname);
  if (streamID < 0)
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return (1);
    }

  /* Get the variable list of the dataset */
  vlistID = streamInqVlist(streamID);

  /* Get the Time axis from the variable list */
  taxisID = vlistInqTaxis(vlistID);

  gridsize = vlistGridsizeMax(vlistID);
  data = (double *) malloc(gridsize * sizeof(double));
#ifdef PARIO
  parIO.array = (double *) malloc(gridsize * sizeof(double));
  parIO.array_size = gridsize;
#endif

  /* Loop over all time steps */
  tsID = 0;
  while ((nrecs = streamInqTimestep(streamID, tsID)))
    {
      /* Get the verification date and time */
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);

      /* Read all records */
      for (recID = 0; recID < nrecs; recID++)
        {
          parIO.recID = recID;
          parIO.nrecs = nrecs;
#ifdef PARIO
          stream_read_record_par(streamID, &varID, &levelID, data, &nmiss, &parIO);
#else
          stream_read_record(streamID, &varID, &levelID, data, &nmiss);
#endif
          code = vlistInqVarCode(vlistID, varID);
          gridsize = gridInqSize(vlistInqVarGrid(vlistID, varID));
          fmin = 1.e33;
          fmax = -1.e33;
          fmean = 0;
          for (int j = 0; j < 2; ++j)
            for (i = 0; i < gridsize; ++i)
              {
                if (data[i] < fmin) fmin = data[i];
                if (data[i] > fmax) fmax = data[i];
                fmean += data[i];
              }
          fmean /= gridsize;
          fprintf(stdout, "%3d %3d %3d %d %g %g %g\n", code, varID, levelID, gridsize, fmin, fmean, fmax);
        }

      tsID++;
    }

  /* Close the input stream */
  streamClose(streamID);

  free(data);

  return 0;
}
