#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#include "async_worker.h"

#include "cdi.h"
#include "error.h"

#include <stdbool.h>

#ifdef HAVE_LIBPTHREAD
#ifdef __APPLE__
#include <dispatch/dispatch.h>
#else
#include <errno.h>
#include <semaphore.h>
#endif

typedef struct sema
{
#ifdef __APPLE__
  dispatch_semaphore_t sem;
#else
  sem_t sem;
#endif
} sema_t;
#endif

struct AsyncJob
{
  bool inUse;
#ifdef HAVE_LIBPTHREAD
  sema_t request, completion;
#endif
  int (*work)(void *data);
  void *data;
  int result;
};

struct AsyncManager
{
  int workerCount, idleWorkerCount;
  AsyncJob *communicators;
};

#ifdef HAVE_LIBPTHREAD
static inline int
sema_init(sema_t *s, int pshared, uint32_t value)
{
  int status = 0;
#ifdef __APPLE__
  dispatch_semaphore_t *sem = &s->sem;

  (void) pshared;
  *sem = dispatch_semaphore_create(value);
#else
  status = sem_init(&s->sem, pshared, value);
#endif
  return status;
}

static inline int
sema_wait(sema_t *s)
{
#ifdef __APPLE__
  dispatch_semaphore_wait(s->sem, DISPATCH_TIME_FOREVER);
#else
  int r;

  do
    {
      r = sem_wait(&s->sem);
    }
  while (r == -1 && errno == EINTR);
#endif
  return 0;
}

static inline int
sema_post(sema_t *s)
{
#ifdef __APPLE__
  dispatch_semaphore_signal(s->sem);
#else
  sem_post(&s->sem);
#endif
  return 0;
}

static void *
workerMain(void *arg)
{
  AsyncJob *communicator = (AsyncJob *) arg;

  while (true)
    {
      while (sema_wait(&communicator->request))
        ;
      if (communicator->work)
        {
          communicator->result = communicator->work(communicator->data);
          if (sema_post(&communicator->completion)) xabort("sema_post() failed");
        }
      else
        {
          if (sema_post(&communicator->completion)) xabort("sema_post() failed");
          break;
        }
    }

  return NULL;
}

static void
startWorker(AsyncJob *communicator)
{
  communicator->inUse = false;
  communicator->work = NULL;
  communicator->data = NULL;
  communicator->result = 0;
  if (sema_init(&communicator->request, 0, 0)) xabort("sema_init() failed");
  if (sema_init(&communicator->completion, 0, 0)) xabort("sema_init() failed");

  pthread_t worker;
  if (pthread_create(&worker, NULL, workerMain, communicator)) xabort("pthread_create() failed");
  if (pthread_detach(worker)) xabort("pthread_detach() failed");
}
#endif

int
AsyncWorker_init(AsyncManager **jobManager, int threadCount)
{
  if (threadCount <= 0)
    {
      xabort("CPU core count discovery not implemented yet");
      return CDI_EINVAL;  // TODO: discover CPU core count, and set threadCount to a sensible positive value
    }

  if (*jobManager) return CDI_NOERR;

#ifdef HAVE_LIBPTHREAD
  *jobManager = (AsyncManager *) malloc(sizeof(AsyncManager));
  if (!*jobManager) return CDI_ESYSTEM;
  (*jobManager)->workerCount = threadCount;
  (*jobManager)->communicators = (AsyncJob *) malloc(threadCount * sizeof(AsyncJob));
  if (!(*jobManager)->communicators) xabort("memory allocation failure");

  for (int i = 0; i < threadCount; i++) startWorker(&((*jobManager)->communicators[i]));
  (*jobManager)->idleWorkerCount = threadCount;
#else

  Error("pthread support not compiled in!");
#endif

  return CDI_NOERR;
}

AsyncJob *
AsyncWorker_requestWork(AsyncManager *jobManager, int (*work)(void *data), void *data)
{
  if (!jobManager) xabort("AsyncWorker_requestWork() called without calling AsyncWorker_init() first");
  if (!work)
    xabort("AsyncWorker_requestWork() called without a valid function pointer");  // need to catch this condition to stop users from
                                                                                  // terminating our worker threads

  // find an unused worker
  if (!jobManager->idleWorkerCount) return NULL;

  AsyncJob *worker = NULL;
  for (int i = 0; i < jobManager->workerCount; i++)
    {
      if (!jobManager->communicators[i].inUse)
        {
          worker = &jobManager->communicators[i];
          break;
        }
    }
  if (!worker) xabort("internal error: idleWorkerCount is not in sync with the worker states, please report this bug");

  // pass the request to that worker
  jobManager->idleWorkerCount--;
  worker->inUse = true;
  worker->work = work;
  worker->data = data;
  worker->result = 0;
#ifdef HAVE_LIBPTHREAD
  if (sema_post(&worker->request)) xabort("sema_post() failed");
#endif
  return worker;
}

int
AsyncWorker_wait(AsyncManager *jobManager, AsyncJob *job)
{
  if (!jobManager) xabort("AsyncWorker_wait() called without calling AsyncWorker_init() first");
  if (job < jobManager->communicators) return CDI_EINVAL;
  if (job >= jobManager->communicators + jobManager->workerCount) return CDI_EINVAL;
  if (!job->inUse) return CDI_EINVAL;

#ifdef HAVE_LIBPTHREAD
  while (sema_wait(&job->completion))
    ;
#endif
  int result = job->result;

  // reset the communicator
  job->work = NULL;
  job->data = NULL;
  job->result = 0;
  job->inUse = false;
  jobManager->idleWorkerCount++;

  return result;
}

int
AsyncWorker_availableWorkers(AsyncManager *jobManager)
{
  if (!jobManager) return 0;
  return jobManager->idleWorkerCount;
}

int
AsyncWorker_finalize(AsyncManager *jobManager)
{
  int result = CDI_NOERR;
  if (!jobManager) return CDI_NOERR;

  for (int i = 0; i < jobManager->workerCount; i++)
    {
      AsyncJob *curWorker = &jobManager->communicators[i];

      // finish any pending job
      if (curWorker->inUse)
        {
          AsyncWorker_wait(jobManager, curWorker);
          if (curWorker->result) result = curWorker->result;
        }

      // send the teardown signal
      curWorker->inUse = true;
      curWorker->work = NULL;
      curWorker->data = NULL;
      curWorker->result = 0;
#ifdef HAVE_LIBPTHREAD
      if (sema_post(&curWorker->request)) xabort("sema_post() failed");
#endif
      // wait for the worker to exit
      AsyncWorker_wait(jobManager, curWorker);
    }

  free(jobManager->communicators);
  free(jobManager);

  return result;
}
