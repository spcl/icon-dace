#ifndef ASYNC_WORKER_H
#define ASYNC_WORKER_H

typedef struct AsyncJob AsyncJob;
typedef struct AsyncManager AsyncManager;

// a negative threadCount gives the number of cores that should remain unused by the worker threads, returns an error code
int AsyncWorker_init(AsyncManager **jobManager, int threadCount);

// executes work(data) in a worker thread, must be followed by a call to AsyncWorker_wait()
AsyncJob *AsyncWorker_requestWork(AsyncManager *jobManager, int (*work)(void *data), void *data);

// waits for the async job to finish and returns its result (or some other error code)
int AsyncWorker_wait(AsyncManager *jobManager, AsyncJob *job);

// return the number of workers that are currently idle
int AsyncWorker_availableWorkers(AsyncManager *jobManager);

// waits for all pending jobs to finish, stops all workers, returns a non-zero error code from a pending job if there were any
int AsyncWorker_finalize(AsyncManager *jobManager);

#endif
