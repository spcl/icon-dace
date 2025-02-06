#ifndef _SERVICE_H
#define _SERVICE_H

typedef struct
{
  int checked;
  int byteswap;
  int header[8];
  int hprec; /* header precision */
  int dprec; /* data   precision */
  size_t datasize;
  size_t buffersize;
  void *buffer;
} srvrec_t;

const char *srvLibraryVersion(void);

void srvDebug(int debug);

int srvCheckFiletype(int fileID, int *swap);

void *srvNew(void);
void srvDelete(void *srv);

int srvRead(int fileID, void *srv);
void srvWrite(int fileID, void *srv);

int srvInqHeader(void *srv, int *header);
int srvInqDataSP(void *srv, float *data);
int srvInqDataDP(void *srv, double *data);

int srvDefHeader(void *srv, const int *header);
int srvDefDataSP(void *srv, const float *data);
int srvDefDataDP(void *srv, const double *data);

#endif /* _SERVICE_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
