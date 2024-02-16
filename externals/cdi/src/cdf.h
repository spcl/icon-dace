#ifndef CDF_H
#define CDF_H

void cdfDebug(int debug);

extern int CDF_Debug;

const char *cdfLibraryVersion(void);

int cdfOpen(const char *filename, const char *mode, int filetype);
int cdf4Open(const char *filename, const char *mode, int *filetype);
void cdfClose(int fileID);

#endif /* CDF_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
