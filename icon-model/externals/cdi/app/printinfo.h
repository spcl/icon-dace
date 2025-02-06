#ifndef PRINTINFO_H
#define PRINTINFO_H

#include <cdi.h>

void datetime2str(CdiDateTime dt, char *datetimestr, int maxlen);
void date2str(CdiDate date, char *datestr, int maxlen);
void time2str(CdiTime time, char *timestr, int maxlen);

const char *comptype_to_name(int comptype);

void printFiletype(int streamID, int vlistID);
void printGridInfo(int vlistID);
void printZaxisInfo(int vlistID);
void printSubtypeInfo(int vlistID);
void printTimesteps(int streamID, int taxisID, int verbose);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
