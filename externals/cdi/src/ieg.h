#ifndef _IEG_H
#define _IEG_H

// clang-format off

/* Level Types */
#define  IEG_LTYPE_SURFACE               1
#define  IEG_LTYPE_99                   99
#define  IEG_LTYPE_ISOBARIC            100
#define  IEG_LTYPE_MEANSEA             102
#define  IEG_LTYPE_ALTITUDE            103
#define  IEG_LTYPE_HEIGHT              105
#define  IEG_LTYPE_SIGMA               107
#define  IEG_LTYPE_HYBRID              109
#define  IEG_LTYPE_HYBRID_LAYER        110
#define  IEG_LTYPE_LANDDEPTH           111
#define  IEG_LTYPE_LANDDEPTH_LAYER     112
#define  IEG_LTYPE_SEADEPTH            160

/*
 *  Data representation type (Grid Type) [Table 6]
 */
#define  IEG_GTYPE_LATLON             0  /*  latitude/longitude                       */
#define  IEG_GTYPE_LATLON_ROT        10  /*  rotated latitude/longitude               */

#define  IEG_P_CodeTable(x)   (x[ 5])  /*  Version number of code table                 */
#define  IEG_P_Parameter(x)   (x[ 6])  /*  Parameter indicator                          */
#define  IEG_P_LevelType(x)   (x[ 7])  /*  Type of level indicator                      */
#define  IEG_P_Level1(x)      (x[ 8])  /*  Level 1                                      */
#define  IEG_P_Level2(x)      (x[ 9])  /*  Level 2                                      */
#define  IEG_P_Year(x)        (x[10])  /*  Year of century (YY)                         */
#define  IEG_P_Month(x)       (x[11])  /*  Month (MM)                                   */
#define  IEG_P_Day(x)         (x[12])  /*  Day (DD)                                     */
#define  IEG_P_Hour(x)        (x[13])  /*  Hour (HH)                                    */
#define  IEG_P_Minute(x)      (x[14])  /*  Minute (MM)                                  */

/*
 *  Macros for the grid definition section ( Section 2 )
 */
#define  IEG_G_Size(x)        (x[ 0])
#define  IEG_G_NumVCP(x)      (x[3] == 10 ? (x[0]-42)/4 : (x[0]-32)/4)
#define  IEG_G_GridType(x)    (x[ 3])  /*  Data representation type */
#define  IEG_G_NumLon(x)      (x[ 4])  /*  Number of points along a parallel (Ni)       */
#define  IEG_G_NumLat(x)      (x[ 5])  /*  Number of points along a meridian (Nj)       */
#define  IEG_G_FirstLat(x)    (x[ 6])  /*  Latitude of the first grid point             */
#define  IEG_G_FirstLon(x)    (x[ 7])  /*  Longitude of the first grid point            */
#define  IEG_G_ResFlag(x)     (x[ 8])  /*  Resolution flag: 128 regular grid            */
#define  IEG_G_LastLat(x)     (x[ 9])  /*  Latitude of the last grid point              */
#define  IEG_G_LastLon(x)     (x[10])  /*  Longitude of the last grid point             */
#define  IEG_G_LonIncr(x)     (x[11])  /*  i direction increment                        */
#define  IEG_G_LatIncr(x)     (x[12])  /*  j direction increment                        */
#define  IEG_G_ScanFlag(x)    (x[13])
#define  IEG_G_LatSP(x)       (x[16])  /*  Latitude of the southern pole of rotation    */
#define  IEG_G_LonSP(x)       (x[17])  /*  Longitude of the southern pole of rotation   */
#define  IEG_G_ResFac(x)      (x[18])  /*  Resolution factor                            */

// clang-format on

typedef struct
{
  int checked;
  int byteswap;
  int dprec; /* data   precision */
  int ipdb[37];
  double refval;
  int igdb[22];
  double vct[100];
  size_t datasize;
  size_t buffersize;
  void *buffer;
} iegrec_t;

const char *iegLibraryVersion(void);

void iegDebug(int debug);
int iegCheckFiletype(int fileID, int *swap);

void *iegNew(void);
void iegDelete(void *ieg);
void iegInitMem(void *ieg);

int iegRead(int fileID, void *ieg);
int iegWrite(int fileID, void *ieg);

void iegCopyMeta(void *dieg, void *sieg);
int iegInqHeader(void *ieg, int *header);
int iegInqDataSP(void *ieg, float *data);
int iegInqDataDP(void *ieg, double *data);

int iegDefHeader(void *ieg, const int *header);
int iegDefDataSP(void *ieg, const float *data);
int iegDefDataDP(void *ieg, const double *data);

#endif /* _IEG_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
