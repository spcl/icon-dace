#ifndef CDI_DATETIME_H
#define CDI_DATETIME_H

#include <stdbool.h>
#include <stdint.h>

// clang-format off

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int   year;      // year of date
  short month;     // month of date
  short day;       // day of date
} CdiDate;

typedef struct
{
  short hour;      // hour part of time
  short minute;	   // minute part of time
  short second;	   // second part of time
  short ms;        // milli-second part of time. 0<=ms<=999
} CdiTime;

typedef struct
{
  CdiDate date;    // date elements
  CdiTime time;    // time elements
} CdiDateTime;

CdiDateTime cdiDateTime_set(int64_t date, int time);
CdiDate cdiDate_set(int64_t date);
CdiTime cdiTime_set(int time);
int64_t cdiDate_get(CdiDate cdiDate);
int cdiTime_get(CdiTime cdiTime);

CdiDate cdiDate_encode(int year, int month, int day);
void cdiDate_decode(CdiDate cdiDate, int *year, int *month, int *day);
CdiTime cdiTime_encode(int hour, int minute, int second, int ms);
void cdiTime_decode(CdiTime cdiTime, int *hour, int *minute, int *second, int *ms);

void cdiDate_init(CdiDate *cdiDate);
void cdiTime_init(CdiTime *cdiTime);
void cdiDateTime_init(CdiDateTime *cdiDateTime);

bool cdiDate_isEQ(CdiDate cdiDate1, CdiDate cdiDate2);
bool cdiTime_isEQ(CdiTime cdiTime1, CdiTime cdiTime2);
bool cdiDateTime_isEQ(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2);
bool cdiDateTime_isNE(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2);
bool cdiDateTime_isLT(CdiDateTime cdiDateTime1, CdiDateTime cdiDateTime2);
bool cdiDateTime_isNull(CdiDateTime cdiDateTime);

const char *CdiDateTime_string(CdiDateTime cdiDateTime);

#ifdef __cplusplus
}
#endif

// clang-format on

#endif /* CDI_DATETIME_H */
