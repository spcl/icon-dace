// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

/* Portable CPU-timer (User + Sys); also WALL CLOCK-timer */

// Fortran interface to the following functions is
// implemented in ../src/shared/mo_util_timer.f90

#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>
#include <stdio.h>

#define _HIGH_RESOLUTION_TIMER

int util_cputime(double *user_time, double *system_time) {
    struct tms tbuf;
    const double clock_ticks = (double) sysconf(_SC_CLK_TCK);

    if (times(&tbuf) == -1) return -1;

    *user_time = ((double) tbuf.tms_utime) / clock_ticks;
    *system_time = ((double) tbuf.tms_stime) / clock_ticks;

    return 0;
}

double util_walltime() {
    static double time_init = 0.;
    double time_in_secs;
    struct timeval tbuf;

    if (gettimeofday(&tbuf, NULL) == -1) perror("UTIL_WALLTIME");

    if (time_init == 0.)
        time_init = (double) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0);

    time_in_secs
        = (double) tbuf.tv_sec + (tbuf.tv_usec / 1000000.0) - time_init;

    return time_in_secs;
}

double util_gettimeofday() {
    struct timeval tbuf;

    if (gettimeofday(&tbuf, NULL) == -1) perror("util_gettimeofday");
    return (double) tbuf.tv_sec + tbuf.tv_usec * (1.0 / 1000000.0);
}

/****************************************************************************/

#ifdef _AIX

/*
 * High-Resolution Time
 * pwr4: measurement (util_read_real_time) overhead ~ 0.07 us
 *       conversion  (util_diff_real_time) overhead ~ 0.3 us
 */

#include <sys/systemcfg.h>

static double aix_rt_tconv, aix_rt_fhigh;

void util_init_real_time() {
    double tb_top, tb_bot;

    if (_system_configuration.implementation == POWER_RS2) {
        aix_rt_fhigh = 1.0;
        aix_rt_tconv = 1.0;
    } else {
        /* powerpc family */
        aix_rt_fhigh = 4.294967296;
        tb_top = (double) _system_configuration.Xint;
        tb_bot = (double) _system_configuration.Xfrac;
        aix_rt_tconv = tb_top / tb_bot;
    }
}

void util_get_real_time_size(int *rt_size) {
    /* fortran out:  integer*4:: rt_size(4) */
    *rt_size = (int) sizeof(struct timebasestruct);
}

void util_read_real_time(void *it) {
    /* fortran out:  integer*4:: it(4)
     * raw values - not yet scaled to real time
     */
    read_real_time((struct timebasestruct *) it, TIMEBASE_SZ);
}

void util_diff_real_time(void *it1, void *it2, double *t) {
    /* fortran in:  integer*4:: it1(4), it2(4)
     * fortran out: real*8:: t
     * t is the real time diff between measurements it1 and it2
     */

    struct timebasestruct *tb1, *tb2;

    tb1 = (struct timebasestruct *) it1;
    tb2 = (struct timebasestruct *) it2;

    *t = aix_rt_tconv
         * (aix_rt_fhigh * ((double) (tb2->tb_high - tb1->tb_high))
            + 1.0e-9 * ((double) tb2->tb_low - (double) tb1->tb_low));
}

#elif defined(_HIGH_RESOLUTION_TIMER) && (defined(SX) || defined(ES))

#define CPU_CLOCK                                                         \
    6.25e-10 /* SX-9: 3.2/2 GHz is the increment frequency of the counter \
                register */

long long int irtc(void); /* provided as assembler routine */

void util_init_real_time() {}

void util_get_real_time_size(int *rt_size) { *rt_size = (int) sizeof(double); }

void util_read_real_time(void *it) {
    double *t;
    t = (double *) it;
    *t = CPU_CLOCK * irtc();
}

void util_diff_real_time(void *it1, void *it2, double *t) {
    double *t1, *t2;
    t1 = (double *) it1;
    t2 = (double *) it2;
    *t = *t2 - *t1;
}

#elif defined(_HIGH_RESOLUTION_TIMER) && defined(LINUX)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <asm/msr.h>

static double cpu_clock;

void util_init_real_time() {
    double freq = 0.0;
    FILE *cpuinfo;
    char line[256];

    if ((cpuinfo = fopen("/proc/cpuinfo", "r")) == NULL) {
        fprintf(stderr, "Couldn't open /proc/cpuinfo ...\n");
        exit(-1);
    }

    while (!feof(cpuinfo)) {
        (void) fgets(line, 255, cpuinfo);
        if (strstr(line, "cpu MHz") != NULL) {
            sscanf(line, "cpu MHz         : %lf", &freq);
            break;
        }
    }
    if (freq == 0.0) {
        fprintf(stderr, "Couldn't find cpu MHz in /proc/cpuinfo ...\n");
        exit(-1);
    }

    cpu_clock = freq * 10e6;
}

void util_get_real_time_size(int *rt_size) { *rt_size = (int) sizeof(double); }

void util_read_real_time(void *it) {
    long long llt;
    double *t;
    t = (double *) it;
    rdtscll(llt);
    *t = (double) llt / cpu_freq;
}

void util_diff_real_time(void *it1, void *it2, double *t) {
    double *t1, *t2;
    t1 = (double *) it1;
    t2 = (double *) it2;
    *t = *t2 - *t1;
}

#else

/* fall back to util_walltime */

void util_init_real_time() {}

void util_get_real_time_size(int *rt_size) { *rt_size = (int) sizeof(double); }

void util_read_real_time(void *it) {
    double *t;
    t = (double *) it;
    *t = util_walltime();
}

void util_diff_real_time(void *it1, void *it2, double *t) {
    double *t1, *t2;
    t1 = (double *) it1;
    t2 = (double *) it2;
    *t = *t2 - *t1;
}

#endif
