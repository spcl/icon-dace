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

#define _GNU_SOURCE
#ifdef HAVE_LINK_H
#include <link.h>
#include <mcheck.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <pwd.h>
#include <sys/types.h>
#include <unistd.h>

#include <sys/utsname.h>
#include <sys/resource.h>

#ifdef _AIX
#include <sys/systemcfg.h>
#endif

/* function implementations */

void util_user_name(char *name, int *actual_len) {
    struct passwd *current;

#if defined(__XT3__)
    current = NULL;
#else
    current = getpwuid(getuid());
#endif
    if (current == NULL) {
#if defined(__XT3__)
        strcpy(name, "user name not available");
#else
        strcpy(name, "unknown user name");
#endif
    } else {
        if (strlen(current->pw_name) == 0) {
            strcpy(name, "unknown user name");
        } else {
            strcpy(name, current->pw_gecos);
            strcat(name, " (");
            strcat(name, current->pw_name);
            strcat(name, ")");
        }
    }

    *actual_len = strlen(name);
}

void util_os_system(char *name, int *actual_len) {
    struct utsname system_info;
    int iret;
#ifdef _AIX
    int imp;
#endif

    iret = uname(&system_info);

    /* Cray doesn't handle operating system information proper */

#ifdef _UNICOSMP
    strcpy(name, "UNICOS");
#else
    strcpy(name, system_info.sysname);
#endif
    strcat(name, " ");
#ifdef _AIX
    strcat(name, system_info.version);
    strcat(name, ".");
    strcat(name, system_info.release);
#else
    strcat(name, system_info.release);
#endif
    strcat(name, " ");
#ifdef _AIX
    imp = _system_configuration.implementation;
    if (imp == POWER_4)
        strcat(name, "Power4");
    else if (imp == POWER_5)
        strcat(name, "Power5");
    else if (imp == POWER_6)
        strcat(name, "Power6");
#ifdef _AIX61
    else if (imp == POWER_7)
        strcat(name, "Power7");
#endif
    else {
        strcat(name, "unknown");
    }
#else
    strcat(name, system_info.machine);
#endif

    *actual_len = strlen(name);
}

void util_node_name(char *name, int *actual_len) {
    struct utsname host_info;
    int iret;
    char *fqdn;

    iret = uname(&host_info);
    if (iret != EFAULT) {
        strcpy(name, host_info.nodename);

    } else {
        strcpy(name, "unknown");
    }

    /* if FQDN is returned remove domainname */
    if ((fqdn = strchr(name, '.')) != NULL) {
        *fqdn = '\0';
    }

    *actual_len = strlen(name);

    return;
}


/* Get the maximum resident set size used, in kilobytes. That is, the maximum
 * number of kilobytes of physical memory that processes used
 * simultaneously.
 */
void util_get_maxrss(int *maxrss) {
#if defined HAVE_GETRUSAGE
    struct rusage usage;

    /* get resource usage */
    usage.ru_maxrss = 0;
    if (getrusage(RUSAGE_SELF, &usage) == -1)
        *maxrss = 0; /* Error */
    else
        *maxrss = (int) (usage.ru_maxrss / 1024);

#else
    /* do nothing */
#endif
}


void util_compiler_release(char *release_str, int *rstr_len) {
#ifdef _CRAYC
    strcpy(release_str, _RELEASE_STRING);
#else
    strcpy(release_str, "unknown");
#endif
    *rstr_len = strlen(release_str);
    return;
}

#ifdef HAVE_LINK_H
static int dump_dl(struct dl_phdr_info *info, size_t size, void *data) {
    FILE *f = (FILE *) data;

    fprintf(f, "%s: dlpi_addr = 0x%016lx\n", info->dlpi_name,
            (unsigned long) info->dlpi_addr);

    for (ElfW(Half) i = 0; i < info->dlpi_phnum; ++i) {
        const ElfW(Phdr) *phdr = info->dlpi_phdr + i;
        if (phdr->p_type != PT_LOAD) continue;

        fprintf(f, "\tLOAD: 0x%016lx 0x%016lx 0x%lx\n",
                (unsigned long) phdr->p_vaddr, (unsigned long) phdr->p_paddr,
                (unsigned long) phdr->p_memsz);
    }

    return 0;
}

static void dump_dls(void) {
    char fname[256];
    FILE *f;

    snprintf(fname, sizeof(fname), "dl.map.%lu", (unsigned long) getpid());

    f = fopen(fname, "w");

    if (f == NULL) return;

    dl_iterate_phdr(dump_dl, f);

    fclose(f);

    return;
}


void util_set_mtrace(void) {
    char malloc_trace[256];

    snprintf(malloc_trace, sizeof(malloc_trace), "m.trace.%lld",
             (long long) getpid());

    setenv("MALLOC_TRACE", malloc_trace, 1);

    printf("MALLOC_TRACE:%s\n", getenv("MALLOC_TRACE"));


    dump_dls();

    mtrace();
}


void util_unset_mtrace(void) { muntrace(); }
#endif


void util_c_getpid(long int *pid) {
    *pid = (long int) getpid();
    return;
}
