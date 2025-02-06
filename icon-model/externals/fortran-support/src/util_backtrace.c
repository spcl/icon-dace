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

// Additional copyright for stacktrace
// Copyright (C) 2011 by Aivars Kalvans <aivars.kalvans@gmail.com>
// SPDX-License-Identifier: MIT

#include "config.h"

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>

#define FS_BACKTRACE_NONE 0
#define FS_BACKTRACE_UNWIND 1
#define FS_BACKTRACE_EXECINFO 2
#define FS_BACKTRACE_NEC 100

#ifndef FS_BACKTRACE
#if defined(__NEC__)
#define FS_BACKTRACE FS_BACKTRACE_NEC
#elif defined(HAVE_UNWIND_H) && defined(__ELF__)                  \
    && (defined(__linux) || defined(__linux__) || defined(linux)) \
    && _POSIX_C_SOURCE >= 2
#define FS_BACKTRACE FS_BACKTRACE_UNWIND
#elif defined(HAVE_EXECINFO_H)
#define FS_BACKTRACE FS_BACKTRACE_EXECINFO
#else
#define FS_BACKTRACE FS_BACKTRACE_NONE
#endif
#endif

#if FS_BACKTRACE == FS_BACKTRACE_UNWIND
#include <unwind.h>
#elif FS_BACKTRACE == FS_BACKTRACE_EXECINFO
#include <execinfo.h>
#endif

#include "util_backtrace.h"

#if FS_BACKTRACE == FS_BACKTRACE_NEC
void util_backtrace(void) {
    /* Requirements:
    Code must be compiled and linked with -traceback=verbose
    Code must be executed with environment variable VE_TRACEBACK=VERBOSE
    */
    unsigned long *frame_address = __builtin_frame_address(0);
    __builtin_traceback(frame_address);
    fflush(stdout);
}
#elif FS_BACKTRACE == FS_BACKTRACE_UNWIND
struct stacktrace_frame {
    void *addr;
    char *file;
    char *func;
    int line;
};

struct file_map {
    unsigned long long start;
    unsigned long long end;
    unsigned long long offset;
    int is_pic;
    char *file;
};

struct stacktrace {
    char *exe;  /* Name of executable file */
    char *maps; /* Process memory map for this snapshot */

    size_t frames_size;
    size_t frames_len;
    struct stacktrace_frame *frames;

    size_t files_len;
    struct file_map *files;

    unsigned skip;
};

static char *read_whole_file(char *fname) {
    /* procfs files don't have size, so read chunks until EOF is reached */
    char *data = NULL;

    FILE *f = fopen(fname, "r");
    if (f != NULL) {
        int n, len = 0;
        int size = 1024;

        data = (char *) malloc(size);
        for (;;) {
            int max = size - len;

            n = fread(data + len, 1, max, f);
            if (n > 0) {
                len += n;
            }

            if (n != max) {
                break;
            }
            size *= 2;
            data = (char *) realloc(data, size);
        }
        data[len] = '\0';
        fclose(f);
    }
    return data;
}

static _Unwind_Reason_Code collect(struct _Unwind_Context *ctx, void *p) {
    struct stacktrace *trace = (struct stacktrace *) p;

    if (trace->skip > 0) {
        trace->skip--;
        return _URC_NO_REASON;
    }

    struct stacktrace_frame frame;
    frame.addr = (void *) _Unwind_GetIP(ctx);
    frame.file = NULL;
    frame.func = NULL;
    frame.line = 0;

    if (trace->frames_len == trace->frames_size) {
        trace->frames_size = trace->frames_size * 2;
        trace->frames = (struct stacktrace_frame *) realloc(
            trace->frames,
            sizeof(struct stacktrace_frame) * trace->frames_size);
    }
    trace->frames[trace->frames_len++] = frame;
    return _URC_NO_REASON;
}

static struct stacktrace *stacktrace_get(unsigned skip) {
    struct stacktrace *trace
        = (struct stacktrace *) malloc(sizeof(struct stacktrace));
    trace->skip = skip + 1;
    trace->maps = NULL;
    trace->exe = NULL;

    trace->frames_size = 128;
    trace->frames_len = 0;
    trace->frames = (struct stacktrace_frame *) malloc(
        sizeof(struct stacktrace_frame) * trace->frames_size);

    trace->files_len = 0;
    trace->files = NULL;

    char procf[512];
    snprintf(procf, sizeof(procf), "/proc/%d/exe", (int) getpid());

    for (size_t len = 512;; len *= 2) {
        trace->exe = (char *) realloc(trace->exe, len);

        size_t n = readlink(procf, trace->exe, len);
        if (n == -1) {
            break;
        }
        if (n < len) {
            break;
        }
    }

    snprintf(procf, sizeof(procf), "/proc/%d/maps", (int) getpid());
    trace->maps = read_whole_file(procf);

    _Unwind_Backtrace(collect, trace);

    return trace;
}

static void stacktrace_free(struct stacktrace *trace) {
    if (trace != NULL) {
        free(trace->exe);
        free(trace->maps);
        if (trace->frames != NULL) {
            for (int i = 0; i < trace->frames_len; i++) {
                free(trace->frames[i].func);
                free(trace->frames[i].file);
            }
        }
        free(trace->frames);
        free(trace->files);
        free(trace);
    }
}

static void read_map(struct stacktrace *trace) {
    if (trace->files_len > 0) {
        return;
    }

    size_t files_size = 1;
    trace->files
        = (struct file_map *) malloc(sizeof(struct file_map) * files_size);

    char *saveptr;
    char *line = strtok_r(trace->maps, "\n", &saveptr);
    while (line != NULL) {
        char *p;
        char *saveptr2;
        unsigned long long start, end, offset;
        char *name;

        if (trace->files_len >= files_size) {
            files_size *= 2;
            trace->files = (struct file_map *) realloc(
                trace->files, sizeof(struct file_map) * files_size);
        }

        /* sscanf requires different format strings for 32/64 bits :( */
        p = strtok_r(line, "-", &saveptr2);
        if (p != NULL) {
            start = strtoull(p, NULL, 16);
            p = strtok_r(NULL, " ", &saveptr2);
            if (p != NULL) {
                end = strtoull(p, NULL, 16);
                p = strtok_r(NULL, " ", &saveptr2);
                if (p != NULL) {
                    p = strtok_r(NULL, " ", &saveptr2);
                    if (p != NULL) {
                        offset = strtoull(p, NULL, 16);
                        p = strtok_r(NULL, " ", &saveptr2);
                        if (p != NULL) {
                            p = strtok_r(NULL, " ", &saveptr2);
                            if (p != NULL) {
                                p = strtok_r(NULL, " ", &saveptr2);
                                if (p != NULL) {
                                    name = p;
                                } else {
                                    name = trace->exe;
                                }

                                trace->files[trace->files_len].start = start;
                                trace->files[trace->files_len].end = end;
                                trace->files[trace->files_len].offset = offset;
                                trace->files[trace->files_len].file = name;
                                trace->files_len++;
                            }
                        }
                    }
                }
            }
        }

        line = strtok_r(NULL, "\n", &saveptr);
    }

    if (trace->files_len > 0) {
        // To avoid reading the ELF tables, we assume that the executable IS
        // position-independent if the first memory region is NOT mapped at the
        // absolute virtual addresses 0x400000 (64-bit) and 0x08048000 (32-bit):
        trace->files[0].is_pic = !(trace->files[0].start == 0x400000
                                   || trace->files[0].start == 0x08048000);
        // The rest of the memory regions are position-independent if they are
        // loaded from other files, i.e. shared libraries (we do not care about
        // regions that are not associated with real files, e.g. [heap], [stack]
        // or anonymous):
        for (size_t i = 1; i < trace->files_len; i++) {
            trace->files[i].is_pic
                = trace->files[0].is_pic
                  || (strcmp(trace->files[0].file, trace->files[i].file) != 0);
        }
    }
}

static struct file_map *_find_file(struct stacktrace *trace,
                                   unsigned long long addr) {
    for (size_t i = 0; i < trace->files_len; i++) {
        if (trace->files[i].start <= addr && addr <= trace->files[i].end) {
            return &trace->files[i];
        }
    }
    return NULL;
}

static void _addr2line(struct stacktrace_frame *frame, struct file_map *file) {
    // We need to calculate the address to pass it to the addr2line utility. For
    // non-PIC (non-PIE) objects, which use the absolute addressing, it is
    // simply the address of the frame:
    unsigned long long addr = (unsigned long long) frame->addr;

    // For PIC (PIE) objects, which use the relative addressing, we subtract the
    // starting address of the memory region and add the file offset of the
    // mapping (the best we can do without parsing the ELF table):
    if (file->is_pic) {
        addr += file->offset - file->start;
    }

    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "addr2line -C -f -e %s 0x%llx", file->file,
             addr);
    /*printf("CMD: %s\n", cmd);*/

    FILE *f = popen(cmd, "r");
    if (f == NULL) {
        return;
    }

    char line[1024];
    if (fgets(line, sizeof(line), f) != NULL) {
        char *p = strchr(line, '\n');
        if (p != NULL) {
            *p = '\0';
        }
        frame->func = strdup(line);
        if (fgets(line, sizeof(line), f) != NULL) {
            p = strchr(line, ':');

            if (p != NULL) {
                *p++ = '\0';
                frame->line = atoi(p);
            }
            if (strcmp(line, "??") == 0) {
                frame->file = strdup(file->file);
            } else {
                frame->file = strdup(line);
            }
        }
    }
    pclose(f);
}

static void stacktrace_resolve(struct stacktrace *trace) {
    read_map(trace);

    for (size_t i = 0; i < trace->frames_len; i++) {
        struct file_map *file;
        file = _find_file(trace, (unsigned long long) trace->frames[i].addr);
        if (file == NULL) {
            continue;
        }
        _addr2line(&trace->frames[i], file);
    }
}

static void stacktrace_fprint(struct stacktrace *trace, FILE *f) {
    stacktrace_resolve(trace);

    for (size_t i = 0; i < trace->frames_len; i++) {
        struct stacktrace_frame *frame = &trace->frames[i];
        fprintf(f, "[util_backtrace]: #%ld %p - %s in %s:%d\n", i, frame->addr,
                frame->func ? frame->func : "??",
                frame->file ? frame->file : "??", frame->line);
    }
}

static void stacktrace_print(struct stacktrace *trace) {
    stacktrace_fprint(trace, stdout);
}

void util_backtrace(void) {
    struct stacktrace *trace = stacktrace_get(0);
    stacktrace_print(trace);
    stacktrace_free(trace);
    fflush(stdout);
}
#elif FS_BACKTRACE == FS_BACKTRACE_EXECINFO
void util_backtrace(void) {
    void *array[32];
    // TODO: we should detect the actual type of the size argument of the
    // backtrace and backtrace_symbols functions and use it here. The type is
    // int on most platforms but size_t on some other (e.g. FreeBSD). See how
    // the check is implemented in libcdi.
    int size = backtrace(array, 32);
    char **strings = backtrace_symbols(array, size);

    for (int i = 0; i < size; i++) {
        fprintf(stdout, "[util_backtrace]: %s\n", strings[i]);
    }

#if defined(__ELF__)
    fprintf(stdout, "[util_backtrace]: use addr2line for addresses to line "
                    "number conversion.\n");
#endif

    free(strings);
    fflush(stdout);
}
#else
void util_backtrace(void) {
    printf("[util_backtrace]: No backtrace available\n");
    fflush(stdout);
}
#endif
