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

#define XOPEN_SOURCE 500

#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

int util_islink(char *path) {
    char *buf;
    ssize_t len;
    long max_buf_len;
    int iret;

    iret = -1;

    max_buf_len = pathconf("/", _PC_NAME_MAX);

    buf = (char *) malloc(max_buf_len * sizeof(char));

    if ((len = readlink(path, buf, max_buf_len - 1)) != -1) {
        buf[len] = '\0';
        /* file is a link ..., return C true */
        iret = 1;
    } else {
        if (errno == EINVAL) {
            /* file is not a link ... */
            iret = 0;
        }
    }

    free(buf);

    /* something else went wrong ... */
    return iret;
}

static int direxists(const char *path) {
    struct stat statbuf;
    return stat(path, &statbuf) == 0 && S_ISDIR(statbuf.st_mode);
}

static const char *find_tmpdir(void) {
    if (direxists(P_tmpdir))
        return P_tmpdir;
    else if (strcmp(P_tmpdir, "/tmp") != 0 && direxists("/tmp"))
        return "/tmp";
    else
        return NULL;
}

int util_create_tmpfile(char *filename, const int max_len) {
    const char *dir = find_tmpdir();
    if (dir == NULL) return -1;

    const char *basename_prefix = "icon";
    const int pid_infix_len = 8;
    const char *templ_suffix = "XXXXXX";

    int filename_len = strlen(dir) + 1 + strlen(basename_prefix) + pid_infix_len
                       + strlen(templ_suffix) + 1;

    if (max_len < filename_len) return -1;

    const pid_t pid = getpid();
    char pid_infix[16];
    sprintf(pid_infix, "%0*ld", pid_infix_len, (long) pid);
    pid_infix[pid_infix_len] = 0;

    sprintf(filename, "%s/%s%s%s", dir, basename_prefix, pid_infix,
            templ_suffix);

    int fd = mkstemp(filename);
    if (fd < 0) return -1;

    close(fd);

    return filename_len;
}

long int util_filesize(char *filename) {
    struct stat statbuf;

    if (stat(filename, &statbuf) == -1) {
        return 0;
    }

    return ((long int) statbuf.st_size);
}

int util_file_is_writable(char *filename) {
    int result = 0;
    int rval = access(filename, W_OK);
    if (rval == 0)
        return 1;
    else if ((errno == EACCES) || (errno == EROFS))
        return 0;
    return 0;
}

// This wrapper to symlink() does not fail if there is already a symlink at
// linkName; it will simply overwrite the existing symlink. However, if
// something exists at linkName which is not a symlink, the error EEXIST is
// returned.
//
// May return EEXIST or any error returned by lstat(), unlink(), or symlink(),
// returns zero on success.
int createSymlink(const char *targetPath, const char *linkName) {
    errno = 0;
    struct stat fileInfo;
    if (!lstat(linkName, &fileInfo)) {
        if ((fileInfo.st_mode & S_IFMT) == S_IFLNK) {
            if (unlink(linkName)) return errno;
        } else {
            return EEXIST;  // Something exists at linkName, and it's not a
                            // symlink. Bail out.
        }
    } else if (errno != ENOENT)  // ENOENT is the only non-fatal error, it just
                                 // means that there is no such link there yet
    {
        return errno;
    }

    // At this point we know that there is no file at linkName.
    int err = symlink(targetPath, linkName);
    return err != 0 ? errno : err;
}
