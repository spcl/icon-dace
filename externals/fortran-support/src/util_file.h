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

#ifndef UTIL_FILE_H
#define UTIL_FILE_H

#ifdef __cplusplus
extern "C" {
#endif

int util_islink(char *path);
long int util_filesize(char *filename);
int util_file_is_writable(char *filename);

#ifdef __cplusplus
}
#endif

#endif  // UTIL_FILE_H
