#!/usr/bin/env python3

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import os
import subprocess
import sys

HEADER = """\
#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <string.h>

#if defined(__cplusplus)
extern "C"
{
#  define LANG_PREFIX cxx_
#else
#  define LANG_PREFIX c_
#endif

#define GLUE_helper(x, y) x##y
#define GLUE(x, y) GLUE_helper(x, y)

#define CHARS1(x) ('0' + ((x) % 10))
#define CHARS2(x) ('0' + (((x) / 10) % 10)), CHARS1(x)
#define CHARS3(x) ('0' + (((x) / 100) % 10)), CHARS2(x)
#define CHARS4(x) ('0' + (((x) / 1000) % 10)), CHARS3(x)
#define CHARS5(x) ('0' + (((x) / 10000) % 10)), CHARS4(x)

#define DEFINE_STRING_GETTER_helper(x)                             \\
  void pvcs_get_##x(const char **const val, size_t *const val_len) \\
  {                                                                \\
    *val = pvcs_##x;                                               \\
    *val_len = (*val == NULL) ? 0 : strlen(*val);                  \\
  }
#define DEFINE_STRING_GETTER(x) DEFINE_STRING_GETTER_helper(x)

#if !defined(COMPILER_VERSION_ONLY)

struct pvcs_info
{
  const char *const name;
  const char *const revision;
  const char *const remote_url;
  const char *const local_branch;
};

static const struct pvcs_info pvcs_infos[] = {"""

ELEMENT_TEMPLATE = """\
  // {comment}
  {{ .name = {name},
    .revision = {revision},
    .remote_url = {remote_url},
    .local_branch = {local_branch}, }},"""

FOOTER = """\
};

static const struct pvcs_info *
pvcs_find_info(const char *key, const size_t key_len)
{
  static const size_t count = sizeof(pvcs_infos) / sizeof(struct pvcs_info);
  for (size_t ii = 0; ii < count; ++ii)
    {
      const struct pvcs_info *info = &pvcs_infos[ii];
      if (strlen(info->name) == key_len
          && strncmp(info->name, key, key_len) == 0)
        {
          return info;
        }
    }
  return NULL;
}

#  define DEFINE_INFO_GETTER(x)                                          \\
    int pvcs_get_##x(const char *const key, const size_t key_len,        \\
                     const char **const val, size_t *const val_len)      \\
    {                                                                    \\
      const struct pvcs_info *const info = pvcs_find_info(key, key_len); \\
      *val = (info == NULL) ? NULL : info->x;                            \\
      *val_len = (*val == NULL) ? 0 : strlen(info->x);                   \\
      return (info == NULL) ? 1 : 0;                                     \\
    }

DEFINE_INFO_GETTER(revision)
DEFINE_INFO_GETTER(remote_url)
DEFINE_INFO_GETTER(local_branch)

static const char *const pvcs_icon_version =
#  if defined(PACKAGE_VERSION)
  PACKAGE_VERSION
#  else
  NULL
#  endif
  ;

DEFINE_STRING_GETTER(icon_version)

#endif

// Intel
#if defined(__INTEL_LLVM_COMPILER)
#  define COMPILER_NAME "Intel"
#  if __INTEL_LLVM_COMPILER < 1000000L
#    define COMPILER_VERSION_MAJOR __INTEL_LLVM_COMPILER / 100
#    define COMPILER_VERSION_MINOR __INTEL_LLVM_COMPILER / 10 % 10
#    undef COMPILER_VERSION_PATCH
#  else
#    define COMPILER_VERSION_MAJOR __INTEL_LLVM_COMPILER / 10000
#    define COMPILER_VERSION_MINOR __INTEL_LLVM_COMPILER / 100 % 100
#    define COMPILER_VERSION_PATCH __INTEL_LLVM_COMPILER % 100
#  endif
// Intel Classic
#elif defined(__INTEL_COMPILER)
#  define COMPILER_NAME "Intel Classic"
#  if __INTEL_COMPILER < 2021
#    define COMPILER_VERSION_MAJOR __INTEL_COMPILER / 100
#    define COMPILER_VERSION_MINOR __INTEL_COMPILER / 10 % 10
#    if __INTEL_COMPILER_BUILD_DATE == 20181018 \\
      || __INTEL_COMPILER_BUILD_DATE == 20200306
#      define COMPILER_VERSION_PATCH 1
#    elif defined(__INTEL_COMPILER_UPDATE)
#      define COMPILER_VERSION_PATCH __INTEL_COMPILER_UPDATE
#    else
#      undef COMPILER_VERSION_PATCH
#    endif
#  else
#    define COMPILER_VERSION_MAJOR __INTEL_COMPILER
#    define COMPILER_VERSION_MINOR __INTEL_COMPILER_UPDATE
#    if __INTEL_COMPILER_BUILD_DATE == 20201208
#      define COMPILER_VERSION_PATCH 2
#    else
#      define COMPILER_VERSION_PATCH 0
#    endif
#  endif
// Cray
#elif defined(__cray__)
#  define COMPILER_NAME "Cray"
#  define COMPILER_VERSION_MAJOR __cray_major__
#  define COMPILER_VERSION_MINOR __cray_minor__
#  define COMPILER_VERSION_PATCH __cray_patchlevel__
// Cray Classic
#elif defined(_CRAYC)
#  define COMPILER_NAME "Cray Classic"
#  define COMPILER_VERSION_MAJOR _RELEASE_MAJOR
#  define COMPILER_VERSION_MINOR _RELEASE_MINOR
#  define COMPILER_VERSION_PATCH _RELEASE_PATCHLEVEL
// NEC
#elif defined(__NEC__)
#  define COMPILER_NAME "NEC"
#  define COMPILER_VERSION_MAJOR __NEC_VERSION__ / 10000
#  define COMPILER_VERSION_MINOR __NEC_VERSION__ / 100 % 100
#  define COMPILER_VERSION_PATCH __NEC_VERSION__ % 100
// NVHPC
#elif defined(__NVCOMPILER)
#  define COMPILER_NAME "NVHPC"
#  define COMPILER_VERSION_MAJOR __NVCOMPILER_MAJOR__
#  define COMPILER_VERSION_MINOR __NVCOMPILER_MINOR__
#  define COMPILER_VERSION_PATCH __NVCOMPILER_PATCHLEVEL__
// PGI
#elif defined(__PGI)
#  define COMPILER_NAME "PGI"
#  define COMPILER_VERSION_MAJOR __PGIC__
#  define COMPILER_VERSION_MINOR __PGIC_MINOR__
#  define COMPILER_VERSION_PATCH __PGIC_PATCHLEVEL__
// Apple Clang
#elif defined(__apple_build_version__)
#  define COMPILER_NAME "Apple Clang"
#  define COMPILER_VERSION_MAJOR __clang_major__
#  define COMPILER_VERSION_MINOR __clang_minor__
#  define COMPILER_VERSION_PATCH __clang_patchlevel__
// Clang
#elif defined(__clang__)
#  define COMPILER_NAME "Clang"
#  define COMPILER_VERSION_MAJOR __clang_major__
#  define COMPILER_VERSION_MINOR __clang_minor__
#  define COMPILER_VERSION_PATCH __clang_patchlevel__
#else
#  undef COMPILER_NAME
#  undef COMPILER_VERSION_MAJOR
#  undef COMPILER_VERSION_MINOR
#  undef COMPILER_VERSION_PATCH
#endif

// GCC: either the primary or the secondary compiler
#if defined(__GNUC__)
#  if defined(COMPILER_NAME)
#    define COMPILER_SECONDARY_NAME "GNU"
#    define COMPILER_SECONDARY_VERSION_MAJOR __GNUC__
#    define COMPILER_SECONDARY_VERSION_MINOR __GNUC_MINOR__
#    define COMPILER_SECONDARY_VERSION_PATCH __GNUC_PATCHLEVEL__
#  else
#    define COMPILER_NAME "GNU"
#    define COMPILER_VERSION_MAJOR __GNUC__
#    define COMPILER_VERSION_MINOR __GNUC_MINOR__
#    define COMPILER_VERSION_PATCH __GNUC_PATCHLEVEL__
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_name)) =
#if defined(COMPILER_NAME)
  COMPILER_NAME
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_name))

#undef PATCH_CHARS
#undef MINOR_CHARS
#undef MAJOR_CHARS
#if defined(COMPILER_VERSION_MAJOR)
#  if COMPILER_VERSION_MAJOR < 10
#    define MAJOR_CHARS CHARS1(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 100
#    define MAJOR_CHARS CHARS2(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 1000
#    define MAJOR_CHARS CHARS3(COMPILER_VERSION_MAJOR)
#  elif COMPILER_VERSION_MAJOR < 10000
#    define MAJOR_CHARS CHARS4(COMPILER_VERSION_MAJOR)
#  else
#    define MAJOR_CHARS CHARS5(COMPILER_VERSION_MAJOR)
#  endif
#  if defined(COMPILER_VERSION_MINOR)
#    if COMPILER_VERSION_MINOR < 10
#      define MINOR_CHARS CHARS1(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 100
#      define MINOR_CHARS CHARS2(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 1000
#      define MINOR_CHARS CHARS3(COMPILER_VERSION_MINOR)
#    elif COMPILER_VERSION_MINOR < 10000
#      define MINOR_CHARS CHARS4(COMPILER_VERSION_MINOR)
#    else
#      define MINOR_CHARS CHARS5(COMPILER_VERSION_MINOR)
#    endif
#    if defined(COMPILER_VERSION_PATCH)
#      if COMPILER_VERSION_PATCH < 10
#        define PATCH_CHARS CHARS1(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 100
#        define PATCH_CHARS CHARS2(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 1000
#        define PATCH_CHARS CHARS3(COMPILER_VERSION_PATCH)
#      elif COMPILER_VERSION_PATCH < 10000
#        define PATCH_CHARS CHARS4(COMPILER_VERSION_PATCH)
#      else
#        define PATCH_CHARS CHARS5(COMPILER_VERSION_PATCH)
#      endif
#    endif
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_version)) =
#if defined(MAJOR_CHARS)
  (const char[])
  {
    MAJOR_CHARS,
#  if defined(MINOR_CHARS)
    '.', MINOR_CHARS,
#    if defined(PATCH_CHARS)
    '.', PATCH_CHARS,
#    endif
#  endif
    '\\0',
  }
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_version))

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_secondary_name)) =
#if defined(COMPILER_SECONDARY_NAME)
  COMPILER_SECONDARY_NAME
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_secondary_name))

#undef PATCH_CHARS
#undef MINOR_CHARS
#undef MAJOR_CHARS
#if defined(COMPILER_SECONDARY_VERSION_MAJOR)
#  if COMPILER_SECONDARY_VERSION_MAJOR < 10
#    define MAJOR_CHARS CHARS1(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 100
#    define MAJOR_CHARS CHARS2(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 1000
#    define MAJOR_CHARS CHARS3(COMPILER_SECONDARY_VERSION_MAJOR)
#  elif COMPILER_SECONDARY_VERSION_MAJOR < 10000
#    define MAJOR_CHARS CHARS4(COMPILER_SECONDARY_VERSION_MAJOR)
#  else
#    define MAJOR_CHARS CHARS5(COMPILER_SECONDARY_VERSION_MAJOR)
#  endif
#  if defined(COMPILER_SECONDARY_VERSION_MINOR)
#    if COMPILER_SECONDARY_VERSION_MINOR < 10
#      define MINOR_CHARS CHARS1(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 100
#      define MINOR_CHARS CHARS2(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 1000
#      define MINOR_CHARS CHARS3(COMPILER_SECONDARY_VERSION_MINOR)
#    elif COMPILER_SECONDARY_VERSION_MINOR < 10000
#      define MINOR_CHARS CHARS4(COMPILER_SECONDARY_VERSION_MINOR)
#    else
#      define MINOR_CHARS CHARS5(COMPILER_SECONDARY_VERSION_MINOR)
#    endif
#    if defined(COMPILER_SECONDARY_VERSION_PATCH)
#      if COMPILER_SECONDARY_VERSION_PATCH < 10
#        define PATCH_CHARS CHARS1(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 100
#        define PATCH_CHARS CHARS2(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 1000
#        define PATCH_CHARS CHARS3(COMPILER_SECONDARY_VERSION_PATCH)
#      elif COMPILER_SECONDARY_VERSION_PATCH < 10000
#        define PATCH_CHARS CHARS4(COMPILER_SECONDARY_VERSION_PATCH)
#      else
#        define PATCH_CHARS CHARS5(COMPILER_SECONDARY_VERSION_PATCH)
#      endif
#    endif
#  endif
#endif

static const char *const GLUE(pvcs_, GLUE(LANG_PREFIX, compiler_secondary_version)) =
#if defined(MAJOR_CHARS)
  (const char[])
  {
    MAJOR_CHARS,
#  if defined(MINOR_CHARS)
    '.', MINOR_CHARS,
#    if defined(PATCH_CHARS)
    '.', PATCH_CHARS,
#    endif
#  endif
    '\\0',
  }
#else
  NULL
#endif
  ;

DEFINE_STRING_GETTER(GLUE(LANG_PREFIX, compiler_secondary_version))

#if defined(__cplusplus)
}
#endif
"""

EXTRA_FOOTERS = {
    "art": """\
#if !defined(COMPILER_VERSION_ONLY)
#  define ART_WORKAROUND
#  if defined(ART_WORKAROUND)
#    define ART_DEFAULT_VALUE "unknown"
void
art_repository_url(char *name, int *actual_len)
{
  char const *val;
  size_t len;
  pvcs_get_remote_url("art", 3, &val, &len);
  if (val == NULL)
    {
      val = ART_DEFAULT_VALUE;
      len = strlen(ART_DEFAULT_VALUE);
    }
  if (len > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, val);
      *actual_len = len;
    }
}

void
art_branch_name(char *name, int *actual_len)
{
  char const *val;
  size_t len;
  pvcs_get_local_branch("art", 3, &val, &len);
  if (val == NULL)
    {
      val = ART_DEFAULT_VALUE;
      len = strlen(ART_DEFAULT_VALUE);
    }
  if (len > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, val);
      *actual_len = len;
    }
}

void
art_revision_key(char *name, int *actual_len)
{
  char const *val;
  size_t len;
  pvcs_get_revision("art", 3, &val, &len);
  if (val == NULL)
    {
      val = ART_DEFAULT_VALUE;
      len = strlen(val);
    }
  if (len > *actual_len)
    {
      *actual_len = 0;
    }
  else
    {
      strcpy(name, val);
      *actual_len = len;
    }
}
#  endif
#endif
"""
}


def run_cmd(*args, **kwargs):
    input_string = kwargs.get("input_string", None)

    if sys.version_info < (3, 3, 0):
        devnull = open(os.devnull)
    else:
        devnull = subprocess.DEVNULL

    proc = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=devnull,
        stdin=(subprocess.PIPE if input_string else devnull),
    )
    out, _ = proc.communicate(
        input=input_string.encode("utf-8") if input_string else None
    )

    if sys.version_info < (3, 3, 0):
        devnull.close()

    return str(out.decode("utf-8")), proc.returncode == 0


def run_git_cmd(repo_root, *args, **kwargs):
    return run_cmd(
        os.environ.get("PVCS_GIT", "git"),
        "--git-dir",
        os.path.join(repo_root, ".git"),
        "--work-tree",
        repo_root,
        *args,
        **kwargs
    )


def get_remote_urls(repo_root):
    remotes_raw, success = run_git_cmd(
        repo_root, "config", "--get-regexp", r"^remote\.[^.]+\.url"
    )

    remote_urls = {}
    if success and remotes_raw:
        for k, v in (
            line.split()
            for line in remotes_raw.split("\n")
            if line and not line.isspace()
        ):
            # remote.<name>.url
            remote_urls[k[7:-4]] = v

    return remote_urls


def is_in_remote(repo_root, remote_name, ref="HEAD"):
    remote_ref, success = run_git_cmd(
        repo_root,
        "for-each-ref",
        "--format=%(refname)",
        "--count=1",
        "--contains={0}".format(ref),
        "refs/remotes/{0}".format(remote_name),
    )
    return success and bool(remote_ref)


def get_revision(repo_root, tag_patterns):
    git_args = ["describe", "--abbrev=40", "--always", "--long", "--tags"]

    if tag_patterns:
        for p in tag_patterns.split(":"):
            if not p or p.isspace():
                git_args.append("--no-match")
            elif p.startswith("~"):
                p = p[1:]
                if not p or p.isspace():
                    git_args.append("--no-exclude")
                else:
                    git_args.extend(["--exclude", p])
            else:
                git_args.extend(["--match", p])

    revision, success = run_git_cmd(repo_root, *git_args)

    if success and revision:
        revision = revision.rstrip()
    else:
        return None

    # We are not using the --dirty argument of the 'git describe' command
    # because it uses 'git diff-index', which is sensitive to the bare touch of
    # a file in the repository. We also do not want changes in the submodules
    # to affect the "dirty" status of the main repository:
    diff_report, success = run_git_cmd(
        repo_root,
        "diff",
        "--abbrev=40",
        "--full-index",
        "--ignore-submodules=dirty",
        "--raw",
        "HEAD",
    )

    if not success:
        revision += "-broken"
    elif len(diff_report):
        revision += "-dirty"

    return revision


def get_local_branch(repo_root):
    branch, success = run_git_cmd(repo_root, "symbolic-ref", "HEAD")
    return os.path.basename(branch.rstrip()) if success and branch else None


def is_repo_root(repo_root):
    return os.path.exists(os.path.join(repo_root, ".git"))


def get_version_summary(repo_root, tag_patterns=None):
    result = {}

    revision = get_revision(repo_root, tag_patterns=tag_patterns)
    if revision:
        result["revision"] = revision

    remote_urls = get_remote_urls(repo_root)

    # Possible remotes are those that contain a reference to the current HEAD:
    possible_remotes = list(
        filter(
            lambda name_url: is_in_remote(repo_root, name_url[0], "HEAD"),
            remote_urls.items(),
        )
    )

    if len(possible_remotes) == 1:
        result["remote_url"] = possible_remotes[0][1]
    elif len(possible_remotes) > 1:
        possible_remotes.sort(
            key=lambda name_url: (
                # Prioritise remotes with sensible URLs:
                any(name_url[1].startswith(p) for p in ["git@", "https://"]),
                # Prioritise the open repository:
                any(
                    name_url[1] == os_repo
                    for os_repo in [
                        "git@gitlab.dkrz.de:icon/icon-model.git",
                        "https://gitlab.dkrz.de/icon/icon-model.git",
                    ]
                ),
                # Prioritise DKRZ-hosted repositories:
                "gitlab.dkrz.de" in name_url[1],
                # Prioritize the origin:
                name_url[0] == "origin",
                # Prioritize the upstream:
                name_url[0] == "upstream",
                # Sort the rest by the URL:
                name_url[1],
            )
        )
        result["remote_url"] = possible_remotes[-1][1]

    local_branch = get_local_branch(repo_root)
    if local_branch:
        result["local_branch"] = local_branch

    return result


def apply_element_template(name, path, version_summary):
    def prepare_value(value=None):
        if value is None:
            return "NULL"
        else:
            return '"{0}"'.format(value.replace('"', '\\"'))

    return ELEMENT_TEMPLATE.format(
        name=prepare_value(name),
        comment="{0} ({1})".format(name, path),
        revision=prepare_value(version_summary.get("revision", None)),
        remote_url=prepare_value(version_summary.get("remote_url", None)),
        local_branch=prepare_value(version_summary.get("local_branch", None)),
    )


def generate_version_c(srcdir=None, subdirs=None, stream=None):
    srcdir = os.path.abspath(srcdir) if srcdir else os.getcwd()

    icon_version_summary = {}
    if is_repo_root(srcdir):
        icon_version_summary = get_version_summary(
            srcdir, tag_patterns=os.environ.get("PVCS_ICON_TAG_PATTERNS", None)
        )

    version_c_lines = [
        HEADER,
        apply_element_template(
            "icon",
            srcdir,
            icon_version_summary,
        ),
    ]

    subdir_repos = []
    if subdirs:
        subdir_repos.extend(
            (
                # It might lead to name collisions in the future but currently we simply
                # use the basename as the subdirectory name:
                os.path.basename(subdir).lower(),
                os.path.join(srcdir, subdir),
            )
            for subdir in subdirs
        )

    for name, path in subdir_repos:
        # The path might not be a git submodule. We assume that there are no .git
        # directories in between path and srcdir:
        version_summary = icon_version_summary
        if is_repo_root(path):
            version_summary = get_version_summary(path)
        version_c_lines.append(apply_element_template(name, path, version_summary))

    version_c_lines.append(FOOTER)

    for name, _ in subdir_repos:
        if name in EXTRA_FOOTERS:
            version_c_lines.append(EXTRA_FOOTERS[name])

    if stream:
        stream.write("\n".join(version_c_lines))


def _parse_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="ICON version provenance collection tool. Generates a C source "
        "file containing information on the compiled version of ICON as well a its "
        "bundled libraries (repository URLs and revisions).",
        epilog="The tool honours the following environment variables: PVCS_GIT - git "
        "command or the full path to the git executable (default: git); "
        "PVCS_ICON_TAG_PATTERNS - colon-separated list of glob patterns, only tags "
        "matching one of the patterns are considered for the description of the ICON "
        "revision, if a pattern starts with the tilde (~), the tilde is dropped and "
        "the rest of the pattern is considered an exclude one: tags that match it are "
        "excluded from the consideration, an empty or a space-only element of the list "
        "clears and resets the list of match patterns, patterns that start with a "
        "tilde followed by zero or more spaces clear and reset the list of exclude "
        "patterns.",
    )

    parser.add_argument(
        "--srcdir",
        metavar="SRCDIR",
        help="root source directory of ICON (default: %(default)s)",
        default=".",
    )
    parser.add_argument(
        "--subdirs",
        metavar="SUBDIR",
        help="relative (to SRCDIR) paths to git submodules to collect the version "
        "information for",
        nargs="*",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="path to the generated source file; the modification timestamp of the "
        "file is updated only if its contents are updated; if not specified, the "
        "output is emitted to the standard output stream",
    )

    args = parser.parse_args()

    return args


def _main():
    args = _parse_args()

    if not args.output:
        generate_version_c(srcdir=args.srcdir, subdirs=args.subdirs, stream=sys.stdout)
    elif not os.path.exists(args.output):
        with open(args.output, "w") as version_c_stream:
            generate_version_c(
                srcdir=args.srcdir, subdirs=args.subdirs, stream=version_c_stream
            )
    else:
        import io

        version_c_string_stream = io.StringIO()
        generate_version_c(
            srcdir=args.srcdir, subdirs=args.subdirs, stream=version_c_string_stream
        )
        version_c_content = version_c_string_stream.getvalue()
        with open(args.output, "r+") as version_c_file_stream:
            if version_c_content != version_c_file_stream.read():
                version_c_file_stream.seek(0)
                version_c_file_stream.write(version_c_content)
                version_c_file_stream.truncate()


if __name__ == "__main__":
    _main()
