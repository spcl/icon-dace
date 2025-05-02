#!/bin/bash

set -eu
set -o pipefail

script_dir=$(cd "$(dirname "$0")"; pwd)
top_srcdir=$(cd "${script_dir}/../.."; pwd)

test_environment='unknown'
case $(uname -s) in
  Linux)
    node_name=$(uname -n)
    case $(host "${node_name}") in
      mpipc*.mpimet.mpg.de\ * | \
      breeze*.mpimet.mpg.de\ *)
        test_environment='breeze-mpim' ;;
      levante*.lvt.dkrz.de\ *)
        test_environment='levante-dkrz' ;;
      daint*.login.cscs.ch\ *)
        test_environment='daint-cscs' ;;
    esac
esac

test "x${test_environment}" != 'xunknown' || {
  echo 'ERROR: unknown test environment' >&2
  exit 1
}

cd "${top_srcdir}"

git update-index -q --refresh || {
  echo "ERROR: failed to update git index in '${top_srcdir}'" >&2
  exit 1
}

git diff-files --quiet || {
  echo "ERROR: '${top_srcdir}' has unstaged changes" >&2
  exit 1
}

untracked_files=`git ls-files --others` || {
  echo "ERROR: failed to get list of untracked files in '${top_srcdir}'" >&2
  exit 1
}
if test -n "${untracked_files}"; then
  {
    echo "ERROR: '${top_srcdir}' has untracked files:" >&2
    echo "${untracked_files}" >&2
  } >&2
  exit 1
fi

echo "Running tests for '${test_environment}'"

for script in $(find "${script_dir}/${test_environment}" -type f -name 'test.*' -executable); do
  echo "Running '${script}'"
  "${script}"
  git clean -fdx
  git checkout .
done
