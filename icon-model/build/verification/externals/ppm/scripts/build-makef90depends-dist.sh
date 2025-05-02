#! /bin/bash

version=$(sed -n -e '/^AC_INIT/{' -e \
              's/AC_INIT([^,]*,\[\([0-9.]*\)\],.*/\1/;p;q;}' configure.ac)
archive="makef90depends-${version}.tar.gz"
if [[ -f "${archive}" ]]; then
  echo "archive for new version exists: ${archive}" >&2
  exit 1
fi
tar -czf "${archive}" \
  scripts/makef90depends \
  scripts/cpp.pm \
  util/*wrapper \
  doc/makef90depends
