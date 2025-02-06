#!/bin/bash

# Build the documentation using FORD

FORD_MD='project.md'

no_search=""
quiet=""
if [[ $# > 0 ]]; then
  case $1 in
    -f)
      no_search="--no-search"
      ;;
    -q)
      quiet="-q"
      ;;
  esac
fi

if ! hash python3 2>/dev/null; then
  echo "python3 not found"
  exit 1
fi

if grep 'graph: *true' ${FORD_MD}; then
  if ! hash dot; then
    echo "graph: true set in ${FORD_MD} but no installation of graphviz found"
    exit 1
  fi
fi

mkdir 3rdparty >& /dev/null || true
cd 3rdparty
[[ -d ford ]] || git clone https://github.com/RMShur/ford.git
cd ford
git pull --rebase
pip3 install --user -r requirements.txt
pip3 install --user lxml
cd ../..

output_dir=$(grep output_dir project.md | cut -d' ' -f2)

rev=`git rev-parse --short HEAD`
#rev=`git log --pretty='format:%h %D' --first-parent | grep HEAD | cut -f1 -d' '`
[[ -z $(git status -suno) ]] && dirty="" || dirty="(dirty)"
if [[ -z $GITLAB_CI ]]; then
  branch=`git log --pretty='format:%h %D' --first-parent | grep HEAD | cut -f4 -d' ' | sed 's/,//'`
else
  branch=$CI_COMMIT_REF_NAME
fi

echo "Building documentation (Revision:${rev}${dirty} Branch:${branch}) ..."
mkdir pp_temp$$
python3 scripts/dsl4jsb/dsl4jsb.py -n -d src -t pp_temp$$ -k
find pp_temp$$ -type f -name '*.f90' -exec \
    sed -i -e '1s/!> /!> Summary: /' \
           -e '1,20s/icon-model.org/<https:\/\/icon-model.org>  /' \
           -e '1,20s/!> ICON-Land/!>#### ICON-Land/' \
           -e '1,20s/MPI-BGC$/MPI-BGC  /' \
           -e '1,20s/AUTHORS.md$/AUTHORS.md  /' \
           -e '1,20s/license information$/license information  /' \
    '{}' \;
python3 ./ford.py $no_search $quiet -r "Revision:${rev}${dirty} Branch:${branch}" -d pp_temp$$ ${FORD_MD}
find $output_dir -type f -name '*.html' -exec sed -i 's/was developed by/is developed by/' '{}' \;
rm -rf pp_temp$$

