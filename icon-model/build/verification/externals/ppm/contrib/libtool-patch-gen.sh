#! /bin/bash
orig="$HOME/opt/libtool-2.4.2-x64-linux/share/libtool/config/ltmain.sh"
declare -a temp
temp=($(mktemp))
trap "rm \"${temp[@]}\"" EXIT
cp "$orig" "${temp[0]}"
for patchfile in contrib/*/*-libtool.patch ; do
  patch "${temp[0]}" $patchfile
done
diff -u \
  "${temp[0]}" \
  config/ltmain.sh \
  | sed "s|$orig|a/config/ltmain.sh|
s| config/ltmain.sh| b/config/ltmain.sh|"
