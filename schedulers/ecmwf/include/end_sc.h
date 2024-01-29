#==========================================================================
#  Send end-task and exit.
#==========================================================================

smscomplete
trap 0

# Cleanup of link to job output file.
 
if [[ $HOST = @(cc*) ]]; then
  [[ -L $_running_output ]] && rm -f $_running_output
fi

set +e

exit
