#
# Source this file for MakeExperiments! bash utilities
#

# Change to script directory ($SCRIPT_DIR) of a given .config.
# If no .config file is given, the 'update' file is sought and evaluated.
# '-d xxx' will change to $XXX_DIR instead of $SCRIPT_DIR

cdexp () {
    local var=script dir cfg buffer vardir OPTIND
    while getopts :d:h OPTOPT
    do
        case $OPTOPT in
            d)  var="$OPTARG";;
            h)  exec >&2
                echo "Usage: cdexp [-d DIRSPEC] [-h] [CONFIG]"
                echo
                echo "DIRSPEC is converted to upper case and '_DIR' is added"
                echo "If CONFIG is omitted, the 'update' file must be available"
                return 0
                ;;
            \?) echo "Oops: invalid option '$OPTARG'" >&2; return 1;;
        esac
    done
    shift $((OPTIND - 1))
    
    if [[ "$1" ]]
    then
        dir=$(dirname "$1")
        cfg="$1"
        shift
    elif [[ -f update ]]
    then
        eval $(python -c '
import update
u = update.Update("update")
print("dir="+u.get_config_dir()+"\ncfg="+u.get_config_file())
        ')
    else
        echo 'Oops: invalid number of parameters' >&2
        return 1
    fi
    var=${var^^}
    var=${var%_DIR}
    vardir=$(builtin cd "$dir" && getexp -k "${var}_DIR" "$cfg" "$@") &&
        cd "$vardir"
}
