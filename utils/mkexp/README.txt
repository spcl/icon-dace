=================================================================
Make Experiments! - Run-script generation for earth system models
=================================================================

Make Experiments! (mkexp) contains a set of tools for preparing experiments
with the MPI-M's earth system models.

mkexp
    takes a simple configuration file from the user and the templates provided
    by the model to generate a set of shell scripts that - when run - will set
    up directories, prepare input data, define configurations, and launch the
    actual model run and processing.

diffexp
    is a tool that - given two mkexp configuration files - will compare the
    generated scripts with either the standard diff program or any equivalent
    tool, taking care of differences by experiment names etc. It is complemented
    by diffpath which does the same kind of diff but takes its information from
    the command line instead.

rmexp, cpexp, getexp
    More tools to remove experiment data, create a copy under a new name,
    or get experiment info in shell-readable form.

