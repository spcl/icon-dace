#!/usr/bin/env python2

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

import optparse
import sys
import os
import shutil

#==============================================================================
# CONFIG + HELPERS {{{
DEFAULT_DATA_DIR = '/pool/data/ICON'
DEFAULTS_TARGETS = {
        'FORCING' : 'ocean-flux.nc',
        'INIT'    : 'initial_state.nc',
        'REALX'   : 'ocean-relax.nc',
        'DUST'    : 'dust.nc',
        'NITRO'   : 'nitrogen.nc',
        }

parser = optparse.OptionParser()
parser.add_option('-v', '--verbose'         , dest="verbose"                , default=True  , action="store_true"       , )
parser.add_option("-q", '--quiet'           , dest="verbose"                , action="store_false",)
parser.add_option('-n', '--dryrun'          , dest="dryrun"                 , default=False , action="store_true"       , )
parser.add_option('-g', '--grid'            , dest="grid_filename"          , action="store", help="gridfile to be used", )
parser.add_option('-i', '--init'            , dest="initialization_filename", action="store", )
parser.add_option('-f', '--froce'           , dest="forcing_filename"       , action="store", )
parser.add_option('-r', '--relax'           , dest="relaxation_filename"    , action="store", )
parser.add_option('-d', '--data-dir'        , dest="data_dir"               , action="store", default=DEFAULT_DATA_DIR  , )
parser.add_option('-e', '--experiment-dir'  , dest="experiment_dir"         , action="store", )
parser.add_option('-a', '--action'          , dest="action"                 , action="store", default="link"            , help="can be link or copy (default:link)")
# -----------------------------------------------------------------------------
# litte wrapper for link/copy
def action(type,source_dir, dest_dir, source_file,dest_file,dryrun=False,verbose=False):
    if 'link' == type:
        if (dryrun):
            print('ln -s -f %s %s'%('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file])))
        else:
            if os.path.exists(dest_file):
                os.remove(dest_file)
            if verbose:
                print('ln -s -f %s %s'%('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file])))
            os.symlink('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file]))
    else:
        if (dryrun):
            print('cp -f %s %s'%('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file])))
        else:
            if verbose:
                print('cp -f %s %s'%('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file])))
            shutil.copyfile('/'.join([source_dir,source_file]),'/'.join([dest_dir,dest_file]))

#==============================================================================
# MAIN {{{
#
# parse command line options
options, remainder = parser.parse_args()
if options.verbose:
    print('ARGV      :', sys.argv[1:])
    print('OPTIONS   :', options)
    print('REMAINING :', remainder)

# create file targes
sourceDir = options.data_dir
destDir   = options.experiment_dir
action(options.action, sourceDir, destDir, options.grid_filename          , options.grid_filename      , options.dryrun, options.verbose)
action(options.action, sourceDir, destDir, options.initialization_filename, DEFAULTS_TARGETS['INIT']   , options.dryrun, options.verbose)
action(options.action, sourceDir, destDir, options.forcing_filename       , DEFAULTS_TARGETS['FORCING'], options.dryrun, options.verbose)
action(options.action, sourceDir, destDir, options.relaxation_filename    , DEFAULTS_TARGETS['REALX']  , options.dryrun, options.verbose)
action(options.action, sourceDir, destDir, options.dust_filename          , DEFAULTS_TARGETS['DUST']   , options.dryrun, options.verbose)
action(options.action, sourceDir, destDir, options.nitro_filename         , DEFAULTS_TARGETS['NITRO']   , options.dryrun, options.verbose)
# # }}} =======================================================================
