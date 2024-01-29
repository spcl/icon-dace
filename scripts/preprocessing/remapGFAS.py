#!/usr/bin/env python3

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
# Remap GFAS data to the ICON grid using CDO conservative remapping
# 
# Running the script with the option '-h' provides a list of options
# available to control input files, target grid and the calculation of weights.
#
# 01/2023 : D. Rieger (DWD)

import os
import sys
import getopt

inputfiles = []
targetgrids = []
recalcWeights=False

gribNmbDict = {'ocfire' : 'var90',
               'bcfire' : 'var91',
               'so2fire': 'var102'}

def create_gribrules_file(filename):
    f = open(filterfile,'w')
    f.write("""write "gfas_[cfVarName]_[dateTime].grib[editionNumber]";\n""")
    f.close()

def getFilesFromList(this_list):
    out_list = []
    for this_file in this_list:
        if this_file.split("_")[0] == "gfas" and this_file.split(".")[-1] == "grib1":
            out_list.append(this_file)
    return(out_list)

opts, args = getopt.getopt(sys.argv[1:],"hi:g:w",["ifile=","--targetgrid=","--recalc-weights"])
for opt, arg in opts:
    if opt == '-h':
        print ('Example: remapEdgar.py -i <inputfile1> -i <inputfile2> ... -g <targetgrid>')
        print ('Full list of options:')
        print ('-h : help')
        print ('-i : input file (multiple files can be added by separate -i entries')
        print ('-g : target grid (multiple files can be added by separate -g entries)')
        print ('-w : recalculate weightings')
        sys.exit()
    elif opt in ("-i", "--ifile"):
        inputfiles.append(arg)
    elif opt in ("-g", "--targetgrid"):
        targetgrids.append(arg)
    elif opt in ("-w", "--recalc-weights"):
        recalcWeights=True

nifiles=0
for inputfile in inputfiles:
    nifiles=nifiles+1
    print ('Input file '+str(nifiles)+' is', inputfile)
ngfiles=0
for targetgrid in targetgrids:
    ngfiles=ngfiles+1
    print ('Target grid '+str(ngfiles)+' is', targetgrid)



filterfile = "gribrules_file"

create_gribrules_file(filterfile)

ncfileList = []



for inputfile in inputfiles:
    # Split input file based on date and variable name
    print('Applying grib_filter to '+inputfile)
    os.system('grib_filter '+filterfile+' '+inputfile)
    grb1files = getFilesFromList( os.listdir(".") )
    for grb1file in grb1files:
        print("Processing "+grb1file)
        # Convert to netcdf
        ncfile_base  = grb1file.split(".")[0]
        tmpncfile = ncfile_base+'tmp.nc'
        species = grb1file.split("_")[1]
        os.system('cdo -f nc copy '+grb1file+' '+tmpncfile)
        os.system('cdo chname,'+gribNmbDict[species]+','+species+' '+tmpncfile+' '+ncfile_base+'.nc')
        os.system('rm '+grb1file)
        os.system('rm '+tmpncfile)
        ncfileList.append(ncfile_base+'.nc')

for targetgrid in targetgrids:
    # Calculate weightings for each target grid only once
    REMAP_WEIGHTS='remapWeights_'+targetgrid.split('/')[-1]
    if (recalcWeights or not os.path.exists(REMAP_WEIGHTS)):
        if (recalcWeights):
            print('Recalculating weightings as specified...')
        os.system('cdo -P 4 gencon,'+targetgrid+' '+inputfiles[0]+' '+REMAP_WEIGHTS)
    elif ( os.path.exists(REMAP_WEIGHTS) ):
        print('File '+REMAP_WEIGHTS+' exists, no recalculation of remapping weights')

for ncfile in ncfileList:
    OUTFILE=''.join(ncfile.split('.')[:-1])+'_remappedto_'+targetgrid.split('.')[0].split('/')[-1]+'.nc'
    print('Creating file '+OUTFILE)
    os.system('cdo remap,'+targetgrid+','+REMAP_WEIGHTS+' '+ncfile+' '+OUTFILE)
    os.system('rm '+ncfile)
