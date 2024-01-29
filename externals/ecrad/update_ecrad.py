#!/usr/bin/env python3
#
# Update the ecrad code in ICON with files from a branch at ECMWF
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

import traceback
import sys
import getopt
import subprocess

# Default settings (can be overwritten when adding options to script execution
BRANCH_GITECMWF="dwd-master"
WORKDIR_PRE="."
ICONDIR="."

# Settings carved in stone:
SSH_GITECMWF="ssh://git@git.ecmwf.int/escape/dwarf-p-radiation-ecrad.git"
TMPFOLDER="updateecrad"
MODSOURCESTRING = "This file has been modified for the use in ICON"
NOTICESTRING = "For the use in ICON, the following changes have been made to the original ecRad code:\n"+\
               "- OpenACC port for the use of ecRad on GPU architectures\n"+\
               "- Code optimizations to improve performance on vector architectures\n"+\
               "- Replace yomhook with ecradhook to avoid linking conflicts with other ECMWF software used by ICON"

#ACC port and optimizations in the ifsrrtm subfolder and ecradhook

def createCleanWorkdir(workdir):
    # Delete leftovers
    print("Deleting "+workdir+" ...")
    executeCommand("rm -rf "+workdir)
    # Create clean directory
    print("Creating working directory "+workdir+" ...")
    executeCommand("mkdir "+workdir)

def cloneCheckout_gitECMWF(workdir):
    executeCommand("cd "+workdir+" && git clone "+SSH_GITECMWF)
    tmpecradclonedir=workdir+"/dwarf-p-radiation-ecrad"
    executeCommand("cd "+tmpecradclonedir+" && git checkout -b "+BRANCH_GITECMWF+" origin/"+BRANCH_GITECMWF)
    return tmpecradclonedir

def executeCommand(cmd):
    subprocess.run(cmd,shell=True)

def copyFilesFromECMWF(workdir):
    print("Updating ecrad files in ICON...")
    executeCommand("cp "+workdir+"/radiation/*.F90 "+ICONDIR+"/radiation/.")
    executeCommand("cp "+workdir+"/radiation/*.h "+ICONDIR+"/radiation/.")
    executeCommand("cp "+workdir+"/include/*.h "+ICONDIR+"/include/.")
    executeCommand("cp "+workdir+"/ifsrrtm/*.F90 "+ICONDIR+"/ifsrrtm/.")
    executeCommand("cp "+workdir+"/ifsrrtm/AER-BSD3-LICENSE "+ICONDIR+"/ifsrrtm/AER-BSD3-LICENSE")
    executeCommand("cp "+workdir+"/ifsaux/*.F90 "+ICONDIR+"/ifsaux/.")
    executeCommand("cp "+workdir+"/utilities/*.F90 "+ICONDIR+"/utilities/.")
    executeCommand("cp "+workdir+"/data/*.nc "+ICONDIR+"/data/.")
    executeCommand("cp "+workdir+"/NOTICE "+ICONDIR+"/NOTICE")
    # Adjust yomhook to ecradhook
    executeCommand("cd "+ICONDIR+" && sed -i 's/yomhook/ecradhook/g' */*.F90")
    executeCommand("cd "+ICONDIR+" && sed -i 's/YOMHOOK/ecradhook/g' */*.F90")
    # Special rules
    executeCommand("rm "+ICONDIR+"/utilities/print_matrix.F90")

def cleanup(workdir):
    print("Cleaning up...")
    executeCommand("rm -rf "+workdir)

def getModifiedFileList(branch1, branch2, workdir):
    command = "cd "+workdir+" && "+"git diff --name-only "+branch1+" "+branch2
    modFiles = subprocess.run(command, capture_output=True,shell=True,text=True).stdout.splitlines()
    return(modFiles)

def prependStringToFile(commentSign, prependString, filename):
    with open(filename, "r") as f:
        fileContent = f.read()
    with open(filename, "w") as f:
        f.write(commentSign+" "+prependString+"\n\n")
        f.write(fileContent)

def appendStringToFile(appendString,filename):
    with open(filename, "a") as f:
        f.write("\n\n"+appendString)

try:
    opts, remainder = getopt.getopt(sys.argv[1:],"hb:w:i:",["help","branch=","workdir=","icondir="])
    for opt, arg in opts:
        if opt == '-h':
            print ('Example: update_ecrad.py -b <branchname> -w <workdir> -i <icondir>')
            print ('Full list of options:')
            print ('-h : help')
            print ('-b : name of source branch from ECMWF git')
            print ('-w : temporary work directory')
            print ('-i : ICON root directory')
            sys.exit()
        elif opt in ("-b", "--branch"):
            BRANCH_GITECMWF = arg
        elif opt in ("-w", "--workdir"):
            WORKDIR_PRE = arg
        elif opt in ("-i", "--icondir"):
            ICONDIR = arg+"/externals/ecrad"

    print('===Settings===')
    print('Source branch: '+BRANCH_GITECMWF)
    print('Working directory: '+WORKDIR_PRE)
    print('ecRad directory: '+ICONDIR)
    print('==============')

    tmpworkdir=WORKDIR_PRE+"/"+TMPFOLDER
    createCleanWorkdir(tmpworkdir)

    tmpecradclonedir = cloneCheckout_gitECMWF(tmpworkdir)

    # Get list of files that are modified compared to the original ECMWF version
    modFileList = getModifiedFileList('master', BRANCH_GITECMWF, tmpecradclonedir)

    copyFilesFromECMWF(tmpecradclonedir)

    f90list = [f for f in modFileList if '.F90' in f]
    for f90file in f90list:
        prependStringToFile('!', MODSOURCESTRING, ICONDIR+"/"+f90file)
    intfbhlist = [f for f in modFileList if '.intfb.h' in f]
    for intfbhfile in intfbhlist:
        prependStringToFile('!', MODSOURCESTRING, ICONDIR+"/"+intfbhfile)

    appendStringToFile(NOTICESTRING ,ICONDIR+"/NOTICE")

    cleanup(tmpworkdir)

except Exception as e:
    print(e)
    traceback.print_exc(e)
