#! /usr/bin/env python3

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

# -*- coding: utf-8 -*-
#==============================================================================
# model paths classes
#==============================================================================
import os
import glob
import shutil
import sys
#-----------------------------------------------------------------------
class model_paths(object):

  def __init__(self):
    callPath=os.path.dirname(sys.argv[0])
    if not (callPath == "." or callPath == ""):
      self.thisPath   = os.path.dirname(os.path.realpath(__file__))
    else:
      self.thisPath   = os.getcwd()
    #print(callPath, self.thisPath)
    splitPath       = os.path.split(self.thisPath)
    splitPath       = os.path.split(splitPath[0])
    self.basePath   = splitPath[0]
    self.runPath    = self.basePath+"/run"
    self.experimentsListPath = self.thisPath+"/experiment_lists"

  def get_experimentsNames_inPaths(self, pathsInRun):
    os.chdir(self.runPath)
    experimentsNames = []
    for path in pathsInRun:
      experiment = glob.glob(path)
      if not experiment:
        print("No experiment "+path+" found. Stop")
        quit(1)
      experimentsNames.extend(experiment)
    os.chdir(self.thisPath)
    return experimentsNames

  def get_runpath(self):
    return self.runPath

  def thisExperimentExists(self, experimentName):
    return os.path.isfile(self.runPath+"/"+experimentName)

  def get_thisListPath(self, listName):
    return self.experimentsListPath+"/"+listName

  def thisListExists(self, listName):
    return os.path.isfile(self.get_thisListPath(listName))

  def deleteThisList(self, listName):
    if not self.thisListExists(listName):
      print("The list "+listName+" does not exist.")
      quit(1)
    os.remove(self.get_thisListPath(listName))
    
  def copyList(self, fromlist, tolist):
    if not self.thisListExists(fromlist):
      print("The list "+fromlist+" does not exist.")
      quit(1)
    if self.thisListExists(tolist):
      print("The list "+tolist+" exists. Please remove it first.")
      quit(1)
    shutil.copy(self.get_thisListPath(fromlist), self.get_thisListPath(tolist))
    if not self.thisListExists(tolist):
      print("Copy failed")
      quit(1)
    

  def print_paths(self):
    print("Base path:"+self.basePath)
    print("Run path:"+self.runPath)
    print("This path:"+self.thisPath)

  def getPathAndName(self, PathName):
    dirName  = os.path.dirname(PathName)
    fileName = os.path.basename(PathName)
    return dirName, fileName
    
paths=model_paths()

