#! /usr/bin/python

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# -*- coding: utf-8 -*-
#-----------------------------------------------------------------------------

import os
import sys
import glob
import re

plotYearsPerDay=1
plotIterationsCellPerSec=2

totalTimeIndex=6
totalCallIndex=1
minTimeIndex=2
meanTimeIndex=3
maxTimeIndex=4

def Die(arg):
  print arg
  exit(1)
  

def getFileTimer(fileName, timerName, counterIndex):
  
  #print "proccessing ", fileName, timerName, "..."
  timerValue=-1.0
  try:
    infile = open(fileName, 'r')
  except IOError:
    return timerValue
  #print timerName
  
  for line in infile:
    k=line.find(timerName)
    if (k < 0):
      continue
    # we found the timer
    #print line
    s = line[k:].split()
    #print s
    #n = len(s)
#    if (n != 7 and n != 8): continue
    s1=s[counterIndex].rstrip('s')
    try:
      timerValue = float(s[counterIndex].rstrip('s'))
    except ValueError:
      timerValue = -1.0
    
    #print fileName, timerName, s[counterIndex], timerValue
    
    break

  infile.close()
   
  return timerValue


def getTimerValues(folderName, timerName, counterIndex):
  
  meanTime = 0.0
  maxTime  = 0.0
  minTime  = 9999999999999.0
  totalCount = 0.0
  for fileName in glob.glob(folderName+'/timer*'):
    
    timerValue = getFileTimer(fileName, timerName, counterIndex)
    if (timerValue > 0.0 ):
      totalCount = totalCount + 1.0
      meanTime = meanTime + timerValue
      maxTime = max(maxTime, timerValue)
      minTime = min(minTime, timerValue)
      #print timerValue, maxTime

  if (totalCount > 0):
    meanTime = meanTime / totalCount
  else:
    meanTime = 0.0
    maxTime  = 0.0
    minTime  = 0.0
    
  print folderName, timerName, minTime, meanTime, maxTime
     
  return (minTime, meanTime, maxTime)


def sortFile(fileName):
  
  dataFileName=fileName+".dat"
  infile = open(dataFileName, 'r')
  lines = [line for line in infile]
  lines.sort(key=int)

def writeOldTimerProcsDistribution(outName, logFileName, timerName, counterIndex):
  
  outFileName=outName+".dat"
  outFile = open(outFileName, 'w')
  infile = open(logFileName, 'r')
  
  for line in infile:
    s = line.split()
    k=line.find(timerName)
    if (k < 0):
      continue
    # we found the timer
    splitAll = line.split()
    #print line
    #print splitAll
    s = line[k:].split()
    #n = len(s)
#    if (n != 7 and n != 8): continue
    s1=s[counterIndex].rstrip('s')
    try:
      timerValue = float(s[counterIndex].rstrip('s'))
    except ValueError:
      timerValue = -1.0
    
    outFile.write(splitAll[0]+' '+str(timerValue)+"\n")
  
  outFile.close()     
  return 

#logFile="LOG.exp.nat_ape_r2b6_1sendrecv_0initc_78levels_240nodes_8mpi_procs_8threads_11nproma.sdays_1-30.run.o.timers"
#logFile="./run/LOG.exp.nat_ape_r2b5N6.28nodes.896procs.run.o"
#writeOldTimerProcsDistribution("b5N6_28nodes_min_radiation", logFile, " nwp_radiation ", minTimeIndex)
#writeOldTimerProcsDistribution("b5N6_28nodes_mean_radiation", logFile, " nwp_radiation ", meanTimeIndex)
#writeOldTimerProcsDistribution("b5N6_28nodes_max_radiation", logFile, " nwp_radiation ", maxTimeIndex)
#sys.exit()

def writeTimerProcsDistribution(outName, folderName, timerFilePrefix, timerName, counterIndex):
  
  outFileName=outName+".dat"
  outFile = open(outFileName, 'w')
  noOfCores=len(glob.glob(folderName+'/'+timerFilePrefix+'*'))
  for i in range(noOfCores):
    fileName=folderName+'/'+timerFilePrefix+str(i).zfill(4)
   # print "proccessing file:"+fileName
    timerValue = getFileTimer(fileName, timerName, counterIndex)
    if (timerValue > 0.0 ):
      outFile.write(str(i)+' '+str(timerValue)+"\n")
  
  outFile.close()     
  return 


def getTimeMaxInFiles(folderName, timerFilePrefix, timerName, counterIndex):
  
  noOfCores=len(glob.glob(folderName+'/'+timerFilePrefix+'*'))
  maxTime=0.0
  maxProc=-1
  for i in range(noOfCores):
    fileName=folderName+'/'+timerFilePrefix+str(i).zfill(4)
   # print "proccessing file:"+fileName
    timerValue = getFileTimer(fileName, timerName, counterIndex)
    if (timerValue > maxTime ):
      maxTime = timerValue
      maxProc = i
  
  return (maxTime, maxProc)



def processFolderList(foldersList, timersList, counterIndex):
  
  noOfFolders = len(foldersList)
  noOfTimers  = len(timersList)
  
  minTimeList = [None]*noOfFolders
  meanTimeList = [None]*noOfFolders
  maxTimeList = [None]*noOfFolders
  for i in range(noOfFolders):
    minTimeList [i] = [0.0]*noOfTimers
    meanTimeList [i] = [0.0]*noOfTimers
    maxTimeList [i] = [0.0]*noOfTimers
     
  # get mean values
  for i in range(noOfFolders):
    for j in range(noOfTimers):
      print "proccessing ", foldersList[i], timersList[j], counterIndex, "..."
      (minTime, meanTime, maxTime) = getTimerValues(foldersList[i], timersList[j], counterIndex)
      # meanTimeList[i][j] = meanTime, use max time 
      minTimeList[i][j] = minTime
      meanTimeList[i][j] = meanTime
      maxTimeList[i][j] = maxTime
  
  return(minTimeList, meanTimeList, maxTimeList)
    
def writeSingleTimer(name, timer, coresList):
  noOfFolders = len(timeList)
  
  outFileName=name+".dat"
  outFile = open(outFileName, 'w')
  
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    outFile.write("  "+str(timer[i]) )
    outFile.write("\n")      
  outFile.close()

def writeParallelEfficency(name, timerList, timerIndex, coresList):
  noOfFolders = len(timerList)
  
  outFileName=name+"_parallel_efficency.dat"
  outFile = open(outFileName, 'w')

  referenceTime = timerList[0][timerIndex] * coresList[0]
  
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    outFile.write("  "+str( referenceTime / (timerList[i][timerIndex] * coresList[i])) )
    outFile.write("\n")      
  outFile.close()

def writeTimerList(name, timeList, coresList):
  noOfFolders = len(timeList)
  noOfTimers  = len(timeList[0])
  
  outFileName=name+".dat"
  outFile = open(outFileName, 'w')
  
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    for j in range(noOfTimers):
      outFile.write("  "+str(timeList[i][j]) )
    outFile.write("\n")      
  outFile.close()

def writeYearsPerDay(name, timeList, coresList, daysList, timerIndex):
  # write years per day
  noOfFolders = len(coresList)
  
  outFileName=name+"_years_per_day.dat"
  outFile = open(outFileName, 'w')

  for i in range(noOfFolders):
    
    print i, "cores:", coresList[i], " time:", timeList[i][0], " days:", daysList[i]
    #yearsPerDay=( daysList[i] * 24.0 * 3600.0 ) / (timeList[i][timerIndex] * 365.0)
    yearsPerDay=( daysList[i] * 236.712329 ) / (timeList[i][timerIndex])

    outFile.write(str(coresList[i])+" "+str(yearsPerDay)+"\n")

  outFile.close()
  

def writeStatsAllList(name, timeList, coresList, daysList, grid_size, plotScaling):
  # write times 
  noOfFolders = len(coresList)
  noOfTimers  = len(timeList[0])
  
  outFileName=name+"_times.dat"
  outFile = open(outFileName, 'w')
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    for j in range(noOfTimers):
      outFile.write("  "+str(timeList[i][j]) )
    outFile.write("\n")      
  outFile.close()

  # write times per iter/cell
  outFileName=name+"_times_per_cell_iter.dat"
  outFile = open(outFileName, 'w')
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    for j in range(noOfTimers):
      normalizedTime=timeList[i][j] / ( grid_size[i] * daysList[i] )
      outFile.write("  "+str(normalizedTime) )
    outFile.write("\n")      
  outFile.close()

  # write normalized times (cost per node, day, grid size
  outFileName=name+"_scaletimes_normalized.dat"
  outFile = open(outFileName, 'w')
  for i in range(noOfFolders):
    outFile.write(str(coresList[i]))
    for j in range(noOfTimers):
      normalizedTime=( (timeList[i][j] / daysList[i]) * coresList[i] ) / grid_size[i]
      outFile.write("  "+str(normalizedTime) )
    outFile.write("\n")      
  outFile.close()

  # write years per day
  if (plotScaling==plotYearsPerDay):
    outFileName=name+"_years_per_day.dat"
    outFile = open(outFileName, 'w')
  
    for i in range(noOfFolders):
      yearsPerDay=( daysList[i] * 24.0 * 3600.0 ) / (timeList[i][0] * 365.0)

      outFile.write(str(coresList[i])+" "+str(yearsPerDay)+"\n")

    outFile.close()
  
  # write years per day
  if (plotScaling==plotIterationsCellPerSec):
    outFileName=name+"_itercell_per_sec.dat"
    outFile = open(outFileName, 'w')
  
    for i in range(noOfFolders):
      iterPerDay=( daysList[i] * grid_size[i] ) / (timeList[i][0])

      outFile.write(str(coresList[i])+" "+str(iterPerDay)+"\n")

    outFile.close()
      
  return 


def getImbalance(minTimeList, maxTimeList, refTimer):
  # write times 
  noOfFolders = len(minTimeList)
  noOfTimers  = len(minTimeList[0])
  
  imbalanceList = [None]*noOfFolders
  imbalancePerCent = [None]*noOfFolders
  for i in range(noOfFolders):
    imbalanceList [i] = [0.0]*noOfTimers
    imbalancePerCent [i] = [0.0]*noOfTimers
  
  for i in range(noOfFolders):
    for j in range(noOfTimers):
      imbalanceList[i][j]=maxTimeList[i][j] - minTimeList[i][j]
      imbalancePerCent[i][j]=(imbalanceList[i][j] * 100.0) / minTimeList[i][refTimer]
      #print maxTimeList[i][j], minTimeList[i][j], minTimeList[i][refTimer], imbalancePerCent[i][j]

  return (imbalanceList, imbalancePerCent)

def normalizeTimerList(timeList, normValueIndex):
  noOfFolders = len(timeList)
  noOfTimers  = len(timeList[0])
  
  normTimeList = [None]*noOfFolders
  for i in range(noOfFolders):
    normTimeList [i] = [0.0]*noOfTimers
  
  for i in range(noOfFolders):
    for j in range(noOfTimers):
      normTimeList[i][j] = timeList[i][j] / timeList[i][normValueIndex]

  return normTimeList


def addTimers(timerList, timer1, timer2):
  noOfFolders = len(timerList)
  newTimerIndex=len(timerList[0])
  
  for i in range(noOfFolders):
    timerList[i].append(timerList[i][timer1] + timerList[i][timer2])
      
  return newTimerIndex

def subtractTimers(timerList, timer, minus_timer):
  noOfFolders = len(timerList)
  newTimerIndex=len(timerList[0])
  
  for i in range(noOfFolders):
    timerList[i].append(timerList[i][timer] - timerList[i][minus_timer])

  return newTimerIndex

def subtractFromTimer(timerList, timer, minus_timer):
  noOfFolders = len(timerList)  
  for i in range(noOfFolders):
    timerList[i][timer]=timerList[i][timer] - timerList[i][minus_timer]

def newTimer(timerList):
  noOfFolders = len(timerList)
  
  for i in range(noOfFolders):
    timerList[i].append(0.0)

  return len(timerList[0])-1

timer_total=0

timersListTotal=[
" total"
]


timerExchData=[
" exch_data"
]


timersCommunnication=[
#" total",          #2
" exch_data",      #3
" exch_data.wait", #4
" icon_comm_sync", #5
" comm_wait"       #6
]



def stdall_timers_statistics(folders_list, resolution, cores_list, days_list):
  
  (minTimes, meanTimes, maxTimes)=processFolderList(folders_list, timersListNat, totalTimeIndex )
  writeYearsPerDay(resolution, maxTimes, cores_list, days_list, timer_total)
  writeParallelEfficency(resolution, maxTimes, timer_total, cores_list)

  timer_totalcomm_index=addTimers(meanTimes, timer_exch_data, timer_iconcomm_index)
  timer_radcomm_index=subtractTimers(meanTimes, timer_nwp_radiation_index, timer_radiaton_comp_index)
  normTimeList=normalizeTimerList(meanTimes, timer_total)
  writeTimerList(resolution+"_normalizedCost", normTimeList, cores_list)
  
  #(minTimes, meanTimes, maxTimes)=processFolderList(folders_list, timersCommunnication, totalTimeIndex )
  #normTimeList=normalizeTimerList(meanTimes, timer_total)
  #writeTimerList(resolution+"_communication_times", meanTimes, cores_list)
  #writeTimerList(resolution+"_communication_normalizedCost", normTimeList, cores_list)
  
#=========================================================================
timersOcean=[
" total",          #2
" solve_ab",       #3
" gmres",          #4
" gmres_p_sum",    #5
" tracer",         #6
" veloc_prod",     #7
" upd_phys_param", #8
" upd_flx",        #9
" exch_data",      #10
" exch_data.wait", #11
]

timer_exch_data=8
timer_exch_data_wait=9

def getNodesFromFolder(folder):
  folderNumbers=re.findall(r'\d+', folder)
  return float(folderNumbers[2])

def ocean_timers_statistics(resolution, extent):
  
  folders_list=glob.glob('omip_'+resolution+'*'+str(threads)+'threads'+extent)
  folders_list.sort(key=getNodesFromFolder)
  print folders_list
  
  noOfFolders = len(folders_list)
  cores_list = [None] * noOfFolders
  days_list = [10.0] * noOfFolders
  
  for i in range(noOfFolders):
    # print folders_list[i].split('_')
    # folderNumbers=re.findall(r'\d+', folders_list[i])
    # print  folderNumbers
    # cores_list[i]=float(folderNumbers[2]) * 16
    cores_list[i]=float(getNodesFromFolder(folders_list[i])) * 24.0
    
  print cores_list
  prefix=resolution+"_"+str(threads)+"threads"
  
  (minTimes, meanTimes, maxTimes)=processFolderList(folders_list, timersOcean, totalTimeIndex )
  writeParallelEfficency(prefix, maxTimes, timer_total, cores_list)  
  writeYearsPerDay(prefix, maxTimes, cores_list, days_list, timer_total)
      
  normTimeList=normalizeTimerList(meanTimes, timer_total)
  writeTimerList(resolution+"_normalizedCost", normTimeList, cores_list)
  
  #writeTimerList(prefix+"_meantimes", meanTimes, cores_list)
  #writeTimerList(prefix+"_maxtimes", maxTimes, cores_list)
    
def ocean_timer_core_distribution(distibution_folder, timerName):

  timerFilePrefix="timer.ocean." # only needed for cores time distribution
  writeTimerProcsDistribution(timerName.strip()+"_min",   distibution_folder, timerFilePrefix, timerName, minTimeIndex)
  writeTimerProcsDistribution(timerName.strip()+"_mean",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
  writeTimerProcsDistribution(timerName.strip()+"_max",   distibution_folder, timerFilePrefix, timerName, maxTimeIndex)

def oceanAlltests():

  ocean_timers_statistics("10km", "")
  #ocean_timers_statistics("160km", "__1")

threads=12
#oceanAlltests()
ocean_timer_core_distribution("omip_10km_64levels_288nodes_12threads", "exch_data")
ocean_timer_core_distribution("omip_10km_64levels_288nodes_12threads", "tracer")

sys.exit()



#=========================================================================

timersListNat=[
" total",         #2
" nh_solve",      #3
" physics",       #4
" radiaton_comp", #5
" transport",     #6
" exch_data",     #7
" nh_solve.exch", #8
" icon_comm_sync", #9
" nwp_radiation" #10
]
#timer_totalcomm_index 11
#timer_radcomm_index   12
timer_exch_data=5
timer_nhsolveexch_index=6
timer_iconcomm_index=7
timer_nwp_radiation_index=8
timer_radiaton_comp_index=3

def blizzard_bechmarks_01_03_2013_scaling():

  days_list=[6]*20

  #---------------------------------------------------
  resolution="20km"
  folders_list=[
  "nat-ape_6days_iconR2B07-grid_96levels.1decm_0-comm.61nodes.0radsplit.16nproma.1952procs",
  "nat-ape_6days_iconR2B07-grid_96levels.1decm_0-comm.122nodes.0radsplit.16nproma.3904procs",
  #"nat-ape_6days_iconR2B07-grid_96levels.1decm_0-comm.246nodes.0radsplit.21nproma.7872procs",
  "nat-ape_6days_iconR2B07-grid_96levels.104decm_1-comm.246nodes.8radsplit.21nproma.7872procs"
  ]
  cores_list=[1952, 3904, 7872]
  stdall_timers_statistics(folders_list, resolution, cores_list, days_list)
  
  #---------------------------------------------------
  resolution="40km"
  folders_list=[
  "nat-ape_6days_iconR2B06-grid_96levels.1decm_0-comm.61nodes.0radsplit.12nproma.1952procs",
  "nat-ape_6days_iconR2B06-grid_96levels.1decm_0-comm.122nodes.0radsplit.12nproma.3904procs",
  "nat-ape_6days_iconR2B06-grid_96levels.1decm_0-comm.246nodes.0radsplit.21nproma.7872procs"
  #"nat-ape_6days_iconR2B06-grid_96levels.104decm_103-comm.246nodes.6radsplit.21nproma.7872procs",
  #"nat-ape_6days_iconR2B06-grid_96levels.104decm_1-comm.246nodes.6radsplit.21nproma.7872procs"
  ]
  cores_list=[1952, 3904, 7872]
  stdall_timers_statistics(folders_list, resolution, cores_list, days_list)

  
  #---------------------------------------------------
  resolution="26km"
  folders_list=[
  "nat-ape_6days_ico_0026km_hdec07682_dcells096_96levels.1decm_0-comm.61nodes.0radsplit.16nproma.1952procs",
  "nat-ape_6days_ico_0026km_hdec07682_dcells096_96levels.1decm_0-comm.122nodes.0radsplit.16nproma.3904procs",
  #"nat-ape_6days_ico_0026km_hdec07682_dcells096_96levels.1decm_0-comm.241nodes.0radsplit.24nproma.7712procs",
  "nat-ape_6days_ico_0026km_hdec07682_dcells096_96levels.0decm_1-comm.241nodes.0radsplit.24nproma.7682procs"
  #"nat-ape_6days_ico_0026km_hdec07682_dcells096_96levels.0decm_1-comm.241nodes.8radsplit.24nproma.7682procs"
  ]
  cores_list=[1952, 3904, 7682]
  stdall_timers_statistics(folders_list, resolution, cores_list, days_list)

  #---------------------------------------------------
  resolution="53km"
  folders_list=[
  "nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.1decm_0-comm.61nodes.0radsplit.16nproma.1952procs",
  "nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.1decm_0-comm.122nodes.0radsplit.16nproma.3904procs",
  #"nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.1decm_0-comm.241nodes.0radsplit.12nproma.7712procs",
  #"nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.0decm_103-comm.241nodes.4radsplit.12nproma.7682procs",
  #"nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.0decm_1-comm.241nodes.4radsplit.12nproma.7682procs",
  "nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.0decm_1-comm.241nodes.0radsplit.12nproma.7682procs"
  ]
  cores_list=[1952, 3904, 7682]
  stdall_timers_statistics(folders_list, resolution, cores_list, days_list)


  #distibution_folder="nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.1decm_0-comm.241nodes.0radsplit.12nproma.7712procs"
  #distibution_folder="nat-ape_6days_ico_0053km_hdec07682_dcells024_96levels.0decm_1-comm.241nodes.0radsplit.12nproma.7682procs"
  ##timerName=" radiaton_comp "
  #timerName=" nh_solve.p1 "

  #writeTimerProcsDistribution("53km_241nodes_nh_solve_p1_min",   distibution_folder, timerFilePrefix, timerName, minTimeIndex)
  #writeTimerProcsDistribution("53km_241nodes_nh_solve_p1_mean",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
  #writeTimerProcsDistribution("53km_241nodes_nh_solve_p1_max",   distibution_folder, timerFilePrefix, timerName, maxTimeIndex)

def thunder_first_tests():

  distibution_folder="nat-ape_barrier_1days_bind_0msg_mkl_iconR2B06-grid_96levels.1decm_0-comm.32nodes.0radsplit.16nproma.512procs"
  #timerName=" radiaton_comp "
  timerName=" nh_solve.p1 "

  writeTimerProcsDistribution(timerName.strip()+"_min",   distibution_folder, timerFilePrefix, timerName, minTimeIndex)
  writeTimerProcsDistribution(timerName.strip()+"_mean",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
  writeTimerProcsDistribution(timerName.strip()+"_max",   distibution_folder, timerFilePrefix, timerName, maxTimeIndex)


thunder_first_tests()
sys.exit()

blizzard_bechmarks_01_03_2013_scaling()






timersListEcham=[
" total",
" exch_data", 
" echam_phy"
]


timersListRK=[
" total",
" RK_solve",
" echam_phy",
" transport",
" radiation",
" gw_hines_opt",
" hdiff_expl",
" phy2dyn",
" exch_data" 
]

timersListSI=[
" total",
" 2tl_si_solve",
" echam_phy",
" transport",
" radiation",
" gw_hines_opt",
" hdiff_expl",
" phy2dyn",
" exch_data",
" omp_global_sum"
]


timersListNHsolver=[
" nh_solve"
]

timersListComm=[
 " sync_1_1_3dcells_1",
 " sync_1_1_3dcells_2",
 " sync_1_1_3dcells_3",
 " sync_1_1_3dcells_4",
 " sync_1_1_3dedges_1",
 " sync_1_1_3dedges_2",
 " sync_1_1_3dedges_3",
 " sync_1_1_3dedges_4"
]

r2b6_size=[327680]*7



timersList=[
" mpi_barrier",
" icon_comm_sync",
" comm_fillsend",
" comm_fillrecv",
" comm_ircv",
" comm_isend ",
" comm_wait",
" radiaton_recv "
]


#timersListNat=[
#" total",         #2
#" nh_solve",      #3
#" physics",       #4
#" nwp_radiation", #5
#" transport",     #6
#" exch_data",     #7
#" nh_solve.exch", #8
#" mpi_barrier",  #9
#" icon_comm_syn"  #10
#]


distibution_folder="nat-ape_10days_1decm_0-comm_4sendrecv.2threads.64nodes.2048procs.10nproma.iconR2B06-grid.96levels"
#timerName=" radiaton_comp "
timerName=" nwp_radiation "
timerFilePrefix="timer.atmo."

writeTimerProcsDistribution("64nodes_min_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, minTimeIndex)
writeTimerProcsDistribution("64nodes_mean_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
writeTimerProcsDistribution("64nodes_max_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, maxTimeIndex)

sys.exit()

#=========================================================================================================
# old experiments
#nat_days_list=[60]*20
##---------------------------------------------------------------
#nat_folders_list=[
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.2nodes.64procs.12nproma.iconR2B05-grid_dec-362.96levels",
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.4nodes.128procs.18nproma.iconR2B05-grid_dec-362.96levels"
#]
#nat_cores_list=[64, 128]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('140km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()


#---------------------------------------------------------------
# R2B04 160km, 64 levels
#nat_folders_list=[
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.2nodes.64procs.18nproma.iconR2B04-grid.64levels",
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.4nodes.128procs.12nproma.iconR2B04-grid.64levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.8nodes.256procs.16nproma.iconR2B04-grid.64levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.16nodes.512procs.12nproma.iconR2B04-grid.64levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.22nodes.704procs.8nproma.iconR2B04-grid.64levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.27nodes.864procs.12nproma.iconR2B04-grid.64levels"
#]
#nat_cores_list=[64, 128, 256, 512, 704, 864]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('160km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)


#=========================================================================================================
nat_days_list=[10]*20
#---------------------------------------------------------------
# R2B04 160km, 64 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.1nodes.32procs.16nproma.iconR2B04-grid.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.4nodes.128procs.16nproma.iconR2B04-grid.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.8nodes.256procs.20nproma.iconR2B04-grid.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.16nodes.512procs.10nproma.iconR2B04-grid.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.27nodes.864procs.12nproma.iconR2B04-grid.64levels"
#]
##nat_nodes_list=[32, 48, 64, 80]
#nat_cores_list=[32, 128, 256, 512, 864]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('160km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()

#---------------------------------------------------------------
# b5_dec-1082 140km, 64 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.1nodes.32procs.9nproma.iconR2B05-grid_dec-1082.64levels",
##"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.5nodes.160procs.9nproma.iconR2B05-grid_dec-1082.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.12nodes.384procs.9nproma.iconR2B05-grid_dec-1082.64levels",
##"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.23nodes.736procs.9nproma.iconR2B05-grid_dec-1082.64levels",
#"nat-ape_10days_hexdec_24cellspd_0decm_1-comm_4sendrecv.2threads.34nodes.1082procs.12nproma.iconR2B05-grid_dec-1082.64levels"
#]
##nat_cores_list=[32, 160, 384, 736, 1082]
#nat_cores_list=[32, 384, 1082]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('140km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()

#---------------------------------------------------------------
# R2B05_dec-362 122km, 64 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.2nodes.64procs.9nproma.iconR2B05-grid_dec-362.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.6nodes.192procs.9nproma.iconR2B05-grid_dec-362.64levels",
#"nat-ape_10days_hexdec_96cellspd_0decm_1-comm_4sendrecv.2threads.12nodes.362procs.12nproma.iconR2B05-grid_dec-362.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.27nodes.864procs.10nproma.iconR2B05-grid_dec-362.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.45nodes.1440procs.12nproma.iconR2B05-grid_dec-362.64levels"
#]

#nat_cores_list=[64, 192, 362, 864, 1440]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('122km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()

#---------------------------------------------------------------
# b5_dec-1922 105km, 64 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.12nodes.384procs.10nproma.iconR2B05-grid_dec-1922.64levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.24nodes.768procs.10nproma.iconR2B05-grid_dec-1922.64levels",
#"nat-ape_10days_hexdec_24cellspd_0decm_1-comm_4sendrecv.2threads.61nodes.1922procs.12nproma.iconR2B05-grid_dec-1922.64levels"
#]

#nat_cores_list=[384, 768, 1922]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('105km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()


#---------------------------------------------------------------
# R2B04 160km, 96 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.1nodes.32procs.16nproma.iconR2B04-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.4nodes.128procs.16nproma.iconR2B04-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.8nodes.256procs.20nproma.iconR2B04-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.16nodes.512procs.10nproma.iconR2B04-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.27nodes.864procs.12nproma.iconR2B04-grid.96levels"
#]
##nat_nodes_list=[32, 48, 64, 80]
#nat_cores_list=[32, 128, 256, 512, 864]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('160km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()

#---------------------------------------------------------------
# b5_dec-1082 140km, 96 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.1nodes.32procs.9nproma.iconR2B05-grid_dec-1082.96levels",
##"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.5nodes.160procs.9nproma.iconR2B05-grid_dec-1082.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.12nodes.384procs.9nproma.iconR2B05-grid_dec-1082.96levels",
##"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.23nodes.736procs.9nproma.iconR2B05-grid_dec-1082.96levels",
#"nat-ape_10days_hexdec_24cellspd_0decm_1-comm_4sendrecv.2threads.34nodes.1082procs.12nproma.iconR2B05-grid_dec-1082.96levels"
#]
##nat_cores_list=[32, 160, 384, 736, 1082]
#nat_cores_list=[32, 384, 1082]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('140km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()

#---------------------------------------------------------------
# R2B05_dec-362 122km, 96 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.2nodes.64procs.9nproma.iconR2B05-grid_dec-362.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.6nodes.192procs.9nproma.iconR2B05-grid_dec-362.96levels",
#"nat-ape_10days_hexdec_96cellspd_0decm_1-comm_4sendrecv.2threads.12nodes.362procs.12nproma.iconR2B05-grid_dec-362.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.27nodes.864procs.10nproma.iconR2B05-grid_dec-362.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.45nodes.1440procs.12nproma.iconR2B05-grid_dec-362.96levels"
#]

#nat_cores_list=[64, 192, 362, 864, 1440]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('122km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)

#distibution_folder="nat-ape_10days_104decm_1-comm_4sendrecv.2threads.2nodes.64procs.9nproma.iconR2B05-grid_dec-362.96levels"
#timerName=" radiaton_comp "
##timerName=" nwp_radiation "
#filePrefix="timer.atmo."

#writeTimerProcsDistribution("2nodes_min_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, minTimeIndex)
#writeTimerProcsDistribution("2nodes_mean_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
#writeTimerProcsDistribution("2nodes_max_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, maxTimeIndex)
#sys.exit()

#---------------------------------------------------------------
# b5_dec-1922 105km, 96 levels
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.6nodes.192procs.12nproma.iconR2B05-grid_dec-1922.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.12nodes.384procs.10nproma.iconR2B05-grid_dec-1922.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.24nodes.768procs.10nproma.iconR2B05-grid_dec-1922.96levels",
#"nat-ape_10days_hexdec_24cellspd_0decm_1-comm_4sendrecv.2threads.61nodes.1922procs.12nproma.iconR2B05-grid_dec-1922.96levels"
#]

#nat_cores_list=[192, 384, 768, 1922]
#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('105km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()


#---------------------------------------------------------------
# R2B06_dec-2432, 94km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.4nodes.128procs.12nproma.iconR2B06-grid_dec-2432.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.12nodes.384procs.11nproma.iconR2B06-grid_dec-2432.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.42nodes.1344procs.11nproma.iconR2B06-grid_dec-2432.96levels",
#"nat-ape_10days_hexdec_24cellspd_0decm_1-comm_4sendrecv.2threads.76nodes.2432procs.12nproma.iconR2B06-grid_dec-2432.96levels"
#]
#nat_cores_list=[ 128, 384, 1344, 2432 ]

#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('94km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)
#sys.exit()


#---------------------------------------------------------------
# R2B05, 79km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.10nodes.320procs.16nproma.iconR2B05-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.20nodes.640procs.16nproma.iconR2B05-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.40nodes.1280procs.16nproma.iconR2B05-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.80nodes.2560procs.8nproma.iconR2B05-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.107nodes.3424procs.12nproma.iconR2B05-grid.96levels"
#]
#nat_cores_list=[320, 640, 1280, 2560, 3424 ]

#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('80km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", maxTimes, timer_total, nat_cores_list)


# radiation distribution    
#distibution_folder="nat-ape_10days_104decm_1-comm_4sendrecv.2threads.32nodes.1024procs.16nproma.iconR2B06-grid.96levels"
#writeTimerProcsDistribution("32nodes_min_radiaton_comp",  distibution_folder, "timer.atmo.", " radiaton_comp ", minTimeIndex)
#writeTimerProcsDistribution("32nodes_mean_radiaton_comp", distibution_folder, "timer.atmo.", " radiaton_comp ", meanTimeIndex)
#writeTimerProcsDistribution("32nodes_max_radiaton_comp",  distibution_folder, "timer.atmo.", " radiaton_comp ", maxTimeIndex)
#sys.exit()
#---------------------------------------------------------------

#---------------------------------------------------------------
# R2B06-grid_dec-1442, 60km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.20nodes.640procs.12nproma.iconR2B06-grid_dec-1442.96levels",
#"nat-ape_10days_hexdec_96cellspd_0decm_1-comm_4sendrecv.2threads.46nodes.1442procs.12nproma.iconR2B06-grid_dec-1442.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.90nodes.2880procs.12nproma.iconR2B06-grid_dec-1442.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.180nodes.5760procs.12nproma.iconR2B06-grid_dec-1442.96levels"
#]
#nat_cores_list=[640, 1442, 2880, 5760 ]

#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListTotal, totalTimeIndex )
#writeYearsPerDay('60km', maxTimes, nat_cores_list, nat_days_list, 0)
#writeParallelEfficency("parallel_efficency", meanTimes, timer_total, nat_cores_list)
#sys.exit()
#---------------------------------------------------------------

#---------------------------------------------------------------
# R2B07-grid_dec-2432, 47 km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.38nodes.1216procs.12nproma.iconR2B07-grid_dec-2432.96levels",
#"nat-ape_10days_hexdec_96cellspd_0decm_1-comm_4sendrecv.2threads.76nodes.2432procs.12nproma.iconR2B07-grid_dec-2432.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.152nodes.4864procs.12nproma.iconR2B07-grid_dec-2432.96levels"
#]
#nat_cores_list=[1216, 2432, 4864 ]
#---------------------------------------------------------------

    
#---------------------------------------------------------------
# R2B06, 40km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.32nodes.1024procs.16nproma.iconR2B06-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.64nodes.2048procs.10nproma.iconR2B06-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.128nodes.4096procs.10nproma.iconR2B06-grid.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.244nodes.7808procs.11nproma.iconR2B06-grid.96levels"
#]
#nat_cores_list=[1024, 2048, 4096, 7808 ]

# radiation distribution    
#distibution_folder="nat-ape_10days_104decm_1-comm_4sendrecv.2threads.32nodes.1024procs.16nproma.iconR2B06-grid.96levels"
#writeTimerProcsDistribution("32nodes_min_radiaton_comp",  distibution_folder, "timer.atmo.", " radiaton_comp ", minTimeIndex)
#writeTimerProcsDistribution("32nodes_mean_radiaton_comp", distibution_folder, "timer.atmo.", " radiaton_comp ", meanTimeIndex)
#writeTimerProcsDistribution("32nodes_max_radiaton_comp",  distibution_folder, "timer.atmo.", " radiaton_comp ", maxTimeIndex)
#sys.exit()
#---------------------------------------------------------------

#---------------------------------------------------------------
# R2B07-grid_dec-4322, 35 km
#nat_folders_list=[
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.70nodes.2240procs.12nproma.iconR2B07-grid_dec-4322.96levels",
#"nat-ape_10days_hexdec_96cellspd_0decm_1-comm_4sendrecv.2threads.136nodes.4322procs.12nproma.iconR2B07-grid_dec-4322.96levels",
#"nat-ape_10days_104decm_1-comm_4sendrecv.2threads.240nodes.7680procs.9nproma.iconR2B07-grid_dec-4322.96levels"
#]
#nat_cores_list=[2240, 4322, 7680 ]

#writeTimerProcsDistribution("min-physics-distribution_70nodes",  nat_folders_list[0], "timer.atmo.", " physics ", minTimeIndex)
#writeTimerProcsDistribution("mean-physics-distribution_70nodes", nat_folders_list[0], "timer.atmo.", " physics ", meanTimeIndex)
#writeTimerProcsDistribution("max-physics-distribution_70nodes",  nat_folders_list[0], "timer.atmo.", " physics ", maxTimeIndex)
#writeTimerProcsDistribution("min-physics-distribution_136nodes",  nat_folders_list[1], "timer.atmo.", " physics ", minTimeIndex)
#writeTimerProcsDistribution("mean-physics-distribution_136nodes", nat_folders_list[1], "timer.atmo.", " physics ", meanTimeIndex)
#writeTimerProcsDistribution("max-physics-distribution_136nodes",  nat_folders_list[1], "timer.atmo.", " physics ", maxTimeIndex)
#writeTimerProcsDistribution("min-physics-distribution_240nodes",  nat_folders_list[2], "timer.atmo.", " physics ", minTimeIndex)
#writeTimerProcsDistribution("mean-physics-distribution_240nodes", nat_folders_list[2], "timer.atmo.", " physics ", meanTimeIndex)
#writeTimerProcsDistribution("max-physics-distribution_240nodes",  nat_folders_list[2], "timer.atmo.", " physics ", maxTimeIndex)

#(minTimes, meanTimes, maxTimes)=processFolderList(nat_folders_list, timersListNat, totalTimeIndex )

#timer_totalcomm_index=addTimers(meanTimes, timer_exch_data, timer_iconcomm_index)
#timer_radcomm_index=subtractTimers(meanTimes, timer_nwp_radiation_index, timer_radiaton_comp_index)
#normTimeList=normalizeTimerList(meanTimes, timer_total)
#writeTimerList("35km.normalizedCost", normTimeList, nat_cores_list)
#sys.exit()
#---------------------------------------------------------------


#distibution_folder="icon-testbed_jitter_control_20calc.246nodes.7872dec.32nodeprocs.2threads.4nproma.nogrid"
#timerName=" calculate "
#filePrefix="timer.testbed."
distibution_folder="nat-ape_10days_104decm_1-comm_4sendrecv.2threads.2nodes.64procs.9nproma.iconR2B05-grid_dec-362.64levels"
timerName=" radiaton_comp "
#timerName=" nwp_turbulence "
timerFilePrefix="timer.atmo."

(mTime, mProc)=getTimeMaxInFiles(distibution_folder, timerFilePrefix, timerName, minTimeIndex)
print "Max of min ", timerName, ": ", mTime, mProc
(mTime, mProc)=getTimeMaxInFiles(distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
print "Max of mean ", timerName, ": ", mTime, mProc
(mTime, mProc)=getTimeMaxInFiles(distibution_folder, timerFilePrefix, timerName, maxTimeIndex)
print "Max of max ", timerName, ": ", mTime, mProc

writeTimerProcsDistribution("2nodes_min_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, minTimeIndex)
writeTimerProcsDistribution("2nodes_mean_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, meanTimeIndex)
writeTimerProcsDistribution("2nodes_max_radiaton_comp",  distibution_folder, timerFilePrefix, timerName, maxTimeIndex)
 

sys.exit()




timer_totalcomm_index=addTimers(meanTimes, timer_exch_data, timer_iconcomm_index)
timer_remaincomm_index=subtractTimers(meanTimes, timer_exch_data, timer_nhsolveexch_index)

normTimeList=normalizeTimerList(meanTimes, timer_total)
writeTimerList("r2b5_dec-1082.x3.normalizedCost", normTimeList, nat_cores_list)

#writeYearsPerDay('r2b4_96levels', maxTimes, nat_cores_list, nat_days_list, 0)
sys.exit()


grid_size_24cells=coresList_24cells * 96


coresList_strong_scale=[32, 64, 128, 256, 512, 1024]
strong_scale_r2b5_grid_size=[181920]*7


(min_minTimes, mean_minTimes, max_minTimes)=processFolderList(   FoldersList_strong_scale,  timersListComm, minTimeIndex )
(min_meanTimes, mean_meanTimes, max_meanTimes)=processFolderList(FoldersList_strong_scale,  timersListComm, meanTimeIndex )
(min_maxTimes, mean_maxTimes, max_maxTimes)=processFolderList(   FoldersList_strong_scale,  timersListComm, maxTimeIndex )

writeTimerList("mean_minTime",  mean_minTimes,  coresList_strong_scale)
writeTimerList("mean_meanTime", mean_meanTimes, coresList_strong_scale)
writeTimerList("max_meanTime",  max_meanTimes,  coresList_strong_scale)
writeTimerList("mean_maxTime",  mean_maxTimes,  coresList_strong_scale)

sys.exit()

#nodesList_hex_dec=[41, 68]
#nodesList_hex_dec=[162, 642, 1082]
coresList_hex_dec=[64, 184, 324, 1284, 2164]
hex_dec_grid_size=[2880, 8640, 15360, 61440, 103680]
#iter_hex_dec=[16, 16, 16]
iter_hex_dec=[5000]*7

hex_dec_hires_grid_size=[61440, 245760, 414720]
days_hex_dec_hires=[10, 10, 10]


coresList_dsl_strong_scale=[32, 64, 128, 192]
iter_dsl_strong_scale=[1000]*4
r2b4_grid_size=[20480]*4
FoldersList_no_dsl_B4_strong_scale=[
"exp.nat_ape-B4_strong_scale.16.1nodes.4threads.4nroma",
"exp.nat_ape-B4_strong_scale.32.2nodes.4threads.4nroma",
"exp.nat_ape-B4_strong_scale.64.4nodes.4threads.4nroma",
"exp.nat_ape-B4_strong_scale.96.6nodes.4threads.4nroma"
]
FoldersList_dsl_B4_strong_scale=[
"exp.nat_ape-B4_strong_scale.16.1nodes.4threads.32nroma",
"exp.nat_ape-B4_strong_scale.32.2nodes.4threads.32nroma",
"exp.nat_ape-B4_strong_scale.64.4nodes.4threads.32nroma",
"exp.nat_ape-B4_strong_scale.96.6nodes.4threads.32nroma"
]

FoldersList_dsl_weak_scale=[
"exp.nat_ape-hexdec_lowres.32.2nodes.4threads.12nroma",
"exp.nat_ape-hexdec_lowres.92.6nodes.4threads.12nroma",
"exp.nat_ape-hexdec_lowres.162.11nodes.4threads.12nroma",
"exp.nat_ape-hexdec_lowres.642.41nodes.4threads.12nroma",
"exp.nat_ape-hexdec_lowres.1082.68nodes.4threads.4nroma"
]


# weak scaling for 24 cells




#(min_minTimes, mean_minTimes, max_minTimes)=processFolderList(FoldersList_weak_scale_24cells, timersListNHsolver, minTimeIndex )
#(min_meanTimes, mean_meanTimes, max_meanTimes)=processFolderList(FoldersList_weak_scale_24cells, timersListNHsolver, meanTimeIndex )
#(min_maxTimes, mean_maxTimes, max_maxTimes)=processFolderList(FoldersList_weak_scale_24cells, timersListNHsolver, maxTimeIndex )

#writeTimerList("max_min_time", max_minTimes, coresList_24cells)
#writeTimerList("mean_mean_time", mean_meanTimes, coresList_24cells)
#writeTimerList("max_mean_time", max_meanTimes, coresList_24cells)
#writeTimerList("mean_maxTimes", mean_maxTimes, coresList_24cells)

#sys.exit()

#(minTimes, meanTimes, maxTimes)=processFolderList(FoldersList_weak_scale_24cells, timersListNat, totalTimeIndex )
#writeStatsAllList("nat_ape_weak_24cells",
  #maxTimes, coresList_24cells, days_24cells, grid_size_24cells, plotYearsPerDay)

#sys.exit()


#nat_folders_list=[
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.2nodes.64procs.18nproma.iconR2B04-grid.96levels",
#"nat-ape_80km_60days_3-comm_4sendrecv.2threads.4nodes.128procs.12nproma.iconR2B04-grid.96levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.8nodes.256procs.16nproma.iconR2B04-grid.96levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.16nodes.512procs.12nproma.iconR2B04-grid.96levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.22nodes.704procs.8nproma.iconR2B04-grid.96levels",
#"nat-ape_80km_60days_0-comm_4sendrecv.2threads.27nodes.864procs.12nproma.iconR2B04-grid.96levels"
#]
#nat_cores_list=[64, 128, 256, 512, 704, 864 ]
#nat_cores_list=[128, 256, 512, 704, 864 ]
#distibution_folder="nat-ape_80km_10days_8lonlatredrad_104decm_1-comm_4sendrecv.2threads.27nodes.864procs.12nproma.iconR2B04-grid.96levels"

#writeTimerProcsDistribution("min-rad-distribution_27nodes",  distibution_folder, "timer.atmo.", " radiaton_comp ", minTimeIndex)
#writeTimerProcsDistribution("mean-rad-distribution_27nodest", distibution_folder, "timer.atmo.", " radiaton_comp ", meanTimeIndex)
#writeTimerProcsDistribution("max-rad-distribution_27nodes",  distibution_folder, "timer.atmo.", " radiaton_comp ", maxTimeIndex)
#sys.exit()
#---------------------------------------------------------------

#---------------------------------------------------------------
# R2B05-grid_dec-362, 121km
#timer_totalcomm_index=addTimers(meanTimes, timer_exch_data, timer_iconcomm_index) #11
#normTimeList=normalizeTimerList(meanTimes, timer_total)
#writeTimerList("r2b5_dec-362.normalizedCost", normTimeList, nat_cores_list)
#sys.exit()
#---------------------------------------------------------------



FoldersList_jitter=[
"timers.32cores",
"timers.64cores",
"timers.128cores",
"timers.256cores",
"timers.512cores",
"timers.1024cores"
]

coresList_jitter=[
32,
64,
128,
256,
512,
1024
]

timersList_jitter=[
" calculate ",
" mpi_barrier ",
" mpi_barrier_init"
]


calculateTimer=0



(min_minTimes, min_meanTimes, min_maxTimes)=processFolderList(FoldersList_jitter, timersList_jitter, minTimeIndex )
(mean_minTimes, mean_meanTimes, mean_maxTimes)=processFolderList(FoldersList_jitter, timersList_jitter, meanTimeIndex )
(max_minTimes, max_meanTimes, max_maxTimes)=processFolderList(FoldersList_jitter, timersList_jitter, maxTimeIndex )
#(imbalanceTimes, imbalancePerCent)=getImbalance(min_meanTimes,mean_meanTimes, calculateTimer)
#writeTimerList("jitter_time", imbalanceTimes, coresList_jitter)
#writeTimerList("jitter_percent", imbalancePerCent, coresList_jitter)
writeTimerList("minTimes", min_meanTimes, coresList_jitter)
writeTimerList("meanTimes", mean_meanTimes, coresList_jitter)

writeTimerProcsDistribution("minCalculation", "timers", "timer.", " calculate ", minTimeIndex)
writeTimerProcsDistribution("meanCalculation", "timers", "timer.", " calculate ", meanTimeIndex)
writeTimerProcsDistribution("maxCalculation", "timers", "timer.", " calculate ", maxTimeIndex)

sys.exit()

#writeTimerList("max_barrier_l", max_maxTimes, coresList_jitter)

#(maxImbalanceTimes, maxImbalancePerCent)=getImbalance(min_minTimes,max_meanTimes, calculateTimer)
#writeTimerList("heavy_maxImbalance", maxImbalancePerCent, coresList_jitter)


#writeTimerProcsDistribution("minCalculation_512", "timers.512cores", "timer.testbed.", " calculate ", minTimeIndex)
#writeTimerProcsDistribution("meanCalculation_512", "timers.512cores", "timer.testbed.", " calculate ", meanTimeIndex)
#writeTimerProcsDistribution("maxCalculation_512", "timers.512cores", "timer.testbed.", " calculate ", maxTimeIndex)
#writeTimerProcsDistribution("minBarrier_512", "timers.512cores", "timer.testbed.", " mpi_barrier ", minTimeIndex)
#writeTimerProcsDistribution("meanBarrier_512", "timers.512cores", "timer.testbed.", " mpi_barrier ", meanTimeIndex)
#writeTimerProcsDistribution("maxBarrier_512", "timers.512cores", "timer.testbed.", " mpi_barrier ", maxTimeIndex)


#processFolderList("testbed_nh-solver_dsl_weak",
#  FoldersList_dsl_weak_scale, coresList_hex_dec, iter_hex_dec, timersListNHsolver, totalTimeIndex, hex_dec_grid_size, plotIterationsCellPerSec)


#processFolderList("testbed_nh-solver_strong_scale",
#  FoldersList_strong_scale, coresList_strong_scale, iter_hex_dec, timersListNHsolver, totalTimeIndex, strong_scale_grid_size, plotIterationsCellPerSec)

#processFolderList("testbed_nh-solver",
#  FoldersList_hex_dec, coresList_hex_dec, iter_hex_dec, timersListNHsolver, totalTimeIndex, hex_dec_grid_size, plotIterationsCellPerSec)

#processFolderList("comm",
#  FoldersList_hex_comm, coresList_hex_dec, iter_hex_dec, timersListComm, totalTimeIndex, hex_dec_grid_size, plotYearScaling )

#processFolderList("nat_ape_weak",
 # FoldersList_hex_dec, coresList_hex_dec, days_hex_dec_hires, timersListNat, totalTimeIndex, hex_dec_hires_grid_size, plotYearScaling )

#processFolderList("nat.dycore.240",
  #FoldersList_nat_240, nodesList_240, days_nat_240, timersListNat, totalTimeIndex, r2b6_size, plotYearScaling )

#processFolderList("nat.dycore.120",
  #FoldersList_nat_120, nodesList_120, days_nat_240, timersListNat, totalTimeIndex, plotYearScaling )

#processFolderList("nat_rrad.16mpi_procs_4threads",
#  FoldersList_nat_rrad_16mpi_procs_4threads, nodesList_120, days_base_4, timersList, totalTimeIndex, plotYearScaling )

#processFolderList("ssprk_5_4.16mpi_procs_4threads",
  #FoldersList_ssprk_16mpi_procs_4threads, nodesList_120, days_base_1, timersListRK, totalTimeIndex,  r2b6_size, plotYearScaling )

#processFolderList("si.8mpi_procs_8threads",
  #FoldersList_si_8mpi_procs_8threads, nodesList_120, days_base_2, timersListSI, totalTimeIndex,  r2b6_size,  plotYearScaling)

#processFolderList("si.16mpi_procs_4threads",
  #FoldersList_si_16mpi_procs_4threads, nodesList_60, days_base_5, timerTotal, totalTimeIndex, plotYearScaling )

#processFolderList("nat.16mpi_procs_4threads",
  #FoldersList_nat_16mpi_procs_4threads, nodesList_default,  days_base_2, timerTotal, totalTimeIndex, plotYearScaling )
 

#processFolderList("ssprk_5_4.16mpi_procs_4threads",
  #FoldersList_ssprk_16mpi_procs_4threads, nodesList_default, days_base_2, timerTotal, totalTimeIndex, plotYearScaling )

#processFolderList("ssprk_5_4.16mpi_procs_4threads",
  #FoldersList_ssprk_16mpi_procs_4threads, nodesList_default, days_base_2, timersListEcham, totalTimeIndex )

#processFolderList("ssprk_5_4.8mpi_procs_8hreads",
  #FoldersList_ssprk_8mpi_procs_8threads, nodesList_default, days_base_5, timersList, totalTimeIndex )

#processFolderList("si.8mpi_procs_8threads",
  #FoldersList_si_8mpi_procs_8threads, nodesList_default, days_base_2, timersList, totalTimeIndex )

#processFolderList("si.16mpi_procs_4threads",
  #FoldersList_si_16mpi_procs_4threads, nodesList_60, days_base_5, timersList, totalTimeIndex )
 
#processFolderList("nat.16mpi_procs_4threads",
  #FoldersList_nat_16mpi_procs_4threads, nodesList_default,  days_base_2, timersList, totalTimeIndex )
 
#processFolderList("compareCalls",
  #FoldersList_compare_calls, nodesList_compare, days_compare, timerExchData, totalCallIndex )

