#!/usr/bin/env python

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

import os,sys,glob,shutil,json,subprocess,multiprocessing
from cdo import *
from itertools import chain
import matplotlib
matplotlib.use('Agg')
import pylab,numpy,dateutil.parser
import re


# ==============================================================================
# USAGE {{{
def usage():
  return """
  # USAGE =====================================================================
  #
  #   ./archive.py EXP=<expname> FILEPATTERN=<pattern> ARCHDIR=<dir> TAG=<rxxxxx>
  #
  # ===========================================================================
  #
  # EXPNAME     : string for naming things
  # TAG         : additional placeholder, e.g. for revisions, if not used, set it to ''
  # FILEPATTERN : quoted wildcard for model result relative to the 'experiments' dir
  # ARCHDIR     : directory, where the archived data is placed (default: ./archive)
  #
  # more options are (keys area llways uppercase)
  # DEBUG       : False/True
  # CALCPSI     : path, where to find the psi calculation and plotting scipt
  # ICONPLOT    : call for icon_plot.ncl, default is 'nclsh /pool/data/ICON/tools/icon_plot.ncl -altLibDir=/pool/data/ICON/tools'
  #
  # #
  # !!!! THIS SCRIPT IS MEANT TO BE STARTED IN THE MODEL RESULTS DIRECTORY !!!!
  """
# }}} ----------------------------------------------------------------------------------
# OPTION HANDLING {{{ =========================================================
def parseOptions():
  options = {
              'GRID'        : 'global',                    # shortcut for the used gridtype(global,box,channel)
              'ARCHDIR'     : './archive',
              'PLOTDIR'     : './plots',
              'PLOTYEAR'    : None,
              'DOCTYPE'     : 'pdf',                       # target document format
              'OFORMAT'     : 'png',                       # target document format
              'EXP'         : 'oce_mpiom',                 # default experiment name
              'FILEPATTERN' : 'oce_mpiom/oce_mpiom_*.nc*', # default output file pattern
              'DEBUG'       : False,                       # debugging is switched of by default
              'FORCE'       : False,                       # recomputation of verything is switched off by default
                              # the psi processor/plotter
              'CALCPSI'     : '../../scripts/postprocessing/tools/calc_psi.py',
              'TAG'         : '',                          # addition revision information
              'ICONPLOT'    : 'nclsh ../../scripts/postprocessing/tools/icon_plot.ncl -altLibDir=../../scripts/postprocessing/tools -remapOperator=remapycon',
              'PROCS'       : 8,                           # number of threads/procs to be used for parallel operations
              'JOBISRUNNING': True,                        # avoid the last output file/result year by default
              # optional stuff
              'DRYRUN'      : False,                       # with this set to true, the model output is scanned for containing years, only
              'MOCPATTERN'  : 'MOC.*',
              'MOCPLOTTER'  : '../../scripts/postprocessing/tools/calc_moc.ksh',
              # options to select special parts od the script
              #'ACTIONS'     : 'archive,preproc,procMoc,plotMoc,procRegio,plotRegio,plotTf,plotHorz,plotX,procTSR,plotTSR,plotPsi,procIce,plotIce',#finalDoc',
               'ACTIONS'     : 'archive,preproc,procMoc,plotMoc,procRegio,plotRegio,plotTf,plotHorz,plotX,procTSR,plotTSR,plotPsi,procIce,plotIce'.split(','),
             }

  plotActions = 'procMoc,plotMoc,plotX,procTSR,plotTSR,plotPsi'.split(',')
  optsGiven = sys.argv[1:]
  for optVal in optsGiven:
    if ( 0 <= optVal.find('=') ):
      if (optVal.find('=') == optVal.rfind('=')):
        key,value = optVal.split('=')
      else:
        firstSeperatorIndex = optVal.find('=')
        key   = optVal[:(optVal.find('='))]
        value = optVal[firstSeperatorIndex+1:]

      if key in ['FORCE','DEBUG','JOBISRUNNING','DRYRUN']:
        value = value.lower() in ['true','1']
      if 'PROCS' == key:
        value = int(value)
      if 'ACTIONS' == key:
        value = value.split(',')

      options[key] = value
    else:
      print("Cannot parse options '%s'!"%(optVal))
      sys.exit(1)

  # switch i plotting, if PLOTYEAR is given
  if None != options['PLOTYEAR']:
    for plotAction in plotActions:
      if not plotAction in options['ACTIONS']:
        options['ACTIONS'].append(plotAction)

  return options
# }}} ----------------------------------------------------------------------------------
# HELPER METHODS {{{ =========================================================== 
def dbg(obj):
  if options['DEBUG']:
    print(obj)

""" save internal log """
def dumpLog():
  if 0 < len(LOG):
    with open(LOGFILE,"w") as f:
      f.write(json.dumps(LOG))

def add4Cleanup(files):
  if not 'removeThem' in LOG.keys():
    LOG['removeThem'] = []

  map(lambda x : LOG['removeThem'].append(x),files)

""" remove temp files """
def cleanUp():
  dbg(LOG['removeThem'])
  map(lambda x : os.remove(str(x)) if os.path.exists(str(x)) else None, LOG['removeThem'])
  LOG['removeThem'] = []

""" save internal log and exit """
def doExit(value=0):
  # save the restults to omitt reprocessing input files on subsequent calls
  dumpLog()
  cleanUp()
  sys.exit(value)

""" load the last archiving status ir available """
def loadLog(options):
  if options['FORCE']:
    LOG = {}
    if os.path.exists(LOGFILE):
      os.remove(LOGFILE)
  if os.path.exists(LOGFILE):
    dbg("Load last status from LOGFILE")
    LOG = json.load(open(LOGFILE))
  else:
    dbg("Could not find file logile: '#{@logile}'! Create a new one ....")
    options['FORCE'] = True
    LOG = {}

  for key in ['packedYears','removeThem']:
    if not key in LOG.keys():
      LOG[key] = []

  return LOG

""" show output names from the LOG """
def showLogEntries(log,keys):
  for key in keys:
    if key in log.keys():
      print(':'.join([key,log[key]]))

""" append list of image together """ #{{{
def collectImageToMapByRows(images,columns,ofile):
  imageMap     = {}
  imageColumns = {}
  allTempFiles = []
  for i in range(0,columns):
    imageColumns[i] = []

  for i,image in enumerate(images):
    imageColumns[i%columns].append(image)
  for column in imageColumns:
    _tfile  = tempfile.NamedTemporaryFile(delete=True,prefix='IconArchive')
    allTempFiles.append(_tfile)
    rowFile = _tfile.name+'.png'
    cmd     = 'convert +append %s %s'%(' '.join(imageColumns[column]),rowFile)
    dbg(cmd)
    subprocess.check_call(cmd,shell=True,env=os.environ)
    imageMap[column] = rowFile

  cmd = "convert -append %s %s"%(' '.join(imageMap.values()),ofile)
  dbg(cmd)
  subprocess.check_call(cmd,shell=True,env=os.environ)
  map(lambda x: os.remove(x),imageMap.values()) #}}}

""" warpper around cdo commands for multiprocessing """ #{{{
def showyear(file):
  return map(lambda x : int(x) ,cdo.showyear(input = file)[0].split(' '))

def ntime(file):
  return cdo.ntime(input = file)[0]

""" wrapper of 'cdo.yearmean' """
def yearmean(ifile,ofile,forceOutput):
  cdo.yearmean(input = ifile,
               output = ofile,
               forceOutput = forceOutput,
               options = '-f nc4c -z zip_1')
# }}}

""" collect the years/ntime, which aree stored in to files """ #{{{
def scanFilesForTheirYears(fileList,procs,log):

  if not log.has_key('yearsOfFiles'):
    log['yearsOfFiles'] = {}

  pool = multiprocessing.Pool(procs)
  results = []


  for ifile in fileList:
    years = pool.apply_async(showyear,[str(ifile)])
    results.append([ifile,years])

  for result in results:
      f = result[0]
      r = result[1]
      log['yearsOfFiles'][f] = r.get()

  pool.close()
  pool.join()

  dumpLog()

def scanFilesForTheirNTime(fileList,procs,log):

  pool = multiprocessing.Pool(procs)
  results = []

  for ifile in fileList:
    ntimes = pool.apply_async(ntime,[ifile])
    results.append([ifile,ntimes])

  for result in results:
      f = result[0]
      r = result[1]
      log[f] = r.get()

  pool.close()
  pool.join()
#}}}
""" create yearly and yearMean filenames """
def getFileNamesForYears(year,archdir,experimentInfo):
  yearFile     = "%s/%s_%s.nc"%(archdir,experimentInfo,year)
  yearMeanFile = "%s/%s_%s_mean.nc"%(archdir,experimentInfo,year)

  return [yearFile, yearMeanFile]

""" input processing for computing regio mean values """
def _computeRegioMean(depth,varname,ifile,myMask,ofile,newData):
  # mask out region
  # vertical interpolation to target depth
  # mean value computaion
  cdo.fldmean(input = '-div -sellevel,{0} -selname,{1} {2} {3}'.format(depth,varname,ifile,myMask),
              output = ofile,
              forceOutput = newData)

""" compute global mean of bias for t s rho """
def globalTempSalRho1D(ifile,initfile,maskfile,varNames,ofile):
  cdo.timmean(input = '-fldmean -div -sub -selname,{0} {1} {2} {3}'.format(varNames,ifile,initfile,maskfile),
              output = ofile)

""" compute horizontal last 30-year bias for t s rho """
def globalTempSalRho2D(ifile,initfile,varNames,archdir):
  maskName = 'wet_c'
  ofile    = '/'.join([archdir,'TSR_2D_%s'%(os.path.basename(ifile))])
  cdo.timmean(input = '-div -sub -selname,{0} -sellevidx,1 {1} -selname,{2} -sellevidx,1 {3} -sellevidx,1 -selname,{4} -seltimestep,1 {5}'.format(varNames,ifile,varNames,initfile,maskName,ifile),
              output = ofile)

""" write a single file for a given years from a list of input files """ #{{{
def grepYear(ifiles,year,archdir,forceOutput,experimentInfo):
  yearFile, yearMeanFile = getFileNamesForYears(year,archdir,experimentInfo)

  if (not os.path.exists(yearFile) or forceOutput):
    if (os.path.exists(yearFile)):
      os.remove(yearFile)

    catFiles = []
    for ifile in ifiles:
      dbg(ifile)
      ofile = "{0}/_catFile_{1}_{2}".format(archdir,year,os.path.basename(ifile))
      catFiles.append(cdo.selyear(year, input = '{0}'.format(ifile), output = ofile))

    dbg(['catFiles'] + catFiles)

    cdo.cat(input = ' '.join(catFiles), output = yearFile, options = '-f nc4c -z zip_1' if cdo.version() == '1.6.9' else '')

    map(lambda x: os.remove(x),catFiles)
  else:
    print("Use existing yearFile '{0}'".format(yearFile))
#}}}

""" split input files into yearly files """
def splitFilesIntoYears(plotYears,filesOfYears,archdir,forceOutput,expInfo,procs):
  pool      = multiprocessing.Pool(procs)
  yearFiles = []
  for year in plotYears:
    files = filesOfYears[year]
    yearFile               = pool.apply_async(grepYear,[files,str(year),archdir,forceOutput,expInfo])
    yearFile, yearMeanFile = getFileNamesForYears(str(year),archdir,expInfo)
    yearFiles.append(yearFile)

  pool.close()
  pool.join()
  add4Cleanup(yearFiles)
  return yearFiles

"""compute yearly mean files """
def computeYearMeanFiles(log,yearmeanfiles,procs):
  years = log['filesOfYears'].keys()
  years.sort()
  log['years'] = years

  pool = multiprocessing.Pool(procs)
  for year in log['years']:
    pool.apply_async(yearmean,[log[year],yearmeanfiles[year],False])
    log['meansOfYears'][year] = yearmeanfiles[year]

  pool.close()
  pool.join()

""" compute timmean + fldmean of masked temperature, salinity and potential density """
def computeMaskedTSRMeans1D(ifiles,varList,initFile,maskFile,exp,archdir,procs):
  pool    = multiprocessing.Pool(procs)
  results = []

  initFile  = cdo.selname(','.join(varList),input = initFile, output = '/'.join([archdir,'.TSR_{0}_init'.format(exp)]))
  # compute the 1d (vertical) year mean variance to the initial state
  for ifile in ifiles:
    _ofile = '/'.join([archdir,'.TSR_1D_{0}'.format(os.path.basename(ifile))])
    ofile  = pool.apply_async(globalTempSalRho1D,[ifile,initFile,maskFile,','.join(varList),_ofile])
    ofile  = _ofile
    results.append(ofile)

  pool.close()
  pool.join()

  merged = '/'.join([archdir,'TSR_1D_{0}_complete.nc'.format(options['EXP'])])
  if os.path.exists(merged):
    os.remove(merged)
  merged = cdo.cat(input = ' '.join(sorted(results)),
                   output =  merged)

  # rename vertical axis if present - ignore exit status (return code differs depending on the ncp version)
  if subprocess.call('ncrename -d depth_2,depth -v depth_2,depth -O {0}'.format(merged),shell=True,env=os.environ):
    print("ERROR: ncrename failed")
    print('ERROR: ncrename -d depth_2,depth -v depth_2,depth_2 -O {0}'.format(merged))
#   exit(1)

  return merged

""" compute timmean + fldmean of masked temperature, salinity and potential density """
def computeMaskedTSRMeans2D(ifiles,varList,exp,archdir,procs,log):
  # alternative: substract 30year mean values instead of averaging variations
  ofile = '/'.join([archdir,'TSR_2D_%s_complete_timmean.nc'%(options['EXP'])])
  cdo.selname(','.join(varList),
              input  = log['last30YearsMeanBias'],
              output = ofile)
  ofileMasked = '/'.join([archdir,'TSR_2D_{0}_complete_timmean_masked.nc'.format(options['EXP'])])
  applyMask(ofile,LOG['mask'],ofileMasked)

  return ofileMasked
# pool    = multiprocessing.Pool(procs)
# results = []

# initFile  = cdo.selname(','.join(varList),input = '-seltimestep,1 %s'%(ifiles[0]),
#                                           output = '/'.join([archdir,'TSR_%s_init'%(exp)]))
#
# # compute the 2d surface temperature + salinity variance to the initial
# # state; relevant are the last 30 years because PHC itselt is a 30-year
# # climatology
# for ifile in ifiles:
#   ofile = pool.apply_async(globalTempSalRho2D,[ifile,initFile,','.join(varList),archdir])
#   ofile = '/'.join([archdir,'TSR_2D_%s'%(os.path.basename(ifile))])
#   results.append(ofile)
#
# pool.close()
# pool.join()
#
# merged = '/'.join([archdir,'TSR_2D_%s_complete.nc'%(options['EXP'])])
# if os.path.exists(merged):
#   os.remove(merged)
# merged = cdo.cat(input = ' '.join(sorted(results)),
#                  output =  merged)
#
# # compute the timmean at the end
# timmeamMerged  = '/'.join([archdir,'TSR_2D_%s_complete_timmean.nc'%(options['EXP'])])
# cdo.timmean(input = merged, output = timmeamMerged)


""" filter for sorting glob results """
def mtime(filename):
  return os.stat(filename).st_mtime

""" split filename into dirname, basenameWithoutExtension, extension """
def splitFileName(filename):
  dirname                     = os.path.dirname(filename)
  filenameWithoutExt, extname = os.path.splitext(os.path.basename(filename))

  return [dirname,filenameWithoutExt,extname]

""" wrapper for single command line execution """
def call(command):
  proc = subprocess.Popen(command,
      shell  = True,
      stderr = subprocess.PIPE,
      stdout = subprocess.PIPE)
  retvals = proc.communicate()
  print(retvals)

""" execute list of system command in parallel """
def executeInParallel(commandList,procs):

  if not commandList:
    return

  pool = multiprocessing.Pool(procs)
  results = []

  for command in commandList:
    pool.apply_async(call,[command])

  pool.close()
  pool.join()

""" create images form 1d txt output """
def plotOnlineDiagnostics(diagnosticFiles,options):
  plots                       = []
  onlineMeanValueFileName     = options['ARCHDIR']+'/TS_last10YearMean_%s.txt'%(options['EXP'])
  onlineMeanValueFile         = open(onlineMeanValueFileName, 'w')
  onlineMeanValueDocumentName = options['PLOTDIR']+'/TS_last10YearMean_%s.pdf'%(options['EXP'])
  # cat all diagnostics together
  diagnosticJoined = options['ARCHDIR']+'/'+'_'.join([options['EXP'],"diagnostics.txt"])
  dbg(diagnosticFiles)

  # read in all values columnwise
  all_diag = []
  with open(diagnosticJoined,'w') as joinedFile:
    for i, diagFile in enumerate(diagnosticFiles):
      f     = open(diagFile,'r')
      input = f.readlines()
      for lineIndex, line in enumerate(input):
        if (0 == lineIndex and 0 != i):
          continue
        newLine          = re.sub('(\s)+', r'\1', line)
        joinedFile.write(newLine)
        newLine          = newLine.rstrip('\r\n')
        input[lineIndex] = newLine
        all_diag.append(newLine.split(' '))
  all_diag = zip(*all_diag)


  header     = open(diagnosticFiles[0],'r').readline().rstrip().split(' ')
  dateIndex  = header.index('date')

  passages = ['gibraltar',
              'denmark_strait',
              'drake_passage',
              'indonesian_throughflow',
              'scotland_iceland',
              'mozambique',
              'framStrait',
              'beringStrait',
              'barentsOpening',
              'agulhas',
              'agulhas_long',
              'agulhas_longer',
             ]
  ice_fields = ['ice_extent_nh','ice_volume_nh']
  columns2Plot = passages + \
                 ice_fields + \
                 ['total_volume',
#                 'total_salinity',
#                 'total_temperature',
                  ]
  # filter columns according to what is in the input
  columns2Plot = [ x for x in columns2Plot if x in header ]
  dbg(columns2Plot)
  columns2PlotIndeces = {}
  for column in columns2Plot:
    columns2PlotIndeces[column] = header.index(column)

  dbg(columns2PlotIndeces)
  dates = all_diag[dateIndex]
  dates = map(dateutil.parser.parse,dates[1:-1])
  lastCompleteYear = dates[-1].year - 1
  dbg("Online diagnostics:")
  dbg("lastCompleteYear:")
  dbg(lastCompleteYear)

  # write mean values to txt file
  # start with head row
  leftColumnWith   = max(map(len,columns2Plot))
  rightColumnWidth = 8
  onlineMeanValueFile.write('   '.join(["Passage/Parameter".ljust(leftColumnWith),''.rjust(rightColumnWidth)])+"\n")
  # go on with data  
# dbg(all_diag)
  for column in columns2Plot:
    dbg(column)
    dbg(columns2PlotIndeces[column])
    data   = all_diag[columns2PlotIndeces[column]]
    values = map(float,numpy.array(data)[1:-1])
    # compute mean value of last 10 years
    last10YearsDates  = []
    last10YearsValues = []
    for i,date in enumerate(dates):
      if (lastCompleteYear - 10 < date.year) and (date.year <= lastCompleteYear):
        #dbg(date)
        last10YearsDates.append(date)
        if column in passages:
          last10YearsValues.append(values[i]/1E+9)
        else:
          last10YearsValues.append(values[i])
    if column in passages:
      last10YearMean = str(round(numpy.array(last10YearsValues).mean(),2))
    else:
      last10YearMean = '%.2E' % round(numpy.array(last10YearsValues).mean(),2)

    dbgOutput = ' : '.join([column.ljust(leftColumnWith),
                           last10YearMean.rjust(8)])
    dbg(dbgOutput)
    onlineMeanValueFile.write(dbgOutput+"\n")


    # complete timeseries
    ofile = "%s/TS_%s_%s_%s.png"%(options['PLOTDIR'],column,options['EXP'],options['TAG'])
    #if ( not os.path.exists(ofile) or options['FORCE']):
    if column in passages:
      unit = "Sv"
    else:
      unit = ""

    pylab.title("%s , %s : %s TF - all %s years | last 10year mean: %s %s"%(options['EXP'],
                                                                            options['TAG'],
                                                                            column.upper(),
                                                                            len(LOG['years']),
                                                                            str(last10YearMean),
                                                                            unit),
                fontsize=8)
    pylab.grid()
    if column in ice_fields:
      pylab.plot_date(last10YearsDates, last10YearsValues, linestyle='-',marker=',',linewidth=1.0)  
    else:
      pylab.plot_date(dates, values, linestyle='-',marker=',',linewidth=1.0)  
    pylab.savefig(ofile)
    pylab.clf() 
    plots.append(ofile)
#   # first 10 years
#   ofile = "%s/%s_first10Years_%s_%s.png"%(options['ARCHDIR'],column,options['EXP'],options['TAG'])
#   if ( not os.path.exists(ofile) or options['FORCE']):
#     pylab.title("%s , %s : %s TF - first 10y"%(options['EXP'],options['TAG'],column.upper()),fontsize=8)
#     pylab.grid()
#     pylab.plot_date(dates[0:120], values[0:120], linestyle='-',marker=',',linewidth=1.0)  
#     pylab.savefig(ofile)
#     pylab.clf() 
#   plots.append(ofile)
#   # first 20 years
#   ofile = "%s/%s_first20Years_%s_%s.png"%(options['ARCHDIR'],column,options['EXP'],options['TAG'])
#   if ( not os.path.exists(ofile) or options['FORCE']):
#     pylab.title("%s , %s : %s TF - first 20y"%(options['EXP'],options['TAG'],column.upper()),fontsize=8)
#     pylab.grid()
#     pylab.plot_date(dates[0:240], values[0:240], linestyle='-',marker=',',linewidth=1.0)  
#     pylab.savefig(ofile)
#   plots.append(ofile)

  dbg(plots)
  onlineMeanValueFile.flush()
  onlineMeanValueFile.close()

  # create pdf version of the data table in plots directory 
  # start with postscript
   # cmd = "cat %s |  convert -pointsize 24  label:@-  %s"%(onlineMeanValueFileName,onlineMeanValueDocumentName)
  cmd = "enscript %s -o - | ps2pdf - %s"%(onlineMeanValueFileName,onlineMeanValueDocumentName)
  dbg(cmd)
  if subprocess.check_call(cmd, shell=True, env=os.environ):
    print("Creation of data table %s failed!"%(onlineMeanValueDocumentName))
  if options['DEBUG']:
    print("Created %s"%(onlineMeanValueDocumentName))

  return onlineMeanValueDocumentName

""" process model output for regio mean values """
def processRegionMean(options,mask3D,init,yearMeanFile,regioCodes,regioDepths,regioVars,regioMaskVar):
  regioMeanData = {}
  regioPool     = multiprocessing.Pool(options['PROCS'])
  regioLock     = multiprocessing.Lock()

  # create the regio mask first
  regioMasks    = {}
  for location, regioCode in regioCodes.iteritems():
    ofile     = '/'.join([options['ARCHDIR'],'_'.join(['regioMask',location])+'.nc'])
    ofileTemp = '/'.join([options['ARCHDIR'],'_'.join(['_regioMask',location])+'.nc'])
    cdo.eqc(regioCode,input = '-selname,%s %s'%(regioMaskVar,init),output = ofileTemp)
    regioMasks[location] = ofileTemp#cdo.div(input = '%s %s'%(ofileTemp,ofileTemp),output = ofile)

  # create the mask from 3d mask wet_c
  regioVertMasks    = {}
  for depth in regioDepths:
    depth = str(depth)
    ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioVertMask',depth+'m'])+'.nc'])
    regioVertMasks[depth] = cdo.sellevel(depth,input = mask3D,output = ofile)

  # compute the regional mean values
  for location, regioCode in regioCodes.iteritems():
    regioMeanData[location] = {}
    for depth in regioDepths:
      regioMeanData[location][str(depth)] = {}

      # create a single mask file out of regional and vertical mask
      myMask = cdo.mul(input=' '.join([regioMasks[location],regioVertMasks[str(depth)]]),
                       output='/'.join([options['ARCHDIR'],'_'.join(['regioMask',location,str(depth)+'m'])+'.nc']))

      for varname in regioVars:
        regioMeanData[location][str(depth)][varname] = {}
        ofile = '/'.join([options['ARCHDIR'],'_'.join(['regioMean',location,varname,str(depth)+'m'])+'.nc'])
        regioPool.apply_async(_computeRegioMean,[depth,varname,yearMeanFile,myMask,ofile,hasNewFiles])
        regioLock.acquire()
        regioMeanData[location][str(depth)][varname] = ofile 
        regioLock.release()
  regioPool.close()
  regioPool.join()

  return regioMeanData

""" create a single document out of all plots """
def createOutputDocument(plotdir,firstPage,docname,doctype,debug):
  images = glob.glob(plotdir+'/*.png')
  dbg(images)

  ofiles = []
  for image in images:
    ofile = '.'.join([os.path.splitext(image)[0],doctype])
    if subprocess.check_call("convert %s %s"%(image,ofile),shell=True,env=os.environ):
      print("ERROR: Convert of %s to %s failed!"%(image,ofile))
    if debug:
      print("Convert %s %s"%(image,ofile))

    ofiles.append(ofile)

  if ( not '' == firstPage ):
    ofiles.insert(0,firstPage)

  document = '/'.join([plotdir,'.'.join([docname,doctype])])

  if ('ps' == doctype):
    if subprocess.check_call("psmerge -o%s %s"%(document,' '.join(ofiles)),shell=True,env=os.environ):
      print("joining postscripts together into %s failed!")
  elif ('pdf' == doctype):
    if subprocess.check_call("pdftk %s cat output %s"%(' '.join(ofiles),document),shell=True,env=os.environ):
      print("joining pdfs together into %s failed!")
  else:
    print("Unsupported document output type: %s!"%(doctype))
    return

""" consistent yearmean filename """
def yearMeanFileName(archdir,exp,firstYear,lastYear):
    return  '/'.join([archdir,'_'.join([exp,'%s-%s'%(firstYear,lastYear),'yearmean.nc'])])

def last30yearFileName(archdir,exp,firstYear,lastYear):
    return  '/'.join([archdir,'_'.join([exp,'{0}-{1}'.format(firstYear,lastYear),'allData.nc'])])

def yearMonMeanFileName(archdir,exp,firstYear,lastYear):
    return  '/'.join([archdir,'_'.join([exp,'%s-%s'%(firstYear,lastYear),'yearmonmean.nc'])])

""" mask by devision with maskfile """
def applyMask(ifile,maskfile,ofile):
  return cdo.div(input = '%s %s'%(ifile,maskfile), output = ofile)

""" X-Section plot """
def plotXSection(plotConfig,options):
  cmd = [options['ICONPLOT']]
  for key,value in plotConfig.iteritems():
      cmd.append('%s=%s'%(key,value))
#       '-iFile=%s'%(t_s_rho_Output_1D),
#       '-varName=%s'%(varname),
#       '-oType='+options['OFORMAT'],
#       '-hov',
#       '-rStrg="-"',
#       '-bStrg="-"',
#       '-lStrg=%s'%(varname),
#       '-withLineLabels',
#       '+withLines',
#       '-minVar=%s'%(t_s_rho_PlotSetup[varname]['minVar']),
#       '-maxVar=%s'%(t_s_rho_PlotSetup[varname]['maxVar']),
#       '-numLevs=%s'%(t_s_rho_PlotSetup[varname]['numLevs']),
#       '-colormap=BlWhRe',
##      '-tStrg="%s"'%(title),
#       '-oFile=%s'%(oFile)]
#print(cmd)
#dbg(' '.join(cmd))
#if subprocess.check_call(' '.join(cmd),shell=True,env=os.environ):
#  print('CMD: %s has failed!')

def plotHorizontal(plotConfig,options,hasNewFiles):
  plots = []
  for varname in plotConfig['varNames']:
    _dirname, _filename, _extname = splitFileName(plotConfig['iFile'])
    if ( varname in plotConfig['availableVars'] ):
      # surface plot
      oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-atSurface_'+'_'.join([_filename,plotConfig['tag']])])
      cmd = [options['ICONPLOT'],
             '-iFile={0}'.format(plotConfig['iFile']),
             '-varName={0}'.format(varname),
             '-oType='+options['OFORMAT'],plotConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="{0}"'.format(plotConfig['title']),
             '-oFile={0}'.format(oFile)]

      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
      elif ('chanel' == options['GRID']):
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')
      else:
        cmd.append('')

      if ('maskFile' in plotConfig):
        cmd.append('-maskFile={0}'.format(plotConfig['maskFile']))

      if ('mapType' in plotConfig):
        cmd.append('-mapType={0}'.format(plotConfig['mapType']))

      if ('colormap' in plotConfig):
        cmd.append('-colormap={0}'.format(plotConfig['colormap']))

      if ('opts' in plotConfig):
        cmd.append(plotConfig['opts'])

      if ('limits' in plotConfig):
        if varname in plotConfig['limits']:
          limits = plotConfig['limits'][varname]
          if 'plotLevs' in limits:
            cmd.append('-plotLevs=%s '%(limits['plotLevs']))
          else:
            cmd.append('-minVar=%s -maxVar=%s %s'%(limits['minVar'],limits['maxVar'],limits['mode']))

      cmd = ' '.join(cmd)
      print('================================================================================================================================')
      dbg(cmd)
      print('================================================================================================================================')
      plots.append(cmd)

      # plot for roughly 30m,100m and 200m depth
      #   skip sea surface height
      if ( 1 == len(cdo.showlevel(input = "-selname,%s -seltimestep,1 %s"%(varname,plotConfig['iFile'])))):
        continue

      oFile = '/'.join([options['PLOTDIR'],varname+'_Horizontal-at100m_'+'_'.join([_filename,plotConfig['tag']])])
      title = '%s: last year mean '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(plotConfig['iFile']),
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],plotConfig['sizeOpt'],
             '-rStrg="-"',
             '-maskName=wet_c',
             '-withLineLabels',
             '-tStrg="%s"'%(plotConfig['title']),
             '-oFile=%s'%(oFile)]

      if ( 'box' == options['GRID'] ):
        cmd.append('-limitMap')
        cmd.append('+mapLine')
        cmd.append('-levIndex=5')
      elif ('global' == options['GRID']):
        cmd.append('-levIndex=8')
      else:
        cmd.append('-mapLLC=-40,-80 -mapURC=30,-30')

      if ('maskFile' in plotConfig):
        cmd.append('-maskFile=%s'%(plotConfig['maskFile']))

      if ('mapType' in plotConfig):
        cmd.append('-mapType=%s'%(plotConfig['mapType']))

      if ('colormap' in plotConfig):
        cmd.append('-colormap=%s'%(plotConfig['colormap']))

      if ('opts' in plotConfig):
        cmd.append(plotConfig['opts'])

      if ('limits' in plotConfig):
        limits = plotConfig['limits'][varname]
        if 'plotLevs' in limits:
          cmd.append('-plotLevs=%s '%(limits['plotLevs']))
        else:
          cmd.append('-minVar=%s -maxVar=%s %s'%(limits['minVar'],limits['maxVar'],limits['mode']))

      cmd = ' '.join(cmd)
      dbg(cmd)
      plots.append(cmd)

    else:
      print("Could not find variable: %s in %s"%(varname,plotConfig['iFile']))

  executeInParallel(plots,options['PROCS'])
  return

def maskCellVariables(ifile,maskfile,ofile,force,log):
  # get variable lists for each vertical grid type
  # this has to be done, because the masking cannot be done for 2d abd 3d in the same file
  #
  groups = {
      'cellVars_2d'  : cdo.showname(input = '-selzaxis,1 -selgrid,1 -seltimestep,1 {0}'.format(ifile))[0].split(' '),
      'cellVars_3d'  : cdo.showname(input = '-selzaxis,2 -selgrid,1 -seltimestep,1 {0}'.format(ifile))[0].split(' '),
      'cellVars_ice' : cdo.showname(input = '-selzaxis,4 -selgrid,1 -seltimestep,1 {0}'.format(ifile))[0].split(' ')
      }

  dbg(groups)

  # apply the mask for each group of variables
  masked_groups = {}
  for group, groupVars in groups.iteritems():
    if group in ['cellVars_2d','cellVars_ice']:
      sellevel = ' -sellevidx,1'
    else:
      sellevel = ''
    masked_groups[group] = cdo.div(input = '-selname,{0} {1} {2} {3}'.format(','.join(groupVars),ifile,sellevel, maskfile), output = '{0}_{1}_masked'.format(ifile,group))
  # merge groups together
  cdo.merge(input = '{0}'.format(' '.join(masked_groups.values())),output = ofile)
  
  # remove temp files
  add4Cleanup(masked_groups.values())
  add4Cleanup([ofile])

  return list(chain(*groups.values()))

""" collect yearmean values and remove single yearmean files """
def collectAllYearMeanData(years,archdir,exp,log):
  # find out which years have been processed
  availableYears, availableYearMeanFiles = [], []
  for year in years:
    yfile, ymeanfile = getFileNamesForYears(year,archdir,exp)
    if os.path.exists(ymeanfile):
      availableYears.append(year)
      availableYearMeanFiles.append(str(ymeanfile))
  #
  # pack the remaining year mean files
  if (len(availableYears) < 1):
    return None

  collectYearMeanDataFile = '{0}/{1}_{2}-{3}_yearmean.nc'.format(archdir,exp,availableYears[0],availableYears[-1])
  cdo.cat(input = ' '.join(availableYearMeanFiles),output = collectYearMeanDataFile)
  #
  # remove singe mean files for alle years which have just been packed and log this
  for _i, _file in enumerate(availableYearMeanFiles):
    os.remove(_file)
    log['packedYears'].append(availableYears[_i])

  #
  # return the name of file
  return collectYearMeanDataFile

""" create yearmean and yearmonmean files for a given list of years """
def createYmonmeanForYears(log, givenPlotYear,options):
  # compute list of relevant years to process
  plotYears = computePlotYears(log['years'],givenPlotYear)
  firstPlotYear, lastPlotYear = plotYears[0], plotYears[-1]
  #
  # separate model output into files for each year for latet cat
  #info = splitFilesIntoYears()
  plotFiles = splitFilesIntoYears(plotYears,
                                  log['filesOfYears'],                                                                                                                                                                              
                                  options['ARCHDIR'],
                                  options['FORCE'],
                                  options['EXP'],
                                  options['PROCS'])
  #
  # cat together
  yearMonMeanFile       = yearMonMeanFileName(options['ARCHDIR'], options['EXP'],firstPlotYear,lastPlotYear)
  last30yearFile        = last30yearFileName(options['ARCHDIR'], options['EXP'], firstPlotYear,lastPlotYear)
  collectedPlotDataFile = cdo.cat(input   = ' '.join(plotFiles),
                                  output  = str(last30yearFile+'_all'),
                                  options = '-f nc4c -z zip_1')
  cellVars              = maskCellVariables(collectedPlotDataFile,
                                            log['mask'],
                                            last30yearFile,
                                            False,LOG)
  #
  # produce masked yearmonmean file
  LOG['last30YearsMonMean']  = cdo.ymonmean(input = last30yearFile,
                                            output = yearMonMeanFile,
                                            options = "-f nc4c -z zip_1")
  #
  # produce masked yearmean
  LOG['last30YearsMean']     = cdo.timmean(input = str(yearMonMeanFile),
                                           output = '%s/%s_%s-%s_ym.nc'%(options['ARCHDIR'],
                                                                                      options['EXP'],
                                                                                      firstPlotYear,
                                                                                      lastPlotYear),
                                           options = "-f nc4c -z zip_1")
  LOG['last30YearsMeanBias'] = cdo.sub(input = ' {0} -selname,{1} {2}'.format(LOG['last30YearsMean'],','.join(cellVars),LOG['init']),
                                       output = '%s/%s_%s-%s_ymBias.nc'%(options['ARCHDIR'],
                                                                                      options['EXP'],
                                                                                      firstPlotYear,
                                                                                      lastPlotYear),
                                       options = "-Q -f nc4c -z zip_1")
  add4Cleanup([last30yearFile])
  add4Cleanup([collectedPlotDataFile])
  return

""" compute the list of years to select for plotting: last 30,20,10,5 years of run or given plotyear """
def computePlotYears(availableYears,givenPlotYear):
  # use integer to avaid string encoding missmatch
  availableYears = map(lambda x : int(x), availableYears)
  if None != givenPlotYear:
    givenPlotYear  = int(givenPlotYear)
  #
  # compute simulation length L
  simulationLength = len(availableYears)
  #
  # plotPeriod is
  plotPeriod = min(30,simulationLength/2)
  #   30 if L >= 100
  if (simulationLength >= 100):
    plotPeriod = min(30,plotPeriod)
  #   20 if 100 < L < 50
  elif (simulationLength < 100):
    plotPeriod = min(20,plotPeriod)
  #   10 if 50  < L < 30
  elif (simulationLength < 50):
    plotPeriod = min(10,plotPeriod)
  #    5 if 30  < L
  elif (simulationLength < 30):
    plotPeriod = min(5,plotPeriod)
  #
  # if plotYear is given
  plotYear = givenPlotYear if givenPlotYear != None else availableYears[-1]

  #   check if it is available
  if plotYear in availableYears:
    plotYearIndex = availableYears.index(plotYear) + 1
  else:
  #   use last simulation year otherwise
    print("PLOTYEAR:{0} is not available - last year {1} will be used".format(plotYear, availableYears[-1]))
    plotYear = availableYears[-1]
    plotYearIndex = simulationLength

  plotYears = availableYears[(plotYearIndex-plotPeriod):plotYearIndex]

  dbg('allYears');dbg(availableYears)
  dbg('{0}:{1}'.format(simulationLength-plotPeriod,plotYearIndex))
  dbg('simlength    :'+str(simulationLength))
  dbg('plotYear     :'+str(plotYear))
  dbg('plotYearIndex:'+str(plotYearIndex))
  dbg(plotYears)

  return plotYears
# }}} --------------------------------------------------------------------------
#=============================================================================== 
# MAIN =========================================================================
# script for archiving (suprise,surpise):
#
# [a] transform output to
#   yealy (model output does not do this, yet)
#   ymean (collected in chunks of restart interval, but complete years only)
#   ymommean (last 30 years by default, shorter in case of shorter simulation)
# [b] compute diagnostics/plots (optionally):
#   psi from 30ym
#   moc from 30ym
#   throughflows from total ym 
#   glob. vertical mixing form total monmean
#   vertical cross-sections: s,t,rhopot
#   horizontal: s,t,rhopot,h,u,v
#   reginal mean for s,t
#   sea ice mean thickness and concentration (march, september)
#   mixed-layer-depth
#
# most steps can be (de)activated by the ACTIONS key in the options dictionary
#=============================================================================== 
# INTERNALS {{{ ================================================================
options         = parseOptions()
dbg(options)
LOGFILE         = 'archive.log'
LOG             = loadLog(options)
dbg(LOG)
LOG['options']  = options
plotCommands    = []
cdo             = Cdo()
cdo.cdfMod      = 'netcdf4'
cdo.debug       = options['DEBUG']
cdo.forceOutput = options['FORCE']
cdo164          = ['/sw/squeeze-x64/cdo-1.6.4/bin/cdo','/sw/aix61/cdo-1.6.7/bin/cdo']
cdo169          = ['/sw/squeeze-x64/cdo-1.6.9/bin/cdo','/sw/aix61/cdo-1.6.9/bin/cdo']
cdoMyVersion    = ['/scratch/mpi/CC/mh0287/users/m300064/builds/bin/cdo']
for cdopath in cdoMyVersion:
  if os.path.exists(cdopath):
    cdo.setCdo(cdopath)
# }}}
# BASIC PLOTSETUP FOR MAIN VARIABLES {{{ ======================================
PlotConfig =  {
  't_acc'      : {'plotLevs' : '-2,-1,-0,1,2,5,10,15,20,25,30'},
  's_acc'      : {'plotLevs' : '20,25,28,30,32,33,34,34.5,35,35.5,36,36.5,37,38,40'},
  'rhopot_acc' : {'plotLevs' : '20,25,28,30,32,34,36,38,40'},
  'u_acc'      : {'plotLevs' : '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  'v_acc'      : {'plotLevs' : '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  'h_acc'      : {'plotLevs' : '-5,-2,-1,-0.5,-0.2,-0.1,0.0,0.1,0.2,0.5,1,2,5'},
  }
PlotConfigBias =  {
  't_acc'      : {'maxVar' : '3.0', 'minVar' : '-3.0' , 'numLevs' : '20'},
  's_acc'      : {'maxVar' : '0.2', 'minVar' : '-0.2' , 'numLevs' : '16'},
  'rhopot_acc' : {'maxVar' : '0.6', 'minVar' : '-0.6' , 'numLevs' : '24'},
  }
# }}} # =======================================================================
# -----------------------------------------------------------------------------
# INPUT HANDLING: {{{
#   LOGGING requested
#   print out used defined log entries
if 'LOG' in options.keys():
  keys = options['LOG'].split(',')
  showLogEntries(loadLog(options),keys)
  doExit()
#   POST-PROCESSINg
# Check for files accessablility
iFiles = glob.glob(options['FILEPATTERN'])
dbg("All available files: "+','.join(iFiles))

# SKIP the last file, if the job is still running
if options['JOBISRUNNING']:
  iFiles.pop()
dbg("Files to be used: "+','.join(iFiles))

if 0 == len(iFiles):
  print(usage())
  print("Could not find any result files!")
  doExit(1)
else:
  print("will use \n" + "\n".join(iFiles))

# is the archive/plots dir accessable
for location in ['ARCHDIR','PLOTDIR']:
  if not os.path.isdir(options[location]):
    os.makedirs(options[location])

# look for new files (process only those, unless the FORCE options is not set)
if LOG.has_key('yearsOfFiles'):
  processedFiles = LOG['yearsOfFiles'].keys()
else:
  processedFiles = []
newFiles = set(iFiles) - set(processedFiles)
hasNewFiles = (len(newFiles) != 0)
if hasNewFiles:
  dbg("New Files found:")
  dbg(newFiles)

# get the data dir from iFiles for later use
LOG['dataDir'] = os.path.abspath(os.path.dirname(iFiles[0]))
# }}}
#------------------------------------------------------------------------------
# DATA SCANNING {{{ ============================================================
if 'preproc' in options['ACTIONS']:
  if options['FORCE']:
    scanFilesForTheirYears(iFiles,options['PROCS'],LOG)
  else:
    scanFilesForTheirYears(newFiles,options['PROCS'],LOG)
  dbg(LOG['yearsOfFiles'])

  # compute all contributing simulation years
  allYears = set(); 
  for yearlist in LOG['yearsOfFiles'].values():
    allYears.update(yearlist)

  allYears = list(allYears)
  allYears.sort()
  # drop the last one if job is running
  if options['JOBISRUNNING']:
    allYears.pop()
  dbg(allYears)

  # collect all files, which contain data for given years
  filesOfYears = {}
  for year in allYears:
    yearFiles = []
    for _file, _years in LOG['yearsOfFiles'].iteritems():
      if year in _years:
        yearFiles.append(_file)
    dbg(yearFiles)
    filesOfYears[year] = yearFiles

  LOG['filesOfYears']  = filesOfYears
  LOG['years']         = allYears

if options['DRYRUN']:
  dbg(LOG)

dumpLog()
dbg(LOG['years'])
# }}} ==========================================================================
# COMPUTE INITIAL VALUES FILES FOR LATER BIASES {{{ ============================
LOG['init'] = cdo.seltimestep(1,input =  iFiles[0], output = '/'.join([options['ARCHDIR'],'%s_init.nc'%(options['EXP'])]))
# }}} ==========================================================================
# COMPUTE CELL MASK FOR LATER APPLICATION {{{ ==================================
LOG['mask'] = cdo.selname('wet_c',input = '-seltimestep,1 %s'%(iFiles[0]), output = '/'.join([options['ARCHDIR'],'%s_mask.nc'%(options['EXP'])]))
# }}} ==========================================================================
# COMPUTE NUMBER OF VERTICAL LEVELS {{{ ========================================
LOG['depths'] = cdo.showlevel(input=LOG['mask'])[0].split()
# }}} ==========================================================================
# COMPUTE YEARMEAN FILE {{{ ====================================================
#yearMeanCollection = collectAllYearMeanData(LOG['years'], options['ARCHDIR'], options['EXP'],LOG)
#dbg('last yearMeanCollection: {0}'.format(yearMeanCollection))
#TODO doExit()
# }}} ==========================================================================
# COMPUTE SINGLE MEAN FILE FROM THE LAST COMPLETE 30 YEARS {{{ =================
createYmonmeanForYears(LOG,options['PLOTYEAR'],options)
dumpLog()
doExit()
# }}} ==========================================================================
# PREPARE INPUT FOR PSI CALC {{{
# collect the last 20 years if there are more than 40 years, last 10 otherwise
if 'procPsi' in options['ACTIONS'] or 'plotPsi' in options['ACTIONS']:
  if len(LOG['years']) > 40:
    nyears4psi = 20
  else:
    nyears4psi = 10
  years4Psi = LOG['years'][-(nyears4psi+2):-1]
  dbg(LOG['years'])
  dbg(years4Psi)
  yearInfo             = '-'.join([years4Psi[0],years4Psi[-1]])
  psiModelVariableName = "u_vint_acc"
  psiGlobalFile        = '/'.join([options['ARCHDIR'],'_'.join([psiModelVariableName,yearInfo])+'.nc'])
  cdo.timmean(input = '-selname,{0} -selyear,{1}/{2} {3}'.format(psiModelVariableName,years4Psi[0],years4Psi[-1],LOG['last30Years']),
              output = psiGlobalFile)
  add4Cleanup([psiGlobalFile])

# }}} --------------------------------------------------------------------------
# PREPARE INPUT FOR PROFILES {{{
# TODO: SWITCHED OFF
#varNames = ['t_acc','s_acc','u_acc','v_acc']
##  collect the data first (first 10 years)
#varFiles = []
#for year in LOG['years'][0:10]:
## varFiles.append(cdo.selname(' '.join(varNames),
##   input=" -remapnn,lon=-30.0_lat=-65.0 %s"%(LOG[year]),
##   output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))
#  varFiles.append(cdo.remapnn('lon=-30.0_lat=-65.0' , 
#                              options = '-O',
#                              input = ' -selname,%s %s'%(','.join(varNames),LOG[year]),
#                              output='/'.join([options['ARCHDIR'],'4profile_'+os.path.basename(LOG[year])])))
#
#varData = cdo.cat(input=' '.join(varFiles),
#                    output="%s/_data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP']))
## add absolute velocity
#velocityData = cdo.expr("'velocity=sqrt(u_acc*u_acc+v_acc*v_acc);'",input=varData,output='/'.join([options['ARCHDIR'],
#                                                                                                   'vel4profile.nc']))
#profileOutputFile = "%s/data_4profile_%s.nc"%(options['ARCHDIR'],options['EXP'])
#if os.path.exists(profileOutputFile):
#  os.remove(profileOutputFile)
#profileData = cdo.merge(input='%s %s'%(varData,velocityData),output=profileOutputFile)
#if subprocess.check_call("ncrename -d depth_2,depth -O %s"%(profileOutputFile),shell=True,env=os.environ):
#  print("ERROR: ncrename failed")
# }}} --------------------------------------------------------------------------
# PREPARE INPUT FOR REGIONAL MEANS from yearMean output {{{
# for global grid only
if 'procRegio' in options['ACTIONS']:
  print(' procRegio mask START ------------------------------------------------------')
  # setup
  regioCodes    = {'ArcticOcean' : 2, 'NorthAtlantic' : 4,'TropicalAtlantic' : 5, 'SouthernOcean' : 6}
  if 20 == len(LOG['depths']):
    regioDepths   = [110,215,895,2200]
  elif 40 == len(LOG['depths']):
    regioDepths   = [100,220,960,2290]

  else:
    regioDepths = [LOG['depths'][0],
                   LOG['depths'][len(LOG['depths'])/4],
                   LOG['depths'][len(LOG['depths'])/2],
                   LOG['depths'][len(LOG['depths'])-1]]
  regioVars     = ['t_acc','s_acc']
  regioMaskVar  = 'regio_c'
  dbg(regioDepths)
  dbg(LOG['depths'])
  dbg(len(LOG['depths']))
  regioMeanData = processRegionMean(options,LOG['mask'],LOG['init'],ymFile,regioCodes,regioDepths,regioVars,regioMaskVar)
  print(' procRegio mask FINISH -----------------------------------------------------')
# }}} --------------------------------------------------------------------------
# PREPARE INPUT FOR MOC PLOT {{{
mocMeanFile = '/'.join([options['ARCHDIR'],'mocMean'])
mocFiles    = sorted(glob.glob(options['MOCPATTERN']),key = mtime)
if 'procMoc' in options['ACTIONS']:
  print(' procMoc START ----------------------------------------')
  # collect all MOC files
  dbg(options['MOCPATTERN'])
  mocFiles.pop(0)
  # skip the latest one if job is still running
  if options['JOBISRUNNING']:
    mocFiles.pop()
  # filter out emtpy files - this could happen when model finished before output timestep
  _mocFiles = []
  for mocfile in mocFiles:
    if (0 < os.stat(mocfile).st_size and 'Z' == mocfile[-1]):
      _mocFiles.append(mocfile)
  mocFiles = _mocFiles
  dbg(mocFiles)
  # default is to take the mean value ofthe at 10 years as input for the plotscript
  # this means 120 months, with monthly output, this is 120 timesteps
  mocNeededNSteps = 120
  mocLog          = {}
  scanFilesForTheirNTime(mocFiles,options['PROCS'],mocLog)
  dbg(mocLog)
  if not mocFiles:
    print('no MOC files for processing')
  else:
    # check for the numbe rof timesteps in the last moc file
    mocLastNtime    = int(mocLog[mocFiles[-1]]) - 1 # avoid the last one, might be corrupted
    if ( os.path.exists(mocMeanFile) ):
      os.remove(mocMeanFile)
    if mocNeededNSteps <= mocLastNtime:
      # take the last 120 values for timmeaninput
      mocMeanFile = cdo.timmean(input = '-seltimestep,{0}/{1} {2}'.format(mocLastNtime-mocNeededNSteps+1,mocLastNtime,mocFiles[-1]),
                                output = mocMeanFile)
    else:
      mocMeanFile = cdo.timmean(input = mocFiles[-1], output = mocMeanFile)
    dbg(mocMeanFile)
  print(' procMoc FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# PREPARE DATA FOR T,S,RHOPOT BIAS TO INIT {{{ ---------------------------------
# target is a year mean file of fldmean data, but mean value computation
# should come at the very end of the processing chain
if 'procTSR' in options['ACTIONS']:
  print(' procTSR START ----------------------------------------')
  t_s_rho_Input_1D = []
  t_s_rho_Input_2D = []
  for year in LOG['years']:
    t_s_rho_Input_1D.append(LOG[year])
  t_s_rho_Output_1D = computeMaskedTSRMeans1D(t_s_rho_Input_1D,['t_acc','s_acc','rhopot_acc'],
                                              LOG['init'],LOG['mask'],options['EXP'],options['ARCHDIR'],options['PROCS'])
  dbg(t_s_rho_Output_1D)
  t_s_rho_Output_2D = applyMask(' -selname,t_acc,s_acc,rhopot_acc {0}'.format(LOG['last30YearsMeanBias']),
                                LOG['mask'],
                                '/'.join([options['ARCHDIR'],'TSR_2D_{0}_complete_timmean_masked.nc'.format(options['EXP'])]))
  dbg(t_s_rho_Output_2D)

  dumpLog()
  print(' procTSR FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# DIAGNOSTICS ==================================================================
# PSI {{{
if 'plotPsi' in options['ACTIONS']:
  print(' plotPsi START ----------------------------------------')
  psiSelectionConfig = {
          'indonesian_throughflow' : { 'lonlatbox' : '90,150,-20,40',},
          'gibraltar'              : { 'lonlatbox' : '-20,10,25,50',},
          'north_atlantic'         : { 'lonlatbox' : '-60,20,50,80',},
          'drake_passage'          : { 'lonlatbox' : '-90,-30,-80,-40',},
          'beringStrait'           : { 'lonlatbox' : '-180,-100,30,80',},
          'agulhas'                : { 'lonlatbox' : '10,50,-55,-15',},
          }
  psiSelectionConfig = {}

  plotFile = options['PLOTDIR']+'/'+"_".join(["psi",yearInfo,options['EXP'],options['TAG']+'.png'])
  title    = "Bar. Streamfunction for %s\n (file: %s)"%(options['EXP'],psiGlobalFile)
  cmd      = '{0} {1} {2}'.format(options['CALCPSI'], psiGlobalFile, " DEBUG=1 WRITEPSI=true AREA={0} TITLE='{1}' PLOT={2}".format(options['GRID'],title,plotFile))
  dbg(cmd)
  plotCommands.append(cmd)
  if subprocess.check_call(cmd,shell=True,env=os.environ):
    print("ERROR: CALCPSI failed")
  # plot special areas
  for area, selection in psiSelectionConfig.iteritems():
    title    = "Selected Stream function for %s (%s)"%(area,options['EXP'])
    plotFile = options['PLOTDIR']+'/'+"_".join(["psi",area,options['EXP'],options['TAG']+'.png'])
    cmd      = '{0} {1} {2}'.format(options['CALCPSI'], psiGlobalFile, " AREA=local TITLE='{0}' PLOT={1} BOX={2} ".format(title,plotFile,selection['lonlatbox']))
    dbg(cmd)
    plotCommands.append(cmd)
    if subprocess.check_call(cmd,shell=True,env=os.environ):
      print("ERROR: CALCPSI failed")
  print(' plotPsi FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# HORIZONTAL PLOTS: t,s,u,v,abs(velocity) {{{
if 'plotHorz' in options['ACTIONS']:
  print(' plotHorz START ----------------------------------------')
  # A) last year mean
  print('-------------------------------------------------------')
  dbg(LOG[LOG['years'][-2]])
  dbg(cdo.showname(input = str(LOG[LOG['years'][-2]]))[0].split(' '))
  print('-------------------------------------------------------')
  horizontalConfig = {
    'varNames'      : ['t_acc','s_acc','h_acc','u_acc','v_acc'],
    'iFile'         : str(LOG['meansOfYears'][LOG['years'][-2]]),
    'availableVars' : cdo.showname(input = str(LOG[LOG['years'][-2]]))[0].split(' '),
    'sizeOpt'       : '-xsize=1200 -ysize=800',
    'title'         : '%s: last year mean '%(options['EXP']),
    'tag'           : 'lastYearMean',
    'limits'        : PlotConfig,
  }
  dbg(horizontalConfig)
  plotHorizontal(horizontalConfig,options,hasNewFiles)
  # B) last 30-year mean
  print('-------------------------------------------------------')
  dbg(LOG['last30YearsMean'])
  print('-------------------------------------------------------')
  horizontalConfig = {
    'varNames'      : ['t_acc','s_acc','h_acc','u_acc','v_acc'],
    'iFile'         : str(LOG['last30YearsMean']),
    'availableVars' : cdo.showname(input = str(LOG['last30YearsMean']))[0].split(' '),
    'sizeOpt'       : '-xsize=1200 -ysize=800',
    'title'         : '%s: last 30-year-mean '%(options['EXP']),
    'tag'           : 'last30YearMean',
    'limits'        : PlotConfig,
  }
  plotHorizontal(horizontalConfig,options,hasNewFiles)
  horizontalConfig = {
    'varNames'      : ['HeatFlux_Total_acc','FrshFlux_VolumeTotal_acc', 'FrshFlux_Precipitation_acc', 'FrshFlux_Evaporation_acc','FrshFlux_Runoff_acc','HeatFlux_ShortWave_acc', 'HeatFlux_LongWave_acc','HeatFlux_Sensible_acc', 'HeatFlux_Latent_acc'],
    'iFile'         : str(LOG['last30YearsMean']),
    'availableVars' : cdo.showname(input = str(LOG['last30YearsMean']))[0].split(' '),
    'sizeOpt'       : '-xsize=1200 -ysize=800',
    'title'         : '%s: last 30-year-mean '%(options['EXP']),
    'tag'           : 'last30YearMean',
    'limits'        : PlotConfig,
  }
  plotHorizontal(horizontalConfig,options,hasNewFiles)
  print(' plotHorz FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# THROUGH FLOWS / ONLINE DIAGNOSTICS {{{
# for global grid only
if ('plotTf' in options['ACTIONS']):
  print(' plotTf START ----------------------------------------')
  if ( 'global' == options['GRID'] ):
    diagnosticFiles = sorted(glob.glob(os.path.sep.join([LOG['dataDir'],"oce_diagnostics-*txt"])),key=mtime)
    if options['JOBISRUNNING']:
      diagnosticFiles.pop()
    diagnosticTable = plotOnlineDiagnostics(diagnosticFiles,options)
  else:
    diagnosticTable = ''
  print(' plotTf FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# ATLANTIC X-Section: t,s,rhopot  {{{ ==========================================
# for global grid only
XSectionPlotConfig = PlotConfigBias
XSectionPlotConfig['s_acc']['minVar'] = '-1.0'
XSectionPlotConfig['s_acc']['maxVar'] = '1.0'
XSectionPlotConfig['t_acc']['minVar'] = '-5.0'
XSectionPlotConfig['t_acc']['maxVar'] = '5.0'
if 'plotX' in options['ACTIONS']:
  print(' plotX START  ---------------------------------------')
  if ( 'global' == options['GRID'] ):
    for varname in ['t_acc','s_acc','rhopot_acc']:
      # plot last yearmean state
      oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_'+'_'.join([options['EXP'],options['TAG']])])
      iFile4XSection = str(LOG['last30YearsMean']) # take last 30 year mean file
      # mask by devision
      #TODO:
      iFile4XSection = cdo.div(input  = '-selname,{0} {1} {2}'.format(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                               output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
      add4Cleanup([iFile4XSection])
      if ( 'rhopot_acc' == varname ):
        # substract 1000
        iFile4XSection  = cdo.subc(1000.0,
                                   input = '-selname,{0} {1}'.format(varname,iFile4XSection),
                                   output = os.path.splitext(iFile4XSection)[0]+'_{0}_subc1000.nc'.format(varname))
        add4Cleanup([iFile4XSection])

      title = '{0}: last 30 year mean '.format(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],
             '-resolution=r360x180',
             '-selPoints=150',
             '-rStrg="-"',
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      cmd = ' '.join(cmd)
      dbg(cmd)
      plotCommands.append(cmd)
      subprocess.check_call(cmd,shell=True,env=os.environ)
  
      # plot bias to initialization
      oFile = '/'.join([options['PLOTDIR'],varname+'_AtlanticProfile_BiasToInit'+'_'.join([options['EXP'],options['TAG']])])
      iFile4XSection = str(LOG['last30YearsMeanBias']) # take last 30 yearmean bias to init
      # mask by devision
      #TODO:
      iFile4XSection = cdo.div(input  = '-selname,{0} {1} {2}'.format(','.join(['t_acc','s_acc','rhopot_acc']),iFile4XSection,LOG['mask']),
                               output =  os.path.splitext(iFile4XSection)[0]+'_masked.nc')
      title = '%s: last 30 year mean bias to init '%(options['EXP'])
      cmd = [options['ICONPLOT'],
             '-iFile=%s'%(iFile4XSection),
             '-secMode=circle -secLC=-45,-70 -secRC=30,80',
             '-varName=%s'%(varname),
             '-oType='+options['OFORMAT'],
             '-resolution=r360x180',
             '-selPoints=150',
             '-rStrg="-"',
             '-minVar=%s'%(XSectionPlotConfig[varname]['minVar']),
             '-maxVar=%s'%(XSectionPlotConfig[varname]['maxVar']),
             '-numLevs=%s'%(XSectionPlotConfig[varname]['numLevs']),
             '-colormap=BlueWhiteOrangeRed',  # put white in the middle
             '-withLineLabels',
             '-tStrg="%s"'%(title),
             '-oFile=%s'%(oFile)]
      cmd = ' '.join(cmd)
      dbg(cmd)
      plotCommands.append(cmd)
      subprocess.check_call(cmd,shell=True,env=os.environ)

  print(' plotX FINISH ---------------------------------------')
# }}} --------------------------------------------------------------------------
# SOUTH OCEAN t,s,y,v profile at 30w, 65s  {{{ =================================
#  create hovmoeller-like plots
#for varname in ['t_acc','s_acc','u_acc','v_acc','velocity']:
#  # run icon_plot.ncl
#  oFile = '/'.join([options['ARCHDIR'],varname+'_profileAt-30W-65S'+'_'.join([options['EXP'],options['TAG']])])
#  iFile4XSection = profileData
#  if ( not os.path.exists(oFile+'.png') or options['FORCE']):
#    title = '%s: 10y prfile at 30W,65S '%(options['EXP'])
#    cmd = [options['ICONPLOT'],
#           '-iFile=%s'%(iFile4XSection),
#           '-hov',
#           '-varName=%s'%(varname),
#           '-oType='+options['OFORMAT'],
#           '-rStrg="-"',
#           '-withLineLabels',
#           '-tStrg="%s"'%(title),
#           '-oFile=%s'%(oFile)]
#    dbg(' '.join(cmd))
#    subprocess.check_call(' '.join(cmd),shell=True,env=os.environ)
# }}} --------------------------------------------------------------------------
# REGIO MEAN PROFILES {{{ ================================
# for global grids only
if ( 'procRegio' in options['ACTIONS'] and 'plotRegio' in options['ACTIONS']):
  print(' procRegio START ------------------------------------------------------')
  regioPlotNames = {'t_acc' : 'Temperature','s_acc' : 'Salinity'}
  for location, regioCode in regioCodes.iteritems():
    for varname in regioVars:
      imageCollection = []
      for depth in regioDepths:
        ifile = regioMeanData[location][str(depth)][varname]

        data  = cdo.readArray(ifile,varname).flatten()

        dates = re.sub('(\s)+', r'\1', cdo.showdate(input = ifile)[0]).split(' ')
        dates = map(dateutil.parser.parse,dates)

        title = '%s at %s: %s [%s,%s]'%(location,str(depth)+'m',varname,options['EXP'],options['TAG'])
        ofile = '/'.join([options['PLOTDIR'],'_'.join(['.regioMean',location,varname,str(depth)+'m',options['EXP'],options['TAG']])+'.png'])

        unit  = cdo.showunit(input = ifile)[0]

        pylab.title(title,fontsize=9)
        pylab.grid()
        pylab.xlabel("Years")
        pylab.ylabel('%s [%s]'%(regioPlotNames[varname],unit))
        pylab.plot_date(dates, data, linestyle='-',marker='.')  
        pylab.savefig(ofile)
        imageCollection.append(ofile)
        pylab.clf()
      collectImageToMapByRows(imageCollection,
                              2,
                              '/'.join([options['PLOTDIR'],'_'.join(['regioMean',location,varname,options['EXP'],options['TAG']+'.png'])])) 
  print(' procRegio FINISH -----------------------------------------------------')
# }}} --------------------------------------------------------------------------
# MOC PLOT {{{
if 'plotMoc' in options['ACTIONS']:
  mocPlotCmd = '%s %s'%(options['MOCPLOTTER'],mocMeanFile)
  mocDateInfo= re.search('(\d{4}-\d{2}-\d{2})T',mocFiles[-1]).group(1)

  mocPlotSetup = {
    'TITLE' : 'ICON: %s : %s %s'%(options['EXP'],mocFiles[-1],'[last 10y mean]'),
    'IFILE' : mocMeanFile,
    'OFILE' : '/'.join([options['PLOTDIR'],'MOC_%s_%s'%(options['EXP'],mocDateInfo)]),
    'OTYPE' : 'png',
  }
  dbg(mocPlotSetup)
  for k,v in mocPlotSetup.iteritems():
    os.environ[k] = v
  plotCommands.append(mocPlotCmd)
  if subprocess.check_call(mocPlotCmd,shell=True,env=os.environ):
    print("ERROR: MOCPLOT failed")
  # environment cleanup
  for k in mocPlotSetup.keys():
    os.environ.pop(k)
# }}} --------------------------------------------------------------------------
# MIXED LAYER PLOT PLOT {{{
# take march and septempber from the last couple of years and just plot is
if ('procMld' in options['ACTIONS']):
  mldLastYear = LOG[LOG['years'][-2]]
  mldData     = cdo.timmax(input = ' -selname,mld,wet_c %s'%(mldLastYear), output = '/'.join([options['ARCHDIR'],'_'.join(['mld',os.path.basename(mldLastYear)])]))
if ('plotMld' in options['ACTIONS']):
  mldPlot = '/'.join([options['PLOTDIR'],'_'.join(['mld',os.path.basename(mldLastYear)])])
  mldCmd  = [options['ICONPLOT'],
      '-iFile=%s'%(mldLastYear),
      '-varName=mld',
      '-isIcon',
      '-colormap=WhiteBlueGreenYellowRed',
      '-maxVar=5000','-minVar=500',
      '-maskName=wet_c',
      '-oType='+options['OFORMAT'],'-tStrg="%s: ~C~ maximum MLD in the last year"'%(options['EXP']),
      '-oFile=%s'%(mldPlot)]
  print(mldCmd)
  mldCmd = ' '.join(mldCmd)
  if subprocess.check_call(mldCmd,shell=True,env=os.environ):
    print('CMD: %s has failed!'%(mldCmd))
  plotCommands.append(mldCmd)

  # nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$fslo -varName=mld -oType=png -colormap=WhiteBlueGreenYellowRed -timeStep=$t -oFile=mld_LIN_$t -maxVar=5000 -minVar=0 -maskName=wet_c
  # nclsh /pool/data/ICON/tools/icon_plot.ncl -iFile=$fslo -varName=mld -oType=png -colormap=WhiteBlueGreenYellowRed -timeStep=$t -oFile=mld_LOG_$t -maxVar=5000 -minVar=0 -selMode=halflog -maskName=wet_c
# }}} --------------------------------------------------------------------------
# SEA ICE PLOT {{{
# take march and septempber from the last 30 year mean
if ('procIce' in options['ACTIONS']):
  iceData = cdo.selname('hi_acc,conc_acc',
                        input = str(LOG['last30YearsMonMean']),
                        output = '/'.join([options['ARCHDIR'],'sea_ice_{0}'.format(os.path.basename(LOG['last30YearsMean']))]))
  meanIceThicknessData = cdo.setunit('m',
                                     input = " -expr,'hiMean=hi_acc*conc_acc' {0}".format(iceData),
                                     output = '/'.join([options['ARCHDIR'],'meanIceThickness_{0}'.format(os.path.basename(iceData))]))
if ('plotIce' in options['ACTIONS']):
  for mapTypeIndex,mapType in enumerate(['NHps','SHps']):
    month = [3,9][mapTypeIndex]
    iFile = cdo.selmon(month,
                       input = iceData,
                       output = '/'.join([options['ARCHDIR'],'month{0}_{1}'.format(month,os.path.basename(iceData))]))
    icePlotConfig = {
      'varNames'      : ['hi_acc','conc_acc'],
      'iFile'         : iFile,
      'availableVars' : cdo.showname(input = iceData)[0].split(' '),
      'sizeOpt'       : '-xsize=1200 -ysize=800',
      'title'         : '%s:~C~ SEA ICE month %s, last-30-yearmonmean'%(options['EXP'],month),
      'tag'           : 'last30YearMeanICE_%s'%(mapType),
      'mapType'       : mapType,
      'colormap'      : 'WhBlGrYeRe',
      'maskFile'      : LOG['mask'],
    }
    plotHorizontal(icePlotConfig,options,hasNewFiles)
    iFile = cdo.selmon(month,
                       input = meanIceThicknessData,
                       output = '/'.join([options['ARCHDIR'],'month{0}_{1}'.format(month,os.path.basename(meanIceThicknessData))]))
    icePlotConfig = {
      'varNames'      : ['hiMean'],
      'iFile'         : iFile,
      'availableVars' : cdo.showname(input = meanIceThicknessData)[0].split(' '),
      'sizeOpt'       : '-xsize=1200 -ysize=800',
      'title'         : '%s:~C~ SEA ICE, month %s, last-30-yearmonmean ~C~ mean sea ice thickness'%(options['EXP'],month),
      'tag'           : 'last30YearMeanICE_%s'%(mapType),
      'mapType'       : mapType,
      'colormap'      : 'WhBlGrYeRe',
      'maskFile'      : LOG['mask'],
    }
    plotHorizontal(icePlotConfig,options,hasNewFiles)

# }}} --------------------------------------------------------------------------
# take march and septempber from the last couple of years and just plot is
# T S RHOPOT BIAS PLOT {{{
if 'plotTSR' in options['ACTIONS']:
  # global mean bias over depth and time
  t_s_rho_PlotSetup = PlotConfigBias
  for varname in t_s_rho_PlotSetup.keys():
    oFile = '/'.join([options['PLOTDIR'],varname+'_biasToInit_inDepth_overTime'+'_'.join([options['EXP'],options['TAG']])])
    title = '%s: %s bias to init '%(options['EXP'],varname)
    cmd = [options['ICONPLOT'],
           '-iFile=%s'%(t_s_rho_Output_1D),
           '-varName=%s'%(varname),
           '-oType='+options['OFORMAT'],
           '-hov',
           '-rStrg="-"',
           '-bStrg="-"',
           '-lStrg=%s'%(varname),
  #        '-withLineLabels',
           '+withLines',
           '-minVar=%s'%(t_s_rho_PlotSetup[varname]['minVar']),
           '-maxVar=%s'%(t_s_rho_PlotSetup[varname]['maxVar']),
           '-numLevs=%s'%(t_s_rho_PlotSetup[varname]['numLevs']),
           '-colormap=BlWhRe',
    #      '-tStrg="%s"'%(title),
           '-oFile=%s'%(oFile)]
    cmd = ' '.join(cmd)
    dbg(cmd)
    plotCommands.append(cmd)
    if subprocess.check_call(cmd,shell=True,env=os.environ):
      print('CMD: %s has failed!')

# }}} --------------------------------------------------------------------------
# FINAL DOCUMENT CREATION {{{ ==================================================
if 'finalDoc' in options['ACTIONS']:
  createOutputDocument(options['PLOTDIR'],diagnosticTable,'_'.join(['ALL',options['EXP'],options['TAG']]),options['DOCTYPE'],options['DEBUG'])
# }}} --------------------------------------------------------------------------
for plotcmd in plotCommands:
  print(plotcmd)

doExit
#
#
# vim:fdm=marker
# vim:sw=2
