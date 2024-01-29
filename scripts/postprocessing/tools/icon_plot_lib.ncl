; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;---------------------------------------------------------------
; ICON PLOT LIB
; This library aimes to provide an abstract way to access and plot ICON data
; fields from NetCDF files. It should obey some general requirements:
; * Dont-Repeat-Yourself-Principle: Every method has one responsibility, which no
;   other method has. Or in other words: Don't just CopyAndPaste.
; * Use as much information from the input files as possible. If something is
;   missing, create a better file or let the user input this piece of
;   information. Examples are coordinates of data variables and bounds of them.
;   There are attributes for this.
; * naming conventions:
;   # check* methods are procedure which can exit the whole program
;   # get* methods are functions with return, what they describe with their
;     name. This is similar to NCL itself, e.g. getfilevardims()
;   # set* methods are procedures
;---------------------------------------------------------------
; @authors       : Ralf Mueller (ralf.mueller@dkrz.de)
;                 Stephan Lorenz (stephan.lorenz@zmaw.de)
; VERSION = "0.0.9"
;-------------------------------------------------------------------------------
; some global defaults
  DEFAULTLON   = "clon"
  DEFAULTLAT   = "clat"
  RAD2DEG      = 45./atan(1.)
  ABORTMSG     = "Abort Script."
  NUMLEVS      = 10
  CDO          = "cdo"
;  DEFAULTCOLOR = "BlWhRe"             ; ncolors = 103
  DEFAULTCOLOR = "BlueDarkRed18"      ; ncolors = 256
  DEFAULTCOLOR = "BlueWhiteOrangeRed" ; ncolors = 256
;-------------------------------------------------------------------------------
undef("setCDO")
procedure setCDO(path2Cdo)
begin
  CDO=path2Cdo
end
;-------------------------------------------------------------------------------
; signum function
undef("sign")
function sign(val)
begin
  if (val.eq.0) then
    return(1)
  else
    return(val/abs(val))
  end if
end
;---------------------------------------------------------------
; return ncl version string
undef("getNclVersion")
function getNclVersion()
begin
  version_string = systemfunc("ncl -V")
  version_array  = str_split(version_string, ".")
  ; use integers instead
  return stringtointeger(version_array)
end
;---------------------------------------------------------------
; return true if ncl uses groups for netcdf4 input (release 6.1.1 with netcdf4 input)
undef("versionWithNC4GroupNotationSupport")
function versionWithNC4GroupNotationSupport()
begin
  return (/6,1,1/)
end
undef("versionIsGE")
function versionIsGE(otherVersion)
begin
  return True
  diff = otherVersion .gt. getNclVersion()
  print(""+diff)
  do i=0,dimsizes(diff)
    if (.not. diff(i)) then
      return False
    end if
  end do
  return True
end
;-------------------------------------------------------------------------------
; return a filehandle that can access variables names without leading '/' even
; if the input it nc4 and ncl has version 6.1.1 or above
undef("nc4Aware_addfile")
function nc4Aware_addfile(filename,debug)
begin
  filehandle = addfile( filename, "r" )
  filetype   = get_file_version(filehandle)
  if debug then
    print("[nc4Aware_addfile]Read in file:'" + filename +"'")
    print("[nc4Aware_addfile]    Filetype:'" + filetype +"'")
    print("[nc4Aware_addfile]NCL  version:'" + str_join(getNclVersion(),".") +"'")
    print("[nc4Aware_addfile]version > 611'" + versionIsGE((/6,1,1/)) )
  end if
  if ((get_file_version(filehandle) .eq. "NETCDF4") .and. versionIsGE((/6,1,1/))) then
    print("[nc4Aware_addfile]switch to root group '/'")
   ;filehandle = filehandle=>/
  end if
  return filehandle
end
;-------------------------------------------------------------------------------
; Is the 'time' variable a dimesions in the input dimensions array?
undef("has_time")
function has_time(dims)
begin
  return any(dims.eq."time")
end
;-------------------------------------------------------------------------------
; Has the input variable the dimension 'time'?
undef("var_has_time")
function var_has_time(varname,filehandle)
begin
  dims  = getfilevardims(filehandle,varname)
  if ( any(dims.eq."time") ) then
    return(True)
  else
    return(False)
  end if
end
;-------------------------------------------------------------------------------
; has a file a certain variable
undef("has_var")
function has_var(filehandle,varname)
begin
  vNames = getfilevarnames(filehandle)
  n      = dimsizes (vNames)
  do i=0,n-1              ; loop through each variable
    if (vNames(i).eq.varname) return True end if
  end do
  return False
end
;---------------------------------------------------------------
; get vertical dimension
undef("getVertDim")
function getVertDim(filehandle,var)
begin
  dimNames = getvardims(var)
  matches = (/"lev","depth","height","z"/)
  do i=0,dimsizes(matches)-1
    retval = str_match(dimNames,matches(i))
    if (.not.ismissing(retval)) return retval end if
  end do
  return(str_match(dimNames,matches(0)))
end
;---------------------------------------------------------------
; procedure warper for array_append_record
undef("array_append_record_proc")
procedure array_append_record_proc(array,addition)
begin
  arrayCopy = array_append_record(array,addition,0)

  array     = arrayCopy
end
;-------------------------------------------------------------------------------
; Read a horizontal field of the variable 'varname' from the given filehandle.
; If it is time dependant or has a vertical axis, use the given values for
; timstep and levelindex for selecting the field.
;
; Currently this limited to 3 dimensions, so that 4D tracer variables are ignored
undef("selIconField")
function selIconField(varname,filehandle,timestep,levelindex)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 3) then        ; dims: (time,lev,cell)
        var = filehandle->$varname$(timestep,levelindex,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (time,cell)
        var = filehandle->$varname$(timestep,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    else
; TODO: support for 4D tracer fields
;     if (noOfDims .eq. 3) then        ; dims: (lev,cell,tracerID)
;       var = filehandle->$varname$(levIndex,:,tracer)
;       usableVar = True
;     end if
      if (noOfDims .eq. 2) then        ; dims: (lev,cell)
        var = filehandle->$varname$(levelindex,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (cell)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varname + " from input file")
      exit
    end if
    return var
end
;-------------------------------------------------------------------------------
; get full 3d ICON variable
undef("selIconVar")
function selIconVar(varname,filehandle,timestep)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 3) then        ; dims: (time,lev,cells)
        var = filehandle->$varname$(timestep,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (time,cells)
        var = filehandle->$varname$(timestep,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (cells)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    else
      if (noOfDims .eq. 2) then        ; dims: (lev,cells)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (cells)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varname + " from input file")
      exit
    end if
    return var
end
;-------------------------------------------------------------------------------
; see selIconField(), but for fields on a regular lonlat grid
undef("selRegularField")
function selRegularField(varname,filehandle,timestep,levelindex)
begin
    dims      = getfilevardims(filehandle,varname)
    noOfDims  = dimsizes(dims)
    noOfDims  = noOfDims(0)
    usableVar = False
    if ( has_time(dims) ) then
      if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
        var = filehandle->$varname$(timestep,levelindex,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 3) then        ; dims: (time,x,y)
        var = filehandle->$varname$(timestep,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    else
      if (noOfDims .eq. 3) then        ; dims: (lev,x,y)
        var = filehandle->$varname$(levelindex,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (x)
        var = filehandle->$varname$(:)
        usableVar = True
      end if
    end if
    if (.not. usableVar) then
      print("Warning! Could not read " +  varname + " from input file")
      exit
    end if
    return var
end
;-------------------------------------------------------------------------------
; see selRegularField(), but for 3d regular variable
undef("selRegularVar")
function selRegularVar(varname,filehandle,timestep)
begin
  dims      = getfilevardims(filehandle,varname)
  noOfDims  = dimsizes(dims)
  noOfDims  = noOfDims(0)
  print(""+noOfDims)
  usableVar = False
  if ( has_time(dims) ) then
    if (timestep.eq.-1) then
      if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
        var = filehandle->$varname$
        usableVar = True
      end if
      if (noOfDims .eq. 3) then        ; dims: (time,x,y)
        var = filehandle->$varname$(:,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = False
      end if
    else
      if (noOfDims .eq. 4) then        ; dims: (time,lev,x,y)
        var = filehandle->$varname$(timestep,:,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 3) then        ; dims: (time,x,y)
        var = filehandle->$varname$(timestep,:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 2) then        ; dims: (x,y)
        var = filehandle->$varname$(:,:)
        usableVar = True
      end if
      if (noOfDims .eq. 1) then        ; dims: (time)
        var = filehandle->$varname$(:)
        usableVar = False
      end if
    end if
  else
    if (noOfDims .eq. 3) then        ; dims: (lev,x,y)
      var = filehandle->$varname$(:,:,:)
      usableVar = True
    end if
    if (noOfDims .eq. 2) then        ; dims: (x,y)
      var = filehandle->$varname$(:,:)
      usableVar = True
    end if
    if (noOfDims .eq. 1) then        ; dims: (x)
      var = filehandle->$varname$(:)
      usableVar = True
    end if
  end if
  if (.not. usableVar) then
    print("Warning! Could not read " +  varname + " from lonlat input file")
    exit
  end if
  return var
end
;-------------------------------------------------------------------------------
; abstract Field selection
undef("selField")
function selField(varname,filehandle,timestep,levelindex,horizontalgridtype)
begin
  if ( horizontalgridtype .eq. "unstructured") then
    var = selIconField(varname,filehandle,timestep,levelindex)
  else
    var = selRegularField(varname,filehandle,timestep,levelindex)
  end if
  return var
end
;-------------------------------------------------------------------------------
; scale the given variable with the given factor, if it is NOT 1
undef("scaleVar")
procedure scaleVar(var,scalefactor)
local scalefactor
begin
    if (scalefactor .ne. 1) then
      var = var*scalefactor
    end if
end
;-------------------------------------------------------------------------------
; Return the grid type of a given ICON variable (cell, vertex, or edge).
; This is based on the naming of the coordinates.
undef("getGridType")
function getGridType(variable)
begin
  str   = variable@coordinates
  gtype = stringtocharacter(str)
  if (gtype(0) .eq. "c") then
    str = "cell"
  end if
  if (gtype(0) .eq. "e") then
    str = "edge"
  end if
  if (gtype(0) .eq. "v") then
    str = "vertex"
  end if
  return str
end
;-------------------------------------------------------------------------------
; Return the size of coordinates
undef("getCoordSizes")
function getCoordSizes(variable,filehandle)
begin
  lonlat  = str_split(variable@coordinates," ")
  lon     = lonlat(0)
  lat     = lonlat(1)
  lonSize = dimsizes(filehandle->$lon$)
  latSize = dimsizes(filehandle->$lat$)

  return (/lonSize,latSize/)
end
;-------------------------------------------------------------------------------
; Return the coordinates of a given variable in radians
undef("getCoordinates")
function getCoordinates(variable,filehandle)
begin
  lonlat = str_split(variable@coordinates," ")
  lon    = lonlat(0)
  lat    = lonlat(1)
  x      = filehandle->$lon$
  y      = filehandle->$lat$

  return (/x,y/)
end
;-------------------------------------------------------------------------------
; Return the coordinates of a given variable in degrees
undef("getLonLats")
procedure getLonLats(variable,filehandle,x,y,isicon)
begin
  lonlat = str_split(variable@coordinates," ")
  lon    = lonlat(0)
  lat    = lonlat(1)
  x      = filehandle->$lon$
  y      = filehandle->$lat$
  if (isicon) then
    x      = x * RAD2DEG
    y      = y * RAD2DEG
  end if
end
;-------------------------------------------------------------------------------
; get variable from a certain file
undef("getVarFromFile")
function getVarFromFile(varname,filename)
begin 
  f = addfile(filename+".nc","r")
  var = f->$varname$
  return var
end
;-------------------------------------------------------------------------------
; read a coordinates variable from icon file and return its lonlat representation
undef("getCoordinateFromFile")
function getCoordinateFromFile(varname,filename)
begin
  var = getVarFromFile(varname,filename)
  return var*RAD2DEG
end
;-------------------------------------------------------------------------------
; Create an ascending array of numbers with a 'halflog' interval, i.e. for
; every decade the numbers 1, 2 and 5 are present.
; Examples:
;   (/1,2,5,10,20,50,100,200,500/) or
;   (/-1,-5e-1,-2e-1,-1e-1,-5e-2,-2e-2,-1e-2,0,1e-2,2e-2,5e-2,1e-1,2e-1,5e-1,1/)
; scaleLimit determines the number of decades on the positive or negative axis of the output array
undef("createLevels")
function createLevels(minVar, maxVar, scaleLimit)
begin
  delimiter = "|"
  tics4positive     = (/1, 2, 5/)
  tics4negative     = (/5, 2, 1/)

  signMin  = sign(minVar)
  signMax  = sign(maxVar)
  if (minVar.eq.0) then 
    logfirst = 0.0
    first = 0.0
  else
    logfirst = ceil(log10(abs(minVar)))
    first = signMin*10^logfirst
  end if
  if (maxVar.eq.0) then
    loglast = 0.0
    last  = 0.0
  else
    loglast  = ceil(log10(abs(maxVar)))
    last  = signMax*10^loglast
  end if


  retval = (/first, last, logfirst, loglast, signMin, signMax, logfirst*signMin - loglast*signMax/)

  myScale    = logfirst
  levels     = ""
  ;levels     = new(1,float)
  if (first .eq. 0) then
    levels = str_concat((/levels,delimiter,"0"/))
    ;array_append_record_proc(levels,0.0)
    print(levels)
  end if 
  current    = first
  if (current .lt. 0) then
    tics = tics4negative
    do scaleChange=1,scaleLimit
      do ticIndex=0,2
        current = - tics(ticIndex)*10^(logfirst-scaleChange)
        levels = str_concat((/levels,delimiter,flt2string(doubletofloat(current))/))
        ;array_append_record_proc(levels,doubletofloat(current))
      end do
    end do
  end if
  if (first.lt.0 .and. last.gt.0) then
    levels = str_concat((/levels,delimiter,"0"/))
    ;array_append_record_proc(levels,0.0)
  end if
  if (last .gt. 0) then
    tics = tics4positive
    if (first.lt.0) then
      startScale = min((/logfirst, loglast/)) - scaleLimit
    else
      startScale = min((/logfirst, loglast/))
    end if
    j = 0
    do while (current.lt.(last/10))
      myScale = startScale + j
      do ticIndex=0,2
        current = tics(ticIndex)*10^(myScale)
        levels = str_concat((/levels,delimiter,flt2string(doubletofloat(current))/))
        ;array_append_record_proc(levels,current)
      end do
      j = j + 1
    end do
    levels = str_concat((/levels,delimiter,flt2string(doubletofloat(last))/))
    ;array_append_record_proc(levels,last)
  end if
  return(str_split(levels,"|"))
 ;print(levels)
 ;return(levels)
end
;---------------------------------------------------------------
; Print all variable names of a given file
undef("printVarNames")
procedure printVarNames(filehandle)
begin
  filevarnames = getfilevarnames(filehandle)
  print("# === variable names in file: "+filevarnames)
end
;---------------------------------------------------------------
; Print information about dimensions and attributes of a given varaiable
undef("printVar")
procedure printVar(varname, filehandle)
begin
  dims  = getfilevardims(filehandle,varname)
  sizes = filevardimsizes(filehandle,varname)
  var   = filehandle->$varname$
  ; Print info about variable, its dimensions and attributes
  print(dimsizes(sizes) + " Dimensions:")
  if(.not.any(ismissing(dims))) then
    do j = 0, dimsizes(dims) -1
      print( j + ") " + dims(j) + ": " + sizes(j))
    end do
  end if
  atts = getfilevaratts(filehandle,varname)
  if(.not.any(ismissing(atts))) then
    do k = 0, dimsizes(atts) -1
      print(atts(k) + ": " +filehandle->$varname$@$atts(k)$)
      if (atts(k) .eq. "coordinates") then
        coordinates = var@$atts(k)$
      else
        coordinates = DEFAULTLON + " " + DEFAULTLAT
      end if
    end do
  end if
end
;---------------------------------------------------------------
; Check, if a given variable name is present in a file
undef("checkVarname")
procedure checkVarname(varname, filehandle)
begin
  filevarnames = getfilevarnames(filehandle)
  if (.not. any(filevarnames.eq.varname)) then
    print("Could not find variable "+varname+"!")
    exit
  end if
  return True
end
;---------------------------------------------------------------
; Check, if the dimensions of to variables are equal
undef("checkDimsOfVars")
procedure checkDimsOfVars(varname0,varname1,filehandle)
begin
   print(getfilevardims(filehandle,varname0))
   print(getfilevardims(filehandle,varname1))
  if (str_join(getfilevardims(filehandle,varname0),"") .ne. str_join(getfilevardims(filehandle,varname1),"")) then
    print("Variables '"+varname0+"' and '"+varname1+"' must have the same dimensions!")
    print("ABORT!")
    exit
  end if
end
;---------------------------------------------------------------
; Read a variable from the given or an optional file
undef("getMaskVar")
function getMaskVar(maskvarname,filehandle,hasOtherMaskfile,maskfilename,timestep,levelindex,plotmode,horizontalgridtype)
begin
  if (hasOtherMaskfile) then
    mfilehandle = addfile( maskfilename+".nc","r" )
  else
    mfilehandle = filehandle
  end if
  checkVarname(maskvarname, mfilehandle)
    if (plotmode.eq."scalar" .or. plotmode.eq."overlay") then
  if (horizontalgridtype.eq."unstructured") then
      maskvar = selIconField(maskvarname,mfilehandle,timestep,levelindex)
    else
      maskvar = selRegularField(maskvarname,mfilehandle,timestep,levelindex)
    end if
  else
    if (plotmode.eq."vector") then
      maskvar = selRegularField(maskvarname,mfilehandle,timestep,levelindex)
    else ; section
      maskvar = selRegularVar(maskvarname,mfilehandle,timestep)
    end if
  end if

  return maskvar
end
;---------------------------------------------------------------
; Read the dimensions of a variables from the given ot an optional file
;
; GLOBALS: maskFile
undef("getMaskDim")
function getMaskDim(maskvarname,filehandle,hasOtherMaskfile)
begin
  if (hasOtherMaskfile) then
    filehandle   = addfile( maskFile+".nc","r" )
  end if
  checkVarname(maskvarname, filehandle)
  maskdims = getfilevardims(filehandle,maskvarname)

  return maskdims
end
;---------------------------------------------------------------
; set new filename bases on the input filename and intput type
undef("setNewFilename")
function setNewFilename(filename,tag,itype,atmlev)
begin
  atmtag = ""
  if (itype .eq. "atm") atmtag = "atmlev"+str_capital(atmlev)+"_" end if
  newFilename = tag+"_"+atmtag+systemfunc("basename "+filename)
  return newFilename
end
;---------------------------------------------------------------
; set filename for the automatically remapped file
undef("setRemapFilename")
function setRemapFilename(filename,itype,resolution,atmlev,operator)
begin
  return setNewFilename(filename,operator+"_"+resolution,itype,atmlev)
end
undef("setZonmeanFilename")
function setZonmeanFilename(filename,itype,atmlev)
begin
  return setNewFilename(filename,"zonmean",itype,atmlev)
end
;---------------------------------------------------------------
; Check, if two files have the same number of time steps. Exit otherwise
undef("checkRemappedFile")
function checkRemappedFile(infilename,remapfilename,varname,itype,atmlev,atmplevs,atmhlevs)
begin
  retval = True

  if (.not. isfilepresent(remapfilename)) return False end if

  orgFile             = addfile (infilename+".nc" , "r")
  orgTime             = getfilevardimsizes(orgFile,"time")

  if (.not. has_var(orgFile,varname)) then
    print("Cannot find varname:"+varname)
    exit
  end if

  orgvertdimname      = getVertDim(orgFile,orgFile->$varname$)
  if (.not. ismissing(orgvertdimname) ) then
    orgvertdim          = orgFile->$orgvertdimname$
  end if

  remapFile           = addfile(remapfilename+".nc","r")
  remapTime           = getfilevardimsizes(remapFile,"time")
  ; test if the given variable can be found
  if (.not. has_var(remapFile,varname)) then
    print("Cannot find "+varname+" in remapped files:"+remapfilename)
    retval = False
    return retval
  else
    remappedvertdimname = getVertDim(remapFile,remapFile->$varname$)
    if (.not. ismissing(remappedvertdimname)) then
      remappedvertdim     = remapFile->$remappedvertdimname$
    end if
  end if


  if ( orgTime .eq. remapTime ) then
    print("#=====================================================================================")
    print("Remapped File '"+remapfilename+"' seems to have right time axes.")
  else
    print("#=====================================================================================")
    print("Remapped File '"+remapfilename+"' seems to be wrong: Has different number of timesteps")
    print("Original file has "+orgTime(0)+" timesteps")
    print("Remapped File has "+remapTime(0)+" timesteps")
    retval = False
  end if

  if (itype .eq. "atm") then
    if (atmlev .eq. "p" .and. (str_join(atmplevs,",") .ne. str_join(remappedvertdim,",")) .and. (.not. ismissing(remappedvertdimname)))
      print("Vertical pressure levels differ")
      retval = False
    else
      print("Vertical pressure levels ok.")
    end if
    if (atmlev .eq. "h" .and. (str_join(atmhlevs,",") .ne. str_join(remappedvertdim,",")) .and. (.not. ismissing(remappedvertdimname)))
      print("Vertical height levels differ")
      retval = False
    else
      print("Vertical height levels ok.")
    end if
  end if
  if (itype .eq. "oce") then
    if ((.not. ismissing(orgvertdimname)) .and. (str_join(orgvertdim,",") .ne. str_join(remappedvertdim,",")) .and. (.not. ismissing(remappedvertdimname)))
      print("Vertical depth levels differ")
      retval = False
    else
      print("Vertical depth levels ok")
    end if
  end if

  return retval
end
;---------------------------------------------------------------
; Perform a remapping (wich CDO) to a regular grid with a given resolution and return the
; filename of the remapped file. If the file is already present, the name is returned only
undef("remapForVecPlot")
procedure remapForVecPlot(iFile,remapFilename,resolution,useMask,plotMode,debug,addvars,operator)
begin
  if (plotMode.eq."section") then vecVars=(/varName/) end if
  print("#=====================================================================================")
  print("Looking for remapped file: "+remapFilename)
  print("Use CDO to perform remapping: Create intermediate file: "+remapFilename)
  print("Remove this intermediate file, if it was NOT automatically created from your input file "+iFile+"!")
  variables = str_join(vecVars,",")
  if ( addvars(0) .ne. "")
    variables = str_concat((/variables,",",str_join(addvars,",")/))
    print(str_join(addvars,","))
    print("variables:"+variables)
  end if
  if ( useMask ) then
  variables=variables +","+maskName
  end if
  if (plotMode.eq."overlay") variables = variables+","+varName end if

  cmd = CDO+" -P 8 -"+operator+","+resolution+" -selname,"+variables+" "+iFile+" "+remapFilename
  if debug then
    print(cmd)
  end if
  system(cmd)
end
;---------------------------------------------------------------
; perform zonal mean for hoffmueller plot
undef("zonmean4HoffmuellerPlot")
procedure zonmean4HoffmuellerPlot(filename,varname,zonmeanfilename)
begin
  cmd = CDO+" zonmean -selname,"+varname+" "+filename+" "+zonmeanfilename
  if debug then
    print(cmd)
  end if
  system(cmd)
end
;---------------------------------------------------------------
; Return the bounds of the coordinates of a given variable
undef("getBoundsOfCoordinates")
function getBoundsOfCoordinates(variable,filehandle)
begin
    lonlat = str_split(variable@coordinates," ")
    lon    = lonlat(0)
    lat    = lonlat(1)


    boundslonName = filehandle->$lon$@bounds
    boundslatName = filehandle->$lat$@bounds
    boundslon     = filehandle->$boundslonName$ * RAD2DEG
    boundslat     = filehandle->$boundslatName$ * RAD2DEG

    return (/boundslon, boundslat/)
end
;---------------------------------------------------------------
; Return the bounds of the coordinates of a given variable
undef("getBoundsFromFile")
function getBoundsFromFile(lonname,latname,filename)
begin
  f             = addfile(filename+".nc","r")
  boundslonName = f->$lonname$@bounds
  boundslatName = f->$latname$@bounds
  boundslon     = f->$boundslonName$ * RAD2DEG
  boundslat     = f->$boundslatName$ * RAD2DEG

  return (/boundslon, boundslat/)
end
;---------------------------------------------------------------
; (Re)Set end of bounds. Required by the implementation of the grid plot
undef("setBoundsEnds")
procedure setBoundsEnds(blon,blat,lonmin,lonmax,latmin,latmax)
begin
  latmax=max(blat)
  latmin=min(blat)
  lonmax=max(blon)
  lonmin=min(blon)

  if ( lonmin.lt.0 ) then
  ; longitudes are given in the range [-180,180]
    lonmin = -180.
    lonmax =  180.
  else
  ; longitudes are given in the range [0,360]
    lonmin =    0.
    lonmax =  360.
  end if
end
;---------------------------------------------------------------
; Bounds check for the grid plot
undef("checkLongitude")
procedure checkLongitude(x,boundslon,lonmin,lonmax,debug)
begin
  dlon      = 80.
  ncell_tot = dimsizes(x)  ; total # of cells
  nfix      = 0            ; count the # of problematic cells

  do icell = 0, ncell_tot-1
    if ( any(boundslon(icell,:).le.(lonmin+dlon)) .and. \
         any(boundslon(icell,:).gt.(lonmax-dlon))       ) then

      bl1=boundslon(icell,0)
      bl2=boundslon(icell,1)
      bl3=boundslon(icell,2)

      boundslon(icell,:) = where(boundslon(icell,:).gt.(lonmax-dlon), \;
                           boundslon(icell,:)-360.,             \; where true
                           boundslon(icell,:) )                  ; where false

      bn1=boundslon(icell,0)
      bn2=boundslon(icell,1)
      bn3=boundslon(icell,2)

      if (debug) then
        print("cell="+icell+" lon="+bl1+", "+bl2+", "+bl3+" changed into " \
                           +" lon="+bn1+", "+bn2+", "+bn3)
      end if

      nfix = nfix +1
    end if
  end do
  print("Longitude(s) of "+nfix+" cells corrected to avoid plotting problem.")
end
;---------------------------------------------------------------
; Set the Vertical Selection mode of the given NCL resource
undef("setLevels")
procedure setLevels(selmode,resource,minvar,maxvar,scalelimit,numlevs,debug)
begin
  if (selmode .eq. "manual") then
    if ( .not. isvar("plotLevs") ) then
      resource@cnLevelSelectionMode = "ManualLevels"
      resource@cnMinLevelValF       = minvar
      resource@cnMaxLevelValF       = maxvar
  ;    diffLog10                 = log10(abs(maxVar-minVar))
  ;;   if (diffLog10 .lt. 0) then
  ;      resource@cnLevelSpacingF    = 10^(floor(diffLog10))/numLevs
  ;;   else
  ;;     resource@cnLevelSpacingF    = 10^(floor(diffLog10))/numLevs
  ;;   end if
      diffspacing               = abs(maxvar-minvar)/(int2flt(numlevs))
      resource@cnLevelSpacingF      = diffspacing
    else
      resource@cnLevelSelectionMode = "ExplicitLevels"
      plotLevels = plotLevs
      if (debug) print("plotlevels = "+plotLevels) end if
      resource@cnLevels             = plotLevels
    end if
  end if

  if (selmode .eq. "halflog") then
    resource@cnLevelSelectionMode = "ExplicitLevels"
    plotLevels = createLevels(minvar,maxvar,scalelimit)
    if (debug) print("plotlevels = "+plotLevels) end if
    resource@cnLevels             = plotLevels
  end if
  if (debug) then
    print("Manual Plotlevel setting......")
    print("  selMode : " + selmode)
    print("  minVal  : " + minvar)
    print("  maxVal  : " + maxvar)
    print("  numLevs: " + numlevs)
    print("  scalelimit: " + scalelimit)
    print("  plotLevels: " + resource@cnLevels)
  end if
end
;---------------------------------------------------------------
; set the lines and their labels
undef("setLineAndLineLabels")
procedure setLineAndLineLabels(resource,withLines,withLineLabels)
begin
; lines/labels inside the plot ================================================
  resource@cnLinesOn                = False
  resource@cnLineLabelsOn           = False
  if (withLineLabels ) then
    withLines                       = True
    resource@cnLineLabelsOn         = True
    resource@cnLineLabelFontHeightF = 0.008
    resource@cnLineLabelFontThicknessF = 0.006
  end if
  if (withLines) resource@cnLinesOn = True end if
end
;---------------------------------------------------------------
; Create a default NCL resource
undef("setDefaultResource")
function setDefaultResource(verticallabelbar,withLines,withLineLabels,fillmode,maptype,makeYLinear)
begin
  resource                          = True
  resource@gsnMaximize              = False
  resource@gsnFrame                 = False
  resource@gsnDraw                  = True
  resource@gsnSpreadColors          = True
  resource@cnFillOn                 = True

  resource@mapType                  = maptype

; lines/labels inside the plot ================================================
  setLineAndLineLabels(resource,withLines,withLineLabels)
  ; representation of missing values ===========================================
  delete(resource@cnMissingValFillPattern);     no longer filling missing value areas
  delete(resource@cnMissingValFillColor);       no longer filling missing value areas
  resource@cnMissingValPerimOn          = True
  resource@cnMissingValFillPattern      = -1;                set the missing value fill pattern to 5
  resource@cnMissingValFillScaleF       = 0.9;              increase the density of the fill pattern (default   = 1.0)
  resource@cnMissingValPerimColor       = "black";          change the missing value perimeter to black
  resource@cnMissingValPerimDashPattern = 1;           set the dash pattern of the missing value perimeter to 1
  resource@cnMissingValPerimThicknessF  = 1.0;         increase the thickness of the missing value perimeter 3X

  FontHeight0               = 0.012
  FontHeight1               = 0.015
  FontHeight2               = 0.017
  FontHeight3               = 0.020

  resource@tiXAxisFontHeightF   = FontHeight0
  resource@tiYAxisFontHeightF   = FontHeight0
  resource@tmXBLabelFontHeightF = FontHeight0
  resource@tmYLLabelFontHeightF = FontHeight0
  resource@tmXBLabelJust        = "CenterCenter"
  resource@gsnStringFontHeightF = FontHeight0

  resource@mpFillOn               = True
; label settings ==============================================================
  resource@lbLabelBarOn             = True
  if (verticallabelbar .eq. True) then
    resource@lbOrientation            = "vertical"
  else
    resource@pmLabelBarHeightF        = 0.12
    resource@pmLabelBarWidthF         = 0.80
    resource@pmLabelBarOrthogonalPosF = 0.15
  end if
  resource@lbLabelAutoStride        = True
  resource@lbLabelFontHeightF       = FontHeight1     ;   color bar labels

  resource@cnFillMode           = fillmode
  resource@cnRasterSmoothingOn  = True
  resource@mpGreatCircleLinesOn = True
  resource@stMinArrowSpacingF   = 0.001

  ; linearize y axis if desired ===============================================
  resource@gsnYAxisIrregular2Linear = makeYLinear

  return resource
end
;---------------------------------------------------------------
undef("setDefaultSectionResource")
function setDefaultSectionResource(points,secpoints,seclc,secrc,withLines,withLineLabels,makeYLinear)
begin
  resource                              = True;           plot mods desired
  resource@gsnFrame                     = False;                don't turn page yet
  resource@gsnDraw                      = False;                don't draw yet
  resource@tmXBMode                     = "Explicit";     explicitly label x-axis
  latDiff = stringtodouble(secrc(1)) - stringtodouble(seclc(1))
  tic = latDiff/(stringtodouble(secpoints)-1)
  p = (points * tic) + stringtodouble(seclc(1))
 ; resource@tmXBValues                   = points
  resource@tmXBValues                   = (/points(0),points(secpoints-1)/);  points to label values
  labels = new(dimsizes(p),string)
  labels = p
;print(""+labels)
 ;resource@tmXBLabels                   = labels
  resource@tmXBLabels                   = (/ seclc(1) +"N, "+ seclc(0)+"E" , secrc(1)+"N, "+secrc(0)+"E" /)

  resource@tmXBLabelFontHeightF  = 0.01
  resource@cnFillOn                     = True;          turn on color
  resource@lbLabelAutoStride            = True;          nice label bar label stride
  resource@gsnSpreadColors              = True;          use full range of colormap
  resource@lbOrientation                = "vertical";    vertical label bar
  resource@pmLabelBarOrthogonalPosF     = -0.04;         move label bar closer to plot
  resource@lbTitleString      = "kg/m^3"
  resource@lbTitlePosition              = "Bottom"
;   resource@gsnYAxisIrregular2Log      = True

;  delete(resource@cnMissingValFillPattern);     no longer filling missing value areas
;  delete(resource@cnMissingValFillColor);       no longer filling missing value areas
  resource@cnMissingValFillColor        = "gray30"
  resource@cnMissingValPerimOn          = True

  resource@cnMissingValFillPattern      = -1;                set the missing value fill pattern to 5
  resource@cnMissingValFillScaleF       = 0.9;              increase the density of the fill pattern (default   = 1.0)
;         resource@cnMissingValPerimOn  = True;             already turned on above
  resource@cnMissingValPerimColor       = "black";          change the missing value perimeter to black
  resource@cnMissingValPerimDashPattern = 1;           set the dash pattern of the missing value perimeter to 1
  resource@cnMissingValPerimThicknessF  = 3.0;         increase the thickness of the missing value perimeter 3X
  ;resource@gsnYAxisIrregular2Log        = True
  resource@cnFillMode           = "RasterFill"
  resource@cnRasterSmoothingOn  = False

; lines/labels inside the plot ================================================
  setLineAndLineLabels(resource,withLines,withLineLabels)

  ; linearize y axis if desired ===============================================
  resource@gsnYAxisIrregular2Linear = makeYLinear

  return resource
end
;---------------------------------------------------------------
undef("setDefaultOverlayResource")
procedure setDefaultOverlayResource(resource)
begin
  resource@gsnDraw            = False
  resource@gsnFrame           = False
  resource@vcVectorDrawOrder  = "Postdraw"
  resource@vcFillArrowWidthF  = 12.0
  resource@lbOrientation      = "Vertical"
  resource@lbTitleString      = "kg/m^3"
  resource@lbTitlePosition    = "Left"
  resource@cnInfoLabelOn      = False
  resource@lbTitleFontHeightF = 0.02
  resource@lbLabelFontHeightF = 0.015
  resource@lbLeftMarginF      = 0.01
  resource@lbLabelBarOn       = False;  switch of the vector label bar, because the speed is show as vector lenghts
end
;---------------------------------------------------------------
undef("setDefaultColors")
procedure setDefaultColors(resource,minval,maxval,debug)
begin
; white is color 135
  if ((minval*maxval) .LT. 0.0) then ; data goes from positiv to negative or vs.
    if (abs(maxval) .GT. abs(minval)) then
      resource@gsnSpreadColorStart = 135  - toint(126*(abs(minval)/abs(maxval)))
      resource@gsnSpreadColorEnd   = 255
    else
      resource@gsnSpreadColorStart = 2
      resource@gsnSpreadColorEnd   = 135 + toint(126*(abs(maxval)/abs(minval)))
    end if
    resource@gsnContourZeroLineThicknessF = 2.0
    if (debug) then
      print("===================== COMPUTE the COLORS to make zero WHITE:")
      print("minval: "+minval)
      print("maxval: "+maxval)
      print("100*(toint(abs(minval)/abs(maxval))): "+toint(100*(abs(minval)/abs(maxval))))
      print("resource@gsnSpreadColorStart:" + resource@gsnSpreadColorStart)
      print("resource@gsnSpreadColorEnd  :" + resource@gsnSpreadColorEnd)
    end if
  end if
end
;---------------------------------------------------------------
undef("shift4SecMap")
procedure shift4SecMap(resource)
begin
  resource@vpWidthF  = 0.6;            set width of plot
  resource@vpHeightF = 0.4;            set height of plot

  resource@vpXF      = 0.2
  resource@vpYF      = 0.9
end
;---------------------------------------------------------------
; Compute the left shift for string with font hight 0.01, which should be
; displayed in the lower left corner of the plot with gsn_test_ndc()
undef("getLshiftForNDCString")
function getLshiftForNDCString(mystring)
begin
  lshift = (strlen(mystring)/2 - 1)/100.0
  if (strlen(mystring) .gt. 4)
  lshift = lshift - lshift/5 + 0.1/strlen(mystring)
  else
  lshift = lshift + 0.01
  end if

  return lshift
end
;---------------------------------------------------------------
; Display a small string in the lower left corner of the plot
undef("setBaseString")
procedure setBaseString(workstation,plot,resource,basestring)
begin
  resource@txFontHeightF = 0.01 ; if this is changed, please adjust 
                                ; the implementation of getLshiftForNDCString()
  gsn_text_ndc(workstation,basestring,getLshiftForNDCString(basestring),.01,resource)
end
;---------------------------------------------------------------
; see setBaseString, but use an optional global argument
;
; GLOBALS: bStrg
undef("setAutomaticBaseString")
procedure setAutomaticBaseString(workstation,plotmode)
begin
  if(isvar("bStrg")) then
    BaseString = bStrg
  else
    BaseString = "Prgr icon_plot.ncl: " + systemfunc("date") + ","+getenv("USER")
  end if
  res               = True
  res@txFontHeightF = 0.011
  if (plotmode.eq."section") res@txPosYF = 0.01 end if
  gsn_text_ndc(workstation,BaseString,getLshiftForNDCString(BaseString),.12,res)
end
;---------------------------------------------------------------
; Set the Captions of a plot (left,right,center,title)
undef("setPlotCaptions")
procedure setPlotCaptions(resource,leftstring,rightstring,centerstring,titlestring)
begin
  resource@gsnLeftString   = leftstring
  resource@gsnRightString  = rightstring
  resource@gsnCenterString = centerstring
  resource@tiMainString    = titlestring
end
;---------------------------------------------------------------
; Set the capitons of a plot automatically
;
; GLOBALS: lStrg, rStrg, tStrg
undef("setAutomaticPlotCaptions")
procedure setAutomaticPlotCaptions(resource,plotmode,varname,filehandle,filename,timestep,levelindex,itype,atmlev,k2c)
local varname
begin
  ; titles
  if (plotmode.eq."vector") then
    varname = vecVars(0)
    left_str = "velocity"
  else
    left_str = varname
  end if
  if(.not. isvar("lStrg")) then
    hastime = var_has_time(varname, filehandle)
    if (hastime) then
      timeval = filehandle->time(timestep)
      dayfrag = timeval - floor(timeval)
      day     = doubletointeger(floor(timeval))
      hour    = doubletointeger(floor(86400*dayfrag/3600))
      minute  = doubletointeger(floor((86400*dayfrag - hour*3600)/60))
      second  = doubletointeger(86400*dayfrag - hour*3600 - minute*60)
      timeStr = str_join((/day,hour,minute,second/),"|")
      left_str   = left_str + " (ts:"+timestep+"["+timeStr+"] "
    else
      left_str = "("
    end if
    if (plotmode .ne. "section")
      if ( .not.ismissing(getVertDim(filehandle,filehandle->$varname$)) ) then
        vertdimname = getVertDim(filehandle,filehandle->$varname$)
        vertdim     = filehandle->$vertdimname$
        vertval     = vertdim(levelindex)
        left_str    = left_str + "lev:"+vertval
        if (atmlev .ne. "m")
          left_str = left_str+vertdim@units+")"
        else
          left_str = left_str +")"
        end if
      end if
    end if
  else
   left_str = lStrg
  end if

  if(.not. isvar("rStrg")) then
    rStrg      = (/filename/)
  end if
  rightstring  = str_join(rStrg," ")

  if(.not. isvar("tStrg")) then
    tStrg        = "ICON"
  end if
  titlestring  = str_join(tStrg," ")
  centerstring = ""

  resource@gsnLeftString      = left_str
  resource@gsnRightString     = rightstring
  resource@gsnCenterString    = centerstring
  resource@tiMainString       = titlestring
  v                           = filehandle->$varname$
  if ( isatt(v,"units") ) then
    resource@lbTitleString      = "["+v@units+"]"
  else
    resource@lbTitleString      = ""
  end if
  resource@tiMainFontHeightF  = 0.02
  if (itype .eq. "atm" .and. varname .eq. "T" .and. k2c) resource@lbTitleString = "C" end if

  resource@lbTitlePosition    = "Bottom"
  resource@cnInfoLabelOn      = False
  resource@lbTitleFontHeightF = 0.02
  resource@lbLabelFontHeightF = 0.015
  resource@lbLeftMarginF      = 0.01
  resource@cnLabelMasking   = True
  resource@cnLabelMasking   = True
end
;---------------------------------------------------------------
; vertical axes label for section plot
undef("setSectionVertLabel")
procedure setSectionVertLabel(resource,itype,atmlev)
begin
  resource@lbTitleString     = "[kg/m^3]"
  label = ""
  reversYAxis = True
  if (itype .eq. "atm") then
    label = "Index"
    if (atmlev .eq. "p") label = "Pressure (Pa)" end if
    if (atmlev .eq. "h")
      label = "Height (m)"
      reversYAxis = False
    end if
  else
    label = "depth below sea (m)"
  end if
  resource@tiYAxisString = label
  resource@trYReverse    = reversYAxis
end
;---------------------------------------------------------------
; Determine the type of the map
undef("setMapType")
procedure setMapType(resource,maptype,centerlon,centerlat,satdist)
begin
  if (maptype .eq. "ortho") then
    resource@mpProjection              = "Orthographic"
    resource@mpPerimOn                 = False;     turn off box around plot
    resource@mpGridAndLimbOn           = True;      draw grid lines and limb line
    resource@mpGridLineDashPattern     = 2;         choose pattern used to draw the grid
    resource@mpCenterLonF              = centerlon
    resource@mpCenterLatF              = centerlat
    resource@vcRefLengthF              = 0.01
    resource@vcRefAnnoArrowLineColor   = "black"; change ref vector color
    resource@vcRefAnnoArrowUseVecColor = False;  do not use vec color for ref
    resource@vcRefAnnoOn               = False;        turns off ref vector annotation
    resource@vcGlyphStyle              = "CurlyVector"; turn on curly vectors
    resource@vcMonoLineArrowColor      = True;        multiple colors desired
    resource@vcVectorDrawOrder         = "PostDraw"; draw vectors last


    resource@mpGridAndLimbOn        = True               ; turn on lat/lon grid lines
    resource@mpGridMaskMode         = "MaskNotOcean"     ; don't draw over land or
    ; inland water bodies

;    resource@mpLandFillColor        = "tan"
;    resource@mpOceanFillColor       = "LightBlue"
;    resource@mpInlandWaterFillColor = "LightBlue"
  end if

  if (maptype .eq. "lonlat")
    resource@mpProjection         = "CylindricalEquidistant"
    resource@mpLimitMode          = "LatLon"
  end if

  if (maptype .eq. "NHps")
    resource@mpProjection        = "Stereographic"
    resource@mpRelativeCenterLon = True
    resource@mpCenterLonF        = 0
    resource@mpRelativeCenterLat = True
    resource@mpCenterLatF        = 90.
    resource@mpLimitMode         = "Corners"
  end if

  if (maptype .eq. "SHps")
    resource@mpProjection        = "Stereographic"
    resource@mpRelativeCenterLon = True
    resource@mpCenterLonF        = 0
    resource@mpRelativeCenterLat = True
    resource@mpCenterLatF        = -90.
    resource@mpLimitMode         = "Corners"
  end if
  if (maptype .eq. "sat")
    resource@mpProjection     = "Satellite"
    resource@mpMinLatF        = resource@mpMinLatF + 10.
    resource@mpMaxLatF        = resource@mpMaxLatF - 10.
    resource@mpMinLonF        = resource@mpMinLonF + 10.
    resource@mpMaxLonF        = resource@mpMaxLonF - 10.
    resource@mpCenterLonF     = centerlon
    resource@mpCenterLatF     = centerlat
    resource@mpSatelliteDistF = satdist
  end if
  if (maptype .eq. "lambert")
    resource@mpProjection            = "LambertConformal"
    resource@mpProjection            = "Hammer"
    resource@mpProjection            = "Mollweide"
    resource@mpProjection            = "Robinson"
    resource@mpProjection            = "WinkelTripel"
    resource@mpProjection            = "Gnomonic"
    resource@mpProjection            = "Aitoff"
 ;   resource@gsnMaskLambertConformal = True
  end if
end
;---------------------------------------------------------------
; Determine the corners of the map to plot
; mapllc = map-lower-left-corner
; mapurl = map-upper-right-corner
; format: (/lon,lat/)
undef("selMapCut")
procedure selMapCut(resource,mapllc,mapurc)
begin
  ; bounds of plotting area (not used for Orthographic)
  resource@mpMinLatF    = mapllc(1)
  resource@mpMaxLatF    = mapurc(1)
  resource@mpMinLonF    = mapllc(0)
  resource@mpMaxLonF    = mapurc(0)

  ; center of view (for Orthographic and Sattelite only)
;  resource@mpCenterLonF = (resource@mpMaxLonF+resource@mpMinLonF)/2.0
;  resource@mpCenterLatF = (resource@mpMaxLatF+resource@mpMinLatF)/2.0

  ; bounds of plotting area for stereographic plots
  resource@mpLeftCornerLonF = mapllc(0)
  resource@mpLeftCornerLatF = mapllc(1)
  resource@mpRightCornerLonF = mapurc(0)
  resource@mpRightCornerLatF = mapurc(1)
end
;---------------------------------------------------------------
; Set the visibility of the map
undef("setMapVisibility")
procedure setMapVisibility(resource,mapline)
begin
  if (mapline) then
    resource@mpGeophysicalLineColor = "foreground"
  else
    resource@mpGeophysicalLineColor = "transparent"
  end if
end
;---------------------------------------------------------------
; exit if minvar is larger than maxvar
undef("checkMinMaxVar")
procedure checkMinMaxVar(minvar,maxvar)
begin
  if ( (.not. ismissing(minvar)) .and. (.not. ismissing(maxvar)) ) then
    if ( minvar .gt. maxvar ) then
      print("minVar has to be larger than maxVar")
      exit
    end if
  end if
end
;---------------------------------------------------------------
; set the color of masked locations
undef("setMaskColor")
procedure setMaskColor(wks,resource)
begin
  resource@cnMissingValFillColor       = "gray90"
  resource@cnMissingValPerimOn         = True;          turn on missing value outlines
  resource@cnMissingValPerimThicknessF = 2.5;   set the missing value outline thickness
end
undef("setCoordsForResouce")
procedure setCoordsForResouce(x,y,variable,resource)
begin
  if (dimsizes(dimsizes(x)).ne.1) then
     variable@lon2d = x
     variable@lat2d = y
  end if
     resource@sfXArray = x
     resource@sfYArray = y
end
;---------------------------------------------------------------
; set the coordinates for an ICON data field
undef("setCoordinates")
function setCoordinates(resource,filehandle,variable)
begin
  ret = True
  if (isatt(variable,"coordinates")) then
    lonlat = str_split(variable@coordinates," ")
    lon    = lonlat(0)
    lat    = lonlat(1)
  else
    lon = "lon"
    lat = "lat"
    ret = False
  end if
  x      = filehandle->$lon$
  y      = filehandle->$lat$
  if ("radian" .eq. x@units) x = x * RAD2DEG end if
  if ("radian" .eq. y@units) y = y * RAD2DEG end if

  setCoordsForResouce(x,y,variable,resource)

  if (ret) then
  return (/x,y/)
  else
  return (/1,1/)
  end if
end
undef("setCoordinatesFromExtraFile")
function setCoordinatesFromExtraFile(lon_name,lat_name,filename,variable,resource)
begin
  x = getCoordinateFromFile(lon_name,filename)
  y = getCoordinateFromFile(lat_name,filename)

  setCoordsForResouce(x,y,variable,resource)

  return (/x,y/)
end
;---------------------------------------------------------------
; set the coordinates bounds for an ICON data field
undef("setBounds")
procedure setBounds(resource,filehandle,x,y,debug)
begin
  ; add the bounds to the resource if bounds are present
  if (isatt(x,"bounds")) then
    xbounds = filehandle->$x@bounds$
    if (debug) print(""+x@bounds) end if
    resource@sfXCellBounds = xbounds
  end if
  if (isatt(y,"bounds")) then
    ybounds = filehandle->$y@bounds$
    if (debug) print(""+y@bounds) end if
    resource@sfYCellBounds = ybounds
  end if
end
;---------------------------------------------------------------
; Print out some information about the used map
undef("showMapInfo")
procedure showMapInfo(resource,maptype,mapline)
begin
  print("mapType: "+maptype)
  print("mapLine: "+mapline)
  print("CentLon/Lat="+resource@mpCenterLonF+"; "+resource@mpCenterLatF)
  print("Min/MaxLatF="+resource@mpMinLatF+"; "+resource@mpMaxLatF)
  print("Min/MaxLonF="+resource@mpMinLonF+"; "+resource@mpMaxLonF)
end
;---------------------------------------------------------------
; preprocessing for atmosphere data, i.e. put variables on height
; or pressure levels:
; use ml2pl/ml2hl for vertical interpolation
undef("preProc4Atm")
function preProc4Atm(filename,atmlev,atmplevs,atmhlevs,debug)
begin
  atmfilename = "atmPreProc_"+systemfunc("basename "+filename)
  if (atmlev .eq. "p") then
    operator = "ml2plx"
    levels = str_join(atmplevs,",")
  else
    operator = "ml2hlx"
    levels = str_join(atmhlevs,",")
  end if

  cmd = CDO+" "+ operator +","+levels+" -selgrid,1 "+ filename +" "+ atmfilename
  if debug then
    print(cmd)
  end if
  system(cmd)

  return atmfilename
end
;---------------------------------------------------------------
; determine the leveltype: hybrid or non-hybrid
undef("getVertDimType")
function getVertDimType(dim)
begin
  dtype = str_match(dim@long_name,"hybrid")
  if (ismissing(dtype)) then
    return "non-hybrid"
  else
    return "hybrid"
  end if
end
;---------------------------------------------------------------
; determine the input type: oce OR atm
; variables on hybrid layers are belong to atm input, everything
; else is oce
undef("getVarLevelType")
function getVarLevelType(filename,varname)
begin
  f       = addfile(filename,"r")
  vertdim = getVertDim(f,f->$varname$)
  if (ismissing(vertdim)) then
    delete(f)
    delete(vertdim)
    return "none"
  else
    retval  = getVertDimType(f->$vertdim$)
    delete(f)
    delete(vertdim)

    return retval
  end if
end
;---------------------------------------------------------------
; determine the horizontal grid type (unstructured or lonlat)
undef("getHorizGridType")
function getHorizGridType(filename,varname,isIcon)
begin
  if (isIcon) then
    print("DEBUG:"+isIcon)
    gridtype = "unstructured"
  else
    print("CDO:"+CDO)
    print("Call cdo griddes to determine the horizontal gridtype")
    gridtype = systemfunc(CDO+" -griddes -selname,"+varname+" -seltimestep,1 "+filename+" | grep gridtype | cut -d ' ' -f 4")
  end if
  print("Found horizontal grid of type: "+gridtype)
  return gridtype
end
;---------------------------------------------------------------
; compute the straight line between two corners
undef("getStraightPath")
function getStraightPath(startCorner,endCorner,npoints)
begin
  startLat = tofloat(startCorner(1))
  endLat   = tofloat(endCorner(1))
  startLon = tofloat(startCorner(0))
  endLon   = tofloat(endCorner(0))

  if (startLon .lt. 0.0) then
    startLon = startLon + 360.0
    endLon   = endLon   + 360.0
  end if

  lats = fspan(startLat,endLat,npoints)
  lons = fspan(startLon,endLon,npoints)
  ; don't let the lons be greater than 360, because these values will not be in
  ; the data. Otherwise an interpolation error will occur
  lons = where(lons.gt.360, lons-360.0, lons)

  distRes     = True
  distRes@lat = lats
  distRes@lon = lons

  return distRes
end
;---------------------------------------------------------------
; compute the distance between two points of great circle
undef("getCirclePath")
function getCirclePath(startCorner,endCorner,npoints)
begin
  startLat = startCorner(1)
  endLat   = endCorner(1)

  startLon = startCorner(0)
  endLon   = endCorner(0)

  distRes     = True
  dist        = gc_latlon(startLat,startLon,endLat,endLon,npoints,2)
  distRes@lat = dist@gclat
  distRes@lon = dist@gclon

  return distRes
end
;---------------------------------------------------------------
; compute the path of the cross-section
undef("getSectionPath")
function getSectionPath(seclc,secrc,secpoints,secmode)
begin
  if (secmode .eq. "straight") then
    pathRes = getStraightPath(seclc,secrc,secpoints)
  else
    pathRes = getCirclePath(seclc,secrc,secpoints)
  end if
  return pathRes
end
;---------------------------------------------------------------
; compute the path of the cross-section but suitable for plotting in a separate Map
undef("getSectionPathForMap")
function getSectionPathForMap(seclc,secrc,secpoints,secmode)
begin
  pathRes     = getSectionPath(seclc,secrc,secpoints,secmode)
  pathRes@lon = where(pathRes@lon.lt.pathRes@lon(0),pathRes@lon + 360.0,pathRes@lon)
  return pathRes
end
;---------------------------------------------------------------
; setting for vertical cross sections based on remapped icon data
undef("setSection")
function setSection(seclc,secrc,secpoints,var,secmode)
begin
  pathRes = getSectionPath(seclc,secrc,secpoints,secmode)

  trans   = linint2_points(var&lon,var&lat,var,True,pathRes@lon,pathRes@lat,2)

  return trans
end
;---------------------------------------------------------------
; setting for vertical cross sections based on remapped icon data
undef("setSectionFromHybridPress")
function setSectionFromHybridPress(seclc,secrc,npoints,var,filehandle,timestep)
begin
  leftlat  = seclc(1)
  rightlat = secrc(1)

  leftlon  = seclc(0)
  rightlon = secrc(0)

  ;config
  pressName  = "PS"
  geopotName = "PHIS"
  hybCords   = (/"hyam","hybm","hyai","hybi"/)
  hyamName   = hybCords(0)
  hybmName   = hybCords(1)
  hyaiName   = hybCords(2)
  hybiName   = hybCords(3)

  nlevels = 47
  p00     = 1.e3;   reference pressure in hPa



  phis = filehandle->$pressName$(0,:,:)   ; dims: (time,lon,lat)

  psfc = filehandle->$pressName$(timestep,:,:)  ; surface pressure, dims: (time,lon.lat)
  psfc = psfc*1.e-2  ; convert to hPa
  hyam = filehandle->$hyamName$*1.e-2  ; convert to hPa
  hybm = filehandle->$hybmName$
  hyai = filehandle->$hyaiName$*1.e-2  ; convert to hPa
  hybi = filehandle->$hybiName$
  klevi = dimsizes(hyai)                ;
  klevm = dimsizes(hyam)                ;

  ; compute 3d pressure
  dims = dimsizes(psfc)
  print("dims:"+str_join((/nlevels,  dims(0), dims(1)/),"|"))
  pres3d  = new ((/nlevels,  dims(0), dims(1)/),double)
  pres3di = new ((/nlevels+1,dims(0), dims(1)/),double)
  theta3d = new ((/nlevels,  dims(0), dims(1)/),double)
  do k = 0, nlevels-1
    pres3d(k,:,:) = psfc(:,:)*hybm(k) + hyam(k)
  end do
  do k = 0, nlevels
    pres3di(k,:,:) = psfc(:,:)*hybi(k) + hyai(k)
  end do

  if (False) then  ; potential temperature
    var = var*(p00/pres3d)^(2./7.)
  end if

  ; plot pressure levels
  plevs = (ispan(100,1000,25))*1.0

  var_at_p           = new ((/dimsizes(plevs),dimsizes(var&lat),dimsizes(var&lon)/),double)
  copy_VarMeta( theta3d,  var_at_p)
  var_at_p@lon       = var&lon
  var_at_p@lat       = var&lat
  var_at_p@units     = var@units
  var_at_p@long_name = var@long_name

  extrapolate = True ; switch to do extrapolation below the ground:
  intmethod   = 1 ; 1: method for temperature, -1: method for geopotential, 0: other vars
  tlow =  var(nlevels-1,:,:)  ; temperature at lowest model level
                                       ; (use ground temperature if available)

  var_at_p = vinth2p_ecmwf(var,hyam,hybm,plevs,pres3di(nlevels,:,:)*100.,\
                            intmethod,1.,1,extrapolate,1,tlow,phis)

  dist  = gc_latlon(leftlat,leftlon,rightlat,rightlon,npoints,2)
  trans = linint2_points_Wrap(var_at_p&lon,var_at_p&lat,var_at_p,True,dist@gclon,dist@gclat,2)

  return trans
end
;---------------------------------------------------------------
; setting for vertical cross sections based on pure icon input
undef("setSectionFromHybridPressICON")
function setSectionFromHybridPressICON(secLC,secRC,npoints,var,filehandle,timestep)
begin
  leftlat  = secLC(1)
  rightlat = secRC(1)

  leftlon  = secLC(0)
  rightlon = secRC(0)

  ;config
  pressName  = "PS"
  geopotName = "PHIS"
  hybCords   = (/"hyam","hybm","hyai","hybi"/)
  hyamName   = hybCords(0)
  hybmName   = hybCords(1)
  hyaiName   = hybCords(2)
  hybiName   = hybCords(3)

  nlevels = 47
  plevels = 20
  p00     = 1.e3;   reference pressure in hPa
  rad2deg = 45./atan(1.)   ; radians to degrees

  u3d = var
  phis  = filehandle->$pressName$(0,:);    dims: (time,cell)
  psfc  = filehandle->$pressName$(timestep,:);   surface pressure, dims: (time,lon.lat)
  psfc  = psfc*1.e-2;   convert to hPa
  hyam  = filehandle->$hyamName$*1.e-2;   convert to hPa
  hybm  = filehandle->$hybmName$
  hyai  = filehandle->$hyaiName$*1.e-2;   convert to hPa
  hybi  = filehandle->$hybiName$
  klevi = dimsizes(hyai);
  klevm = dimsizes(hyam);

  ; compute 3d pressure
  dims = dimsizes(psfc)
  print("dims:"+str_join((/nlevels,  dims/),"|"))
  pres3d  = new ((/nlevels,  dims/),double)
  pres3di = new ((/nlevels+1,dims/),double)
  theta3d = new ((/nlevels,  dims/),double)
  do k = 0, nlevels-1
    pres3d(k,:) = psfc(:)*hybm(k) + hyam(k)
  end do
  do k = 0, nlevels
    pres3di(k,:) = psfc(:)*hybi(k) + hyai(k)
  end do

  theta3d = var
  theta3d@units="K"
  theta3d@long_name="Temperature"

;  if (False) then  ; potential temperature
;   var = var*(p00/pres3d)^(2./7.)
; end if

  x           = filehandle->clon *rad2deg;    cell center, lon
  y           = filehandle->clat *rad2deg;    cell center, lat

  ylat        = fspan(-90,90,91)
  ylat@units  = "degrees_north"
  xlon        = fspan(-180,180,181)
  xlon@units  = "degrees_east"
  ylat!0      = "ylat"
  xlon!0      = "xlon"

  u3dll       = triple2grid(x,y,u3d,xlon,ylat,False)
  pres3dll    = triple2grid(x,y,pres3d,xlon,ylat,False)
  pres3dill   = triple2grid(x,y,pres3di,xlon,ylat,False)
  theta3dll   = triple2grid(x,y,theta3d,xlon,ylat,False)
  phisll      = triple2grid(x,y,phis,xlon,ylat,False)

  u3dll!1 = "ylat"
  u3dll!2 = "xlon"
  u3dll&ylat  = ylat
  u3dll&xlon  = xlon

  phisll!0    = "ylat"
  phisll!1    = "xlon"
  phisll&ylat = ylat
  phisll&xlon = xlon

  copy_VarMeta( u3dll,  pres3dill)
  copy_VarMeta( u3dll,  pres3dll)
  copy_VarMeta( u3dll,  theta3dll)

 ; now we try to interpolate from model levels to pressure levels
  plevs = (ispan(100,1000,25))*1.0

  var_at_p           = new ((/dimsizes(plevs),dimsizes(ylat),dimsizes(xlon)/),double)
  copy_VarMeta(u3dll,  var_at_p)
;  var_at_p@ylat       =ylat
;  var_at_p@xlon       =xlon
  var_at_p@units     = theta3d@units
  var_at_p@long_name = theta3d@long_name

  extrapolate = True ; switch to do extrapolation below the ground:
  intmethod   = 1 ; 1: method for temperature, -1: method for geopotential, 0: other vars
  tlow =  theta3dll(nlevels-1,:,:)  ; temperature at lowest model level
                                       ; (use ground temperature if available)

  var_at_p = vinth2p_ecmwf(theta3dll, hyam, hybm, plevs, pres3dill(nlevels,:,:)*100.,\
                            intmethod,1.,1,extrapolate,1,tlow,phisll)
                            print(min(var_at_p))

  dist  = gc_latlon(leftlat,leftlon,rightlat,rightlon,npoints,2)
  trans = linint2_points_Wrap(xlon,ylat,var_at_p,True,dist@gclon,dist@gclat,2)

  return trans
end
;---------------------------------------------------------------
; Create flags array the the grid plot
undef("getFlags")
function getFlags(var,boundslon,boundslat,resource)
begin
  flags = new(dimsizes(var),logical,"No_FillValue")

  if (resource@mapType .EQ. "NHps" .OR. resource@mapType .EQ. "SHps") then
    flags = True
    return flags
  end if
  do i = 0,dimsizes(var) - 1
    flags(i) = where(all(boundslon(i,:) .gt. resource@mpMaxLonF) .or. \
                     all(boundslon(i,:) .lt. resource@mpMinLonF) .or. \
                     all(boundslat(i,:) .gt. resource@mpMaxLatF) .or. \
                     all(boundslat(i,:) .lt. resource@mpMinLatF), \
                     False, True)
  end do
  return flags
end
;---------------------------------------------------------------
; Set the default vector plot resource
undef("setDefaultVectorPlot")
procedure setDefaultVectorPlot(resource, refmagnitude, reflength, style, mindistance)
begin
  resource@vcRefMagnitudeF         = refmagnitude ; make vectors larger
  resource@vcRefLengthF            = reflength    ; reference vector length
  resource@vcGlyphStyle            = style        ; turn on curly vectors
  resource@vcMinDistanceF          = mindistance  ; thin out vectors
end
;---------------------------------------------------------------
; Return the Colors of a given plot
undef("getColorsFromPlot")
function getColorsFromPlot(plot)
begin
  getvalues plot@contour
    "cnFillColors" : colors
  end getvalues
  return colors
end
;---------------------------------------------------------------
; Return the levels of a given plot
undef("getLevelsFromPlot")
function getLevelsFromPlot(plot)
begin
  getvalues plot@contour
    "cnLevels"     : levels
  end getvalues
  return levels
end
;---------------------------------------------------------------
; Print out some information about the grid plot
undef("printGridPlotInfo")
procedure printGridPlotInfo(colors,levels,flags)
begin
  print(colors)
  print(levels)
  print ("Outside the plot area:           " + dimsizes(ind(flags .eq. False)) + \
      " triangles - not plotted")
end
;---------------------------------------------------------------
; Create a resource for the grid plot
undef("performGridPlot")
procedure performGridPlot(var,levels,colors,boundslon,boundslat,wks,plot,flags,debug)
begin
  pres             = True
  pres@gsEdgesOn   = True   ; Turn on edges
  pres@gsFillIndex = 0      ; Solid fill, the default

  ; First draw the triangles associated with the lowest level.
  i = 0
  vlow = ind(var .lt. levels(i))
  ; if no index vlow is found with values less than levels(0), vlow is missing value:
  if (.not. ismissing(vlow(0))) then
    do j = 0, dimsizes(vlow)-1
      pres@gsFillColor = colors(i)   ; first color
      if ( flags(vlow(i))) then
      gsn_polygon(wks,plot,boundslon(vlow(j),:),boundslat(vlow(j),:),pres)
      end if
    end do
  end if
  if (debug) then
    print("i=" + i + ", Values smaller than       "+levels(i)+": " + dimsizes(vlow) + \
          " triangles considered, color= " + colors(i))
  end if
  delete(vlow)

  ; Now draw the triangles between the lowest and highest levels.
  do i = 0, dimsizes(levels) -2
    vind = ind(var .ge. levels(i) .and. var .lt. levels(i+1))
    ; if no index is found, vind is missing value:
    if (.not. ismissing(vind(0))) then
    do j = 0, dimsizes(vind)-1
      pres@gsFillColor = colors(i+1)
      ;if ( flags(vind(i))) then
      gsn_polygon( wks,plot,   boundslon(vind(j),:),boundslat(vind(j),:),pres)
      ;end if
    end do
    end if
    if (debug) then
      print ("i=" + i + ", Values between " + levels(i) + " and " + levels(i+1) +": "\
             +dimsizes(vind) + " triangles considered, color= " + colors(i+1))
    end if
    delete(vind)
  end do

  ; Finally draw the triangles associated with the highest level.
  i = dimsizes(levels) -1
  vhig = ind(var .ge. levels(i) )

  ; if no index vhig is found with values larger levels(max-1), vhig is missing value:
  if (.not. ismissing(vhig(0))) then
    do j = 0, dimsizes(vhig) -1
      pres@gsFillColor = colors(i+1)
      if ( flags(vhig(j))) then
      gsn_polygon(wks,plot,boundslon(vhig(j),:),boundslat(vhig(j),:),pres)
      end if
    end do
  end if
  if (debug) then
    print ("i=" + i + ", Values   larger than      "+levels(i) +": " \
           + dimsizes(vhig) + " triangles considered, color= " + colors(i+1))
  end if
  delete(vhig)
end
;---------------------------------------------------------------
; mark the cell location
undef("plotCellMarkers")
procedure plotCellMarkers(cell_lons,cell_lats,wks,plot,debug)
begin
  pres             = True
  pres@gsMarkerIndex = 1
  pres@gsMarkerSizeF = 0.003
  do i = 0,dimsizes(cell_lons)-1
    gsn_polymarker(wks,plot,cell_lons(i),cell_lats(i),pres)
  end do
end
;---------------------------------------------------------------
; Plot the ICON grid of a variable to the plot.
undef("plotGrid")
procedure plotGrid(workstation,plot,variable,x,bounds,filehandle,resource,debug)
begin
  boundslon = bounds(0,:,:)
  boundslat = bounds(1,:,:)

  latmax=flt2dble(0.0)
  latmin=flt2dble(0.0)
  lonmax=flt2dble(0.0)
  lonmin=flt2dble(0.0)
  setBoundsEnds(boundslon,boundslat,lonmin,lonmax,latmin,latmax)
  checkLongitude(x,boundslon,lonmin,lonmax,debug)

  flags = getFlags(variable,boundslon,boundslat,resource)

  colors = getColorsFromPlot(plot)
  levels = getLevelsFromPlot(plot)

  if (debug) printGridPlotInfo(colors,levels,flags) end if
  performGridPlot(variable,levels,colors,boundslon,boundslat,workstation,plot,flags,debug)
end
;---------------------------------------------------------------
; Plot the ICON grid: New version, req. ncl 6.2
undef("plotGrid_62")
function plotGrid_62(wks,plot,var,boundslon,boundslat)
begin
  pres             = True
  pres@gsEdgesOn   = True   ; Turn on edges
  pres@gsFillIndex = 0      ; Solid fill, the default
  ;---Create color array for triangles
  ntri     = dimsizes(var)                  ;-- Number of triangles
  gscolors = new(ntri,integer)
  gscolors = -1                           ;-- Initialize to transparent

  getvalues plot@contour
   "cnLevels"     : levels
   "cnFillColors" : colors
  end getvalues
  ;---All triangles less than lowest level
  vind = ind(var .lt. levels(0))
  if(.not.all(ismissing(vind))) then
    gscolors(vind) = colors(0)
  end if

  ;---All triangles inbetween levels
  do i = 1, dimsizes(levels) - 1
    vind = ind(var .ge. levels(i-1) .and. var .lt. levels(i))
    if(.not.all(ismissing(vind))) then
      gscolors(vind) = colors(i)
    end if
  end do

  ;---All triangles higher than highest level
  delete(vind)
  vind = ind(var .ge. levels(dimsizes(levels)-1))
  if(.not.all(ismissing(vind))) then
    gscolors(vind) = colors(dimsizes(levels)-1)
  end if

  pres@gsColors = gscolors
  pres@gsSegments = ispan(0,dimsizes(var) * 3,3)

  return pres
  ; }}}
end
;---------------------------------------------------------------
; read in mask variable
undef("readMaskVar")
function readMaskVar()
begin
end
;---------------------------------------------------------------
; plot vector/streamlines to resource
undef("plotVecOrStream")
function plotVecOrStream(useScalar,useStreamlines,uvar,vvar,scalarvar,resource,workstation,showMap)
begin
  ; streamline setup
  resource@stArrowStride    = 4;           arrows start every third
  resource@stArrowLengthF   = 0.0015;      default is dynamic
  resource@stMinDistanceF   = 0.008;       streamline dist
  resource@stLineThicknessF = 1.0;         streamline thickness
 ; resource@stCrossoverCheckCount = 100;
  if (useScalar) then
    if (useStreamlines) then
      if (showMap)
        vc = gsn_csm_streamline_contour_map(workstation,uvar,vvar,scalarvar,resource)
      else
        vc = gsn_streamline_scalar(workstation,uvar,vvar,scalarvar,resource)
      end if
    else
      if (showMap)
        vc = gsn_csm_vector_scalar_map(workstation,uvar,vvar,scalarvar,resource)
      else
        vc = gsn_csm_vector_scalar(workstation,uvar,vvar,scalarvar,resource)
      end if
    end if
  else
    if (showMap)
      if (useStreamlines) then
        vc = gsn_csm_streamline_map(workstation,uvar,vvar,resource)
      else
        vc = gsn_csm_vector_map(workstation,uvar,vvar,resource)
      end if
    else
      if (useStreamlines) then
        vc = gsn_csm_streamline(workstation,uvar,vvar,resource)
      else
        vc = gsn_csm_vector(workstation,uvar,vvar,resource)
      end if
    end if
  end if
  return vc
end
;---------------------------------------------------------------
; limit the map to given coords
undef("limitMap2Coords")
procedure limitMap2Coords(resource)
begin
  lons = resource@sfXArray
  lats = resource@sfYArray
  resource@mpMinLatF         = floor(min(lats))
  resource@mpMaxLatF         = ceil(max(lats))
  resource@mpMinLonF         = floor(min(lons))
  resource@mpMaxLonF         = ceil(max(lons))


;   bounds of plotting area for stereographic plots
  resource@mpLeftCornerLonF  = floor(min(lons))
  resource@mpLeftCornerLatF  = floor(min(lats))
  resource@mpRightCornerLonF = ceil(max(lons))
  resource@mpRightCornerLatF = ceil(max(lats))
end

;
; vim:ft=ncl
