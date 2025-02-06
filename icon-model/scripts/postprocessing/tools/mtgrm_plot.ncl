; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;----------------------------------------------------------------------
; mtgrm_plot.ncl 
;----------------------------------------------------------------------
; Meteogram plot script for vertical variables.
;
; Basic Usage:
;   ncl -n mtgrm_plot.ncl 'iFile="path/file.nc"' iStation=1 'varName="T"' 'expnum="${expnum}"'
;
; Required Parameters:
; iFile         : input file (e.g. "METEOGRAM_patch001.nc")
; iStation      : station number (1,2,3,...)
; varName       : name of the variable to (contour) plot
;                 Use varName=FF to plot wind speed, sqrt(U*U+V*V)
;
; Optional Parameter:
; colormap      : string for predefined colormaps of the plot (e.g. 'colormap="BlAqGrYeOrReVi200"')
;                 Default: BlueDarkRed18.
;                 "Red2Blue" reverses "BlueDarkRed18"; "WhiteBlue" reverses "BlueWhite".
; expnum        : Secondary title above left upper corner of plot
; levMin        : minimum level index (default: 0). A negativ number means from level 0 to nlevs + lev0
; oFile         : name of output plot file
; oType         : output graphics format (ps, eps, png, default: png)
; tRange        : range of time axis, date0:date1,  with date in format YYYYmmddHH.
; zAxis         : vertical axis. "ml" = model levels, "h" = height, "logz" is log(H-Hsurf).
;                 Default is "ml".
;
; Info:
; * There is a meteogram plot script available publically, which might be used to add features.
;   http://www.ncl.ucar.edu/Applications/meteo.shtml
; * Defaults: Note that the default escape character is colon ":", which is
;   changed to "~" in our .hluresfile.  Also the default font is set to helvetica.
;
;----------------------------------------------------------------------
; 12/2011 : F. Prill,   DWD (florian.prill@dwd.de)
; 12/2011 : M. Koehler, DWD (martin.koehler@dwd.de)
; 01/2019 : H. Frank,   DWD (helmut.frank@dwd.de)
;----------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;----------------------------------------------------------------------
; Function to find variable in output list
;----------------------------------------------------------------------

function get_var_index (cdf_file, zname : string)
local iname, ivar, nvars
begin
  nvars = dimsizes(cdf_file->var_name)
  do ivar=0,(nvars(0)-1)
    iname = charactertostring(cdf_file->var_name(ivar,:))
    if(zname .eq. iname) then
      return(ivar)
    end if
  end do
  return -1
end


;----------------------------------------------------------------------
; Function to get substring from string
;----------------------------------------------------------------------

function substring (instr[1]:string, first[1]:numeric, last[1]:numeric)
local instr, first, last, main, p2
begin
   main = stringtochar (instr)   
   if (last .ge. first) then
      p2 = last
   else
      p2 = dimsizes (main) - 2
   end if
   return (chartostring (main(first:p2)))
end


;----------------------------------------------------------------------
; Function to retrieve date from iso8601 string
;    Input:  string in the format of "yymmddThhmmssZ"
;    Output: year, month, day, hour, minute, second
;----------------------------------------------------------------------

function convert_date (zstr : string)
begin
  date = True
  date@year   = stringtointeger(substring(zstr, 0, 3))
  date@month  = stringtointeger(substring(zstr, 4, 5))
  date@day    = stringtointeger(substring(zstr, 6, 7))
  date@hour   = stringtointeger(substring(zstr, 9,10))
  date@minute = stringtointeger(substring(zstr,11,12))
  date@second = stringtointeger(substring(zstr,13,14))
  return date
end


;---------------------------------------------------------------------------
; Function to find the first and last index of the dates in the range tRange
;    Input:  date0:date1
;            with the format yyyymmddhh
;            Missing date0 sets t0 = 0
;            Missing date1 sets t1 to the last date index 
;    Output: Indices t0, t1
;---------------------------------------------------------------------------

function get_date_range (cdf_file, tRange : string)
local t0, t1, date0, date1, datum, datum, dates, itr
begin
  t0 = 0
  t1 = dimsizes(cdf_file->time_step(:)) - 1
  if ( tRange .ne. "")
;   Convert date from string to numbers
    datum := tostring( cdf_file->date(:,0:16)) 
    datum = str_sub_str( datum, "T", "")
    datum = str_sub_str( datum, "Z", "") 
    dates = stringtodouble( datum)

    itr  = str_match( tRange, ":")
    if ( ismissing(itr) )
      date0 = dates(t0)
      date1 = stringtodouble( tRange) * 1.d4
    else
      tt = str_split( tRange, ":")
      if ( dimsizes(tt) .eq. 1)
         tR1 = stringtochar( tRange)
        if ( tR1(0) .eq. ":")
          date0 = dates(t0)
          date1 = stringtodouble( tt(0)) * 1.d4
        else
          date0 = stringtodouble( tt(0)) * 1.d4
          date1 = dates(t1)
        end if
      else
        date0 = stringtodouble( tt(0)) * 1.d4
        date1 = stringtodouble( tt(1)) * 1.d4
      end if
    end if
;   Find indices to closest dates
    t0 = closest_val( date0, dates(:) )
    t1 = closest_val( date1, dates(:) )
  end if
  return( (/ t0, t1 /))
end


begin

;----------------------------------------------------------------------
; argument list
;----------------------------------------------------------------------

  if(.not. isvar("oType")) oType = "eps" end if
  if(.not. isvar("oFile")) then
    ext      = get_file_suffix(iFile,0)
    oFile    = ext@fBase
  end if
  if (.not. isvar("colormap")) colormap = "BlueDarkRed18" end if
  if (.not. isvar("tRange") )  tRange   = ""              end if
  if( isvar("levMin"))
    lev0 = stringtointeger( levMin)
  else
    lev0 = 0
  end if

  if(.not. isvar("zAxis")) zAxis = "ml" end if
  zaxis = str_lower( zAxis)
  zax   = substring( zaxis,0,1)
  if ( zax .eq. "h")
    zAxis = "heights"
  else if ( zaxis .eq. "logz")
    zAxis = "logz"
  else
    zAxis = "ml"
  end if
  end if

  iStation     = iStation - 1

  if(.not. isvar("xlabel")) then 
   xlabel = "day"
   mxstride = 1
  end if

  xstride      = 6  ; subdivisions of X-Axis labels (small ticks) 
  if (xlabel .eq. "month" ) mxstride     = 4  end if; minorticks xstride
  lday         =(/1,10,20/) ; for option xlabel=month 

  ystride      = 10 ; subdivisions of right Y-Axis labels

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(iFile,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable
  ivar2 = -1
  if ( varName .eq. "FF") then
    ivar      = get_var_index(cdf_file, "U")
    ivar2     = get_var_index(cdf_file, "V")
  else
    ivar      = get_var_index(cdf_file, varName)
  end if
  nlevs     = cdf_file->var_nlevs(ivar)

  if ( lev0 .ge. 0)
    nlev1 = nlevs-1
  else
    nlev1 = nlevs -1 + lev0
    lev0  = 0
  end if

; height - hsurf = 0 are not possible for log. axis
  if ( zaxis .eq. "logz" .and. cdf_file->heights(nlev1,ivar,iStation) - cdf_file->station_hsurf(iStation) .eq. 0)
    nlev1 = nlev1-1
  end if

; Find dates in range
  date_index = get_date_range( cdf_file, tRange)
  t0 = date_index(0)
  t1 = date_index(1)

  if ( ivar2 .ge. 0 )
    vals = transpose(cdf_file->values(t0:t1, lev0:nlev1, ivar, iStation)) ; (time, level, var, station)
    vals2= transpose(cdf_file->values(t0:t1, lev0:nlev1, ivar2,iStation)) ; (time, level, var, station)
    vals = sqrt( vals*vals + vals2*vals2)
  else
    vals = transpose(cdf_file->values(t0:t1, lev0:nlev1, ivar, iStation)) ; (time, level, var, station)
  end if

;----------------------------------------------------------------------
; plot setup
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
; print("Outputfile '"+str_concat((/oFile,".000001.",oType/))+"' will be created in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens an output file
  if ( colormap .eq. "Red2Blue")
     gsn_define_colormap(wks,"BlueDarkRed18")  
     gsn_reverse_colormap(wks)
  else if ( colormap .eq. "BlueWhite")
     gsn_define_colormap(wks,"WhiteBlue")  
     gsn_reverse_colormap(wks)
  else
     gsn_define_colormap(wks,colormap)  
  end if
  end if

  date = convert_date(charactertostring(cdf_file->date(0,:)))
  yyyymmdd = sprinti("%0.2i", date@year)+sprinti("%0.2i", date@month)+sprinti("%0.2i", date@day)

  res                   = True
  if ( ivar2 .ge. 0) then
    res@gsnCenterString   = "wind speed"
  else
    res@gsnCenterString   = charactertostring(cdf_file->var_long_name(ivar,:))
  end if
  res@gsnCenterStringOrthogonalPosF = 0.06
  res@gsnLeftString     = expnum
  nl = "~C~"             ;  carriage return to write lat, lon on second line
  deg = "~S~o~N~"        ;  degree sign (see https://www.ncl.ucar.edu/Applications/Scripts/text_7.ncl)
  clat = tostring_with_format( cdf_file->station_lat(iStation), "%.2f")
  clon = tostring_with_format( cdf_file->station_lon(iStation), "%.2f")
  res@gsnRightString    = charactertostring(cdf_file->station_name(iStation,:)) + nl + clat + deg + "N, " + clon + deg + "E"
  res@tiXAxisString     = "Time  (from "+yyyymmdd+")" ; x-axis label
  res@sfXArray          = cdf_file->time_step(t0:t1)

  if ( zAxis .eq. "ml")
;   model levels as vertical axis
    res@tiYAxisString     = "Model Level"               ; y-axis label
;   res@sfYArray          = cdf_file->var_levels(ivar,lev0:nlev1)
    res@sfYArray          = ispan( lev0, nlev1, 1)
    res@trYReverse        = True                        ; reverses y axis
    res@vpXF                 = 0.11                  ; viewport left margin (def: 0.2)
    res@vpWidthF             = 0.69                  ; viewport width (def: 0.6)
    res@pmLabelBarWidthF     = 0.24                  ; horizontal shift

  else
;   Do not plot a second axis on the right
    isoil = str_match_ind( varName, "_SO")
    if ( ismissing(isoil)) 
;     height as axis
      res@tiYAxisString     = "Height [m]"              ; y-axis label
      res@trYReverse        = False                      ; do not reverses y axis
    else
;     Soil variables, depth as axis
      res@tiYAxisString     = "Depth [m]"              ; y-axis label
      res@trYReverse        = True                      ; do not reverses y axis
    end if

    if ( zaxis .eq. "logz")
      if ( res@tiYAxisString .eq. "Height [m]")
;       log. axis of H-Hsurf
        res@sfYArray          = cdf_file->heights( lev0:nlev1, ivar, iStation) - cdf_file->station_hsurf(iStation)
        res@tiYAxisString     = " Height - Hsurf [m]"
      else
;       log. axis of depth
        res@sfYArray          = cdf_file->heights( lev0:nlev1, ivar, iStation)
      end if
      res@gsnYAxisIrregular2Log    = True
    else
;     linear height axis
      res@sfYArray          = cdf_file->heights( lev0:nlev1, ivar, iStation)
      res@gsnYAxisIrregular2Linear = True
    end if

    res@vpXF                 = 0.13                  ; viewport left margin (def: 0.2)
;   res@vpWidthF             = 0.69                  ; viewport width (def: 0.6)
    res@vpWidthF             = 0.71                  ; viewport width (def: 0.6)
    res@pmLabelBarWidthF     = 0.18                  ; horizontal shift
  end if
  res@cnFillOn          = True                        ; turn on color fill
  res@gsnSpreadColors   = True                        ; use full colormap
  res@gsnCenterStringFontHeightF = 0.02
  res@gsnLeftStringFontHeightF   = 0.018
  res@gsnRightStringFontHeightF  = 0.018
  res@tiYAxisFontHeightF         = 0.018
  res@tiXAxisFontHeightF         = 0.018

;----------------------------------------------------------------------
; contour and label bar resources
;----------------------------------------------------------------------

  res@lbTitleOn            = True                  ; turn on title
  res@lbTitleString        = "["+ charactertostring(cdf_file->var_unit(ivar,:)) + "]"
  res@lbTitleFontHeightF   = .015                  ; label title font size
 ;res@lbTitlePosition      = "Right"               ; label position
  res@lbOrientation        = "vertical"            ; vertical label bar
; res@pmLabelBarWidthF     = 0.24                  ; horizontal shift
  res@lbBoxMinorExtentF    = 0.15                  ; horizontal width
  res@lbLabelFontHeightF   = 0.015                 ; label font size
 ;res@cnLinesOn            = False                 ; optional lines
 ;res@cnFillMode           = "RasterFill"          ; optional raster
 ;res@cnLevelSelectionMode = "ExplicitLevels"
 ;res@gsnMaximize          = True                  ; use full page
; res@vpXF                 = 0.11                  ; viewport left margin (def: 0.2)
; res@vpWidthF             = 0.69                  ; viewport width (def: 0.6)

;----------------------------------------------------------------------
; create values and date labels for time axis explicitly
;----------------------------------------------------------------------

  ntimesteps = dimsizes(cdf_file->time_step(:))
  expl_labels = new(ntimesteps(0), string)
  expl_values = new(ntimesteps(0), float)

  date = convert_date(charactertostring(cdf_file->date(0,:)))
  last_date = date
  count = 0

; do i=0,(ntimesteps(0)-1)
  do i= t0, t1
      date = convert_date(charactertostring(cdf_file->date(i,:)))
      if ( xlabel .eq. "month" ) then 
       if( (i .eq. 0) .or. ((.not. (date@month .eq. last_date@month))\ 
            .or. (.not. ismissing(ind(date@day .eq. lday) )) .and. (date@hour .eq. 0)) ) then
        lcount="true"
       else
        lcount="false"
       end if
      end if
      if (xlabel .eq. "day") then
       if (mod(i,xstride) .eq. 0) then
        lcount="true"
       else
        lcount="false"
       end if
      end if
    if ( lcount .eq. "true" ) then
      date_str = ""
      if ((i .eq. 0) .or. (.not. (date@$xlabel$ .eq. last_date@$xlabel$)) ) then
        date_str = sprinti("%0.2i", date@month) \
        +          "-" + sprinti("%0.2d", date@day) + " "
      end if
      expl_values(count) = cdf_file->time_step(i)
      if ( xlabel .eq. "day" ) then
       expl_labels(count) = date_str + sprinti("%0.2d", date@hour) \
       +                   ":" + sprinti("%0.2d", date@minute)
      end if
      if ( (xlabel .eq. "month") .and. (.not. ismissing(ind(date@day .eq. lday) ))  ) then
       if (date@day .ne. 1)  then 
        date_str =sprinti("%0.2d", date@day) + " "
       end if
       expl_labels(count) = date_str   \
       +                   sprinti("%0.2d", date@hour) \
       +                   ":" + sprinti("%0.2d", date@minute)

      end if
      count = count + 1
      last_date = date
    end if
  end do


  res@tmXBMode             = "Explicit"            ; Define own tick mark labels.
  res@tmXBMinorOn          = True
  res@tmXBMinorValues      = cdf_file->time_step(t0:t1:mxstride)
  res@tmXBMinorThicknessF  = 0.1
  res@tmXBValues           = expl_values(0:(count-1))  ; location of explicit labels
  res@tmXBLabels           = expl_labels(0:(count-1))
  res@tmXBLabelAngleF      = 90.             ; change label angle
  res@tmXBLabelFontHeightF = .01             ; Font size
  res@tmXBLabelJust        = "CenterRight"
  ; remove mirrored tick marks (top, right):
  res@tmXTOn               = False
  res@tmYROn               = False

  res@tmXMajorGrid         = True
  res@tmXMajorGridLineDashPattern = 2
  res@tmXMajorGridLineColor= "white"
  res@tmXMajorGridThicknessF = 1.0
  res@tmYMajorGrid         = True
  res@tmYMajorGridLineDashPattern = 2
  res@tmYMajorGridLineColor= "white"
  res@tmYMajorGridThicknessF = 1.0

;----------------------------------------------------------------------
; If the vertical axis shows model levels then
; create an independent axis with level heights
;----------------------------------------------------------------------

  if ( zAxis .eq. "ml")
    expl_labels2 = new(nlevs, string)
    do i= 0,nlev1
      expl_labels2 = sprintf("%3.0f m", cdf_file->heights( 0:nlev1, ivar, iStation))
    end do
    xy3 = create "xy3" xyPlotClass wks
    "trYReverse"              : "True"
;   "trYMaxF"                 :  cdf_file->var_levels(ivar, nlev1)
;   "trYMinF"                 :  cdf_file->var_levels(ivar, lev0)
    "trYMaxF"                 :  nlev1
    "trYMinF"                 :  lev0
    "tmYROn"                  : "True"
    "tmYLOn"                  : "False"
    "tmYUseLeft"              : "False"
    "tmYRLabelsOn"            : "True"
    "tmYRMajorLengthF"        : .01
    "tmYRMajorOutwardLengthF" : .01
    "tmXBOn"                  : "False"
    "tmXTOn"                  : "False"
    "tmYLOn"                  : "False"
    "tmYROn"                  : "True"
    "tmYRMode"                : "Explicit"
;   "tmYRValues"              : cdf_file->var_levels(ivar,lev0:nlev1:ystride)
    "tmYRValues"              : ispan( lev0, nlev1, ystride)
    "tmYRLabels"              : expl_labels2(lev0:nlev1:ystride)
    "tmYRLabelJust"           : "CenterRight"
    "tiYAxisSide"             : "Right"
    "tmYRLabelFontHeightF"    : 0.01
    "tiXAxisFontHeightF"      : 0.01
    "tiYAxisFontHeightF"      : 0.03
    "tmYRMinorOn"             : "False"
    "tmYLMinorOn"             : "False"
    "tmXBMinorOn"             : "False"
    "tmXTBorderOn"            : "False"
    "tmXBBorderOn"            : "False"
    "tmYLBorderOn"            : "False"
    end create
    draw(xy3)
  end if

;----------------------------------------------------------------------
; contour plot
;----------------------------------------------------------------------

  plot = gsn_csm_contour(wks,vals,res)
  
end
