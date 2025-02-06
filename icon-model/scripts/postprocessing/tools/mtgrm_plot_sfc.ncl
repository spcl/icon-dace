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
;   ncl -n mtgrm_plot_sfc.ncl 'iFile="path/file.nc"' iStation=1 'varName="SHFL"' 'expnum="${expnum}"'
;
; Required Parameters:
; iFile         : input file (e.g. "METEOGRAM_patch001.nc")
; iStation      : station number (1,2,3,...)
; varName       : name of the variable to (contour) plot
;
; Optional Parameter:
; colormap      : string for predefined colormaps of the plot (e.g. 'colormap="BlAqGrYeOrReVi200"')
; expnum        : Secondary title above left upper corner of plot
; oFile         : plot file without extension (set by the output type: oType)
; oType         : output graphics format (ps, eps, png, default: png)
; tRange        : range of time axis, date0:date1,  with date in format YYYYmmddHH.
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
  nvars = dimsizes(cdf_file->sfcvar_name)
  do ivar=0,(nvars(0)-1)
    iname = charactertostring(cdf_file->sfcvar_name(ivar,:))
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
    t0 = closest_val( date0, dates(:) )
    t1 = closest_val( date1, dates(:) )
  end if
  return( (/ t0, t1 /))
end


begin

;----------------------------------------------------------------------
; argument list
;----------------------------------------------------------------------

  if(.not. isvar("oType")) oType = "png" end if
  if(.not. isvar("oFile")) then
    ext      = get_file_suffix(iFile,0)
    oFile    = ext@fBase
  end if
  if (.not. isvar("colormap")) colormap = "BlueDarkRed18" end if
  if (.not. isvar("tRange") )  tRange   = ""              end if

  iStation     = iStation - 1

  if(.not. isvar("xlabel")) then 
   xlabel = "day"
   mxstride = 1
  end if
  xstride      = 6  ; subdivisions of X-Axis labels (small ticks)
  
  if (xlabel .eq. "month" ) then
   mxstride     = 4  
  end if; minorticks xstride
  lday         =(/1,10,20/) ; for option 

;----------------------------------------------------------------------
; read data
;----------------------------------------------------------------------
  
; open netCDF file
  cdf_file  = addfile(iFile,"r") 
  nstations = dimsizes(cdf_file->station_lon)
; get index corresponding to variable and time range
  date_index = get_date_range( cdf_file, tRange)
  t0 = date_index(0)
  t1 = date_index(1)
  ivar      = get_var_index(cdf_file, varName)
  vals      = cdf_file->sfcvalues(t0:t1, ivar, iStation)   ; (time, var, station)
  time      = cdf_file->time_step(t0:t1)

;----------------------------------------------------------------------
; plot setup
;----------------------------------------------------------------------
  
  print("")
  print("Meteogram  '" + charactertostring(cdf_file->station_name(iStation,:))+"'" )
; print("Outputfile '"+str_concat((/oFile,".",oType/))+"' will be created")   ; in "+systemfunc("dirname "+iFile))

  wks  = gsn_open_wks(oType,oFile)                 ; Opens a ps file
  gsn_define_colormap(wks,colormap)  
  date = convert_date(charactertostring(cdf_file->date(0,:)))
  yyyymmdd = sprinti("%0.2i", date@year)+sprinti("%0.2i", date@month)+sprinti("%0.2i", date@day)

  res                   = True
  res@gsnCenterString   = charactertostring(cdf_file->sfcvar_long_name(ivar,:))
  res@gsnCenterStringOrthogonalPosF = 0.06
  res@gsnLeftString     = expnum
  nl = "~C~"             ;  carriage return to write lat, lon on second line
  deg = "~S~o~N~"        ;  degree sign (see https://www.ncl.ucar.edu/Applications/Scripts/text_7.ncl)
  clat = tostring_with_format( cdf_file->station_lat(iStation), "%.2f")
  clon = tostring_with_format( cdf_file->station_lon(iStation), "%.2f")
  res@gsnRightString    = charactertostring(cdf_file->station_name(iStation,:)) + nl + clat + deg + "N, " + clon + deg + "E"
  res@tiXAxisString     = "Time  (from "+yyyymmdd+")"                                  ; x-axis label
  res@tiYAxisString     = "["+ charactertostring(cdf_file->sfcvar_unit(ivar,:)) + "]"  ; y-axis label
  res@gsnCenterStringFontHeightF = 0.02
  res@gsnLeftStringFontHeightF   = 0.018
  res@gsnRightStringFontHeightF  = 0.018
  res@tiYAxisFontHeightF         = 0.018
  res@tiXAxisFontHeightF         = 0.018
  res@tmXBLabelFontHeightF       = 0.014
  res@vpXF                 = 0.15                  ; viewport left margin (def: 0.2)
  res@vpWidthF             = 0.85                  ; The Width of the plot box
; res@vpHeightF            = 0.58                  ; The height of the plot box
  res@tmXTOn               = False                 ; turn off the top tickmarks
  res@xyLineThicknesses    = 5                     ; increase line thickness
  res@xyLineColor          = "blue"                ; set line color

;----------------------------------------------------------------------
; create values and date labels for time axis explicitly
;----------------------------------------------------------------------

  ntimesteps = dimsizes(cdf_file->time_step(:))
  expl_labels = new(ntimesteps(0), string)
  expl_values = new(ntimesteps(0), float)
 
  date = convert_date(charactertostring(cdf_file->date(0,:)))
  last_date = date
  count = 0  
  do i=0,(ntimesteps(0)-1)
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
  res@tmXMajorGridLineColor= "lightblue"
  res@tmXMajorGridThicknessF = 1.0
  res@tmYMajorGrid         = True
  res@tmYMajorGridLineDashPattern = 2
  res@tmYMajorGridLineColor= "lightblue"
  res@tmYMajorGridThicknessF = 1.0

;----------------------------------------------------------------------
; line plot
;----------------------------------------------------------------------

  plot = gsn_csm_xy(wks,time(t0:t1),vals,res)             ; plot line
  
end
