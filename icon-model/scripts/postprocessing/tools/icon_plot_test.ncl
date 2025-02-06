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
; @author: Ralf Mueller (ralf.mueller@dkrz.de)
;---------------------------------------------------------------
; This script makes contour/vector plots of general ICON data files
; For scalar variables the underlying grid can be used instead of automatic
; contour lines. Both modes are capable of masking the input files before
; plotting, see the command line options below for details.

;------- lib location from the command line ------------------------------------------------------------------
loadscript("icon_plot_lib.ncl")
;------------------------------------------------------------------------------
undef("assert")
function assert(tag,givenResult,expectedResult,errorMessage)
begin
  retval = True
  if (any(givenResult .ne. expectedResult)) then
    print("ERROR found: ================================================")
    print(""+tag)
    print(errorMessage)
    print("Expected Result:"+expectedResult)
    print("Given    Result:"+givenResult)
    retval = False
  else
    print("No Error in "+tag)
    retval = True
  end if

  return retval
end
;------------------------------------------------------------------------------
undef("setSimpleResourceAndWks")
function setSimpleWks(ofile)
begin
  otype               = "pdf"
  otype@wkOrientation = "portrait"
  wks                 = gsn_open_wks(otype,ofile)

  return wks
end
;------------------------------------------------------------------------------
; check version
;;ret = assert("NCL version",getNclVersion(),(/6,1,1/),"wrong version found")
;
; plot cells
undef("plot_cells")
procedure plot_cells(cellList,oFile,iFile)
begin
  File                  = addfile( iFile+".nc", "r" )

  resource              = setDefaultResource(True,True,False,"RasterFill")
  resource@tiMainString = oFile
  wks                   = setSimpleWks(oFile)

  varname               = "t"
  var                   = selField(varname,File,0,0,"unstructured")
  coords                = getCoordinates(var,File)
  setCoordinates(resource,coords(0,:),coords(1,:))
  bounds    = getBoundsOfCoordinates(var,File)
  boundslon = bounds(0,:,:)
  boundslat = bounds(1,:,:)
  plot = gsn_csm_contour_map(wks,var,resource)

  do i=0,dimsizes(cellList)-1
    cell = cellList(i)
    gsn_polygon(wks,plot,boundslon(cell,:),boundslat(cell,:),resource)
  end do
end
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
undef("plot_verts")
procedure plot_verts(vertList,oFile,iFile)
begin
  File                  = addfile( iFile+".nc", "r" )
  resource              = setDefaultResource(True,True,False,"RasterFill")
  resource@tiMainString = oFile
  wks                   = setSimpleWks(oFile)

  varname               = "vort"
  var                   = selField(varname,File,0,0,"unstructured")
  coords                = getCoordinates(var,File)
  setCoordinates(resource,coords(0,:),coords(1,:))
  bounds                = getBoundsOfCoordinates(var,File)
  boundslon             = bounds(0,:,:)
  boundslat             = bounds(1,:,:)
  plot                  = gsn_csm_contour_map(wks,var,resource)

  do i=0,dimsizes(vertList)-1
    vert = vertList(i)-1
    gsn_polygon(wks,plot,boundslon(vert,:),boundslat(vert,:),resource)
  end do
end
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
undef("plot_edges")
procedure plot_edges(edgeList,oFile,iFile)
begin
  File                  = addfile( iFile+".nc", "r" )
  resource              = setDefaultResource(True,True,False,"RasterFill")
  resource@tiMainString = oFile
  wks                   = setSimpleWks(oFile)

  varname               = "vn"
  var                   = selField(varname,File,0,0,"unstructured")
  coords                = getCoordinates(var,File)
  setCoordinates(resource,coords(0,:),coords(1,:))
  bounds                = getBoundsOfCoordinates(var,File)
  boundslon             = bounds(0,:,:)
  boundslat             = bounds(1,:,:)
  plot                  = gsn_csm_contour_map(wks,var,resource)

  do i=0,dimsizes(edgeList)-1
    edge = edgeList(i)-1
    gsn_polygon(wks,plot,boundslon(edge,:),boundslat(edge,:),resource)
  end do
end

; vim:list
