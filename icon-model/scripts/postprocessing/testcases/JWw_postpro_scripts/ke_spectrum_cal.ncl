; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;-----------------------------------------
; Hui Wan (MPI-M, 2009-03)
;-----------------------------------------
begin

  Model    = getenv("Model")
  SpecData = getenv("SpecData") 
  DataPath = getenv("DataPath")
  ExpName  = getenv("ExpName")
  OutFrq   = stringtoint(getenv("output_frequency"))

  vorName  = getenv("vor_coeff_name") 
  divName  = getenv("div_coeff_name") 

  day = (/10,30/)
  ndays = dimsizes(day)

  if (Model.eq."ECHAM") then
     idx = day*OutFrq -1
  else
     idx = day*OutFrq
  end if

  re = 6371229.

;------------------------------------------
; find out the resolution
;------------------------------------------

  File = addfile( SpecData+"_spec.nc","r" )

 ; shape of the arrays in the input file: (ntime,nlev,nspec,2)

  nlev = dimsizes(File->$vorName$(0,:,0,0)) 
  nspc = dimsizes(File->$vorName$(0,0,:,0))

  trunctmp = sqrt(8*nspc+1) -3
  trunc    = floattointeger(trunctmp/2)

;---------------------------------------------------------------
; check whether the simulation is long enough
;---------------------------------------------------------------

  nsteps_in_data_file = dimsizes(File->time)

  if (Model.eq."ECHAM") then
     ndays_in_data_file  = (nsteps_in_data_file -1.)/OutFrq
  else
     ndays_in_data_file  = nsteps_in_data_file/OutFrq
  end if

  nmax_needed = max(day)

  if ( ndays_in_data_file .lt. nmax_needed ) then
     print("Problem encountered: found only "+ ndays_in_data_file +\
           " day(s) of data in the model output.")
     print(nmax_needed +" days are needed for the diagnosis of "  \
          +"the kinetic energy spectrum.")
     exit
  end if

;------------------------------------------
; calculate spectrum
;------------------------------------------
 ; tmp arrays 
 ;------------------------------------------

  svo2d = new( (/nlev,trunc+1,trunc+1,2/), double )
  sd2d  = svo2d

 ;------------------------------------------
 ; time step by time step 
 ;------------------------------------------
  do jday=0,ndays-1

  ;---------------------------------------------------
  ; read vorticity and divergence
  ;---------------------------------------------------

    svo = File->$vorName$(idx(jday),:,:,:)  ; shape of svo: (nlev,nspec,2)
    sd  = File->$divName$(idx(jday),:,:,:)  ; shape of sd : (nlev,nspec,2)

  ;---------------------------------------------------
  ; re-shape the coefficient arrays
  ;---------------------------------------------------
    ii=0

    do im=0,trunc

       if (im.ge.1) then
       svo2d(:,0:im-1,im,:) = 0.
       end if

       do in=im,trunc
          svo2d(:,in,im,:) = svo(:,ii,:)
          sd2d (:,in,im,:) =  sd(:,ii,:)
          ii = ii+1
       end do
    end do

  ;---------------------------------------------------
  ; the array containing the result
  ;---------------------------------------------------

    spec_ke     = new( (/nlev,trunc+1/), double )
    spec_ke!0   = "lev"
    spec_ke&lev = svo&lev
    spec_ke!1   = "wavenumber"

    spec_div    = spec_ke
    spec_vor    = spec_ke

    spec_ke@long_name = "kinetic energy spectrum"
    spec_div@long_name = "kinetic energy spectrum (divergent part)"
    spec_vor@long_name = "kinetic energy spectrum (rotational part)"

  ;---------------------------------------------------
  ; calculate spectra
  ;---------------------------------------------------

    do in=1,trunc

       spec_vor(:,in) = re*re/4./in/(in+1) \
                                                \
                       *(  svo2d(:,in,0,0)^2    \
                         + svo2d(:,in,0,1)^2    \
                                                \
                         + 2*dim_sum( svo2d(:,in,1:,0)^2 \
                                     +svo2d(:,in,1:,1)^2) )

       spec_div(:,in) = re*re/4./in/(in+1) \
                                                \
                       *(  sd2d(:,in,0,0)^2    \
                         + sd2d(:,in,0,1)^2    \
                                                \
                         + 2*dim_sum( sd2d(:,in,1:,0)^2 \
                                     +sd2d(:,in,1:,1)^2) )


       spec_ke (:,in) =  spec_vor(:,in) + spec_div(:,in)

    end do

  ;------------------------------------------------------
  ; output
  ;------------------------------------------------------
   out_file_name = DataPath+ExpName+"_KE-spectrum_day"+day(jday)+".nc"

  ; if an old file with the same name exists, remove it

   istat = systemfunc("test -e "+out_file_name+" ; echo $?")
   if (istat.eq.0) then
       system("rm "+ out_file_name)
   end if

   out_file = addfile(out_file_name,"c") 
   out_file->ke  = spec_ke
   out_file->div = spec_div
   out_file->rot = spec_vor

  end do

end
