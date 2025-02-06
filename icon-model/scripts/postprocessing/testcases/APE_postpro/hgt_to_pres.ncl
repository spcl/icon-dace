; ICON
;
; ------------------------------------------
; Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
; Contact information: icon-model.org
; See AUTHORS.TXT for a list of authors
; See LICENSES/ for license information
; SPDX-License-Identifier: BSD-3-Clause
; ------------------------------------------

;***********************************************************************;
; Procedure : hgt_to_pres                                               ;
;                                                                       ;
; Author: Hui Wan (MPI-M, 2011-04-08)                                   ;
;                                                                       ;
; Returns pressure in mb given array z (geopotential height in km).     ; 
; Z must lie between -0.3 km and 119.7 km                               ;
;                                                                       ;
; This procedure is simply the inverse of procedure gsn_geop_hgt        ;
; in $NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl                   ;
; The height and pressure data used for interpolation are taken from    ;
; "gsn_geop_hgt".                                                       ;
;                                                                       ;
;***********************************************************************;
 
undef("hgt_to_pres")
function hgt_to_pres( z[*]:numeric )
local nsa,psa,zsa,nz,found,i,j
 
begin

  nsa = 53
  psa = new( (/nsa/), float, 1.e36 )
  zsa = new( (/nsa/), float, 1.e36 )

  zsa = (/ -0.3, \                      ; km
            0.0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0, \
            3.5,  4.0,  4.5,  5.0,  5.5,  6.0,  6.5, \
            7.0,  7.5,  8.0,  8.5,  9.0,  9.5, 10.0, \
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, \
           18.0, 19.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, \
           50.0, 60.0, 70.0, 80.0, 84.8, 87.7, 90.6, \
           93.3, 96.1, 97.5,100.4,104.9, \
          110.0,114.2,116.7,119.7/)

  psa = (/ 1050., \                     ; mb (hPa)
           1013.25, 954.61, 898.76, 845.59, 795.01, 746.91, 701.21,    \
           657.80, 616.60, 577.52, 540.48, 505.39, 472.17, 440.75,     \
           411.05, 382.99, 356.51, 331.54, 308.00, 285.84, 264.99,     \
           226.99, 193.99, 165.79, 141.70, 121.11, 103.52, 88.497,     \
           75.652, 64.674, 55.293, 25.492, 11.970, 5.746, 2.871, 1.491,\
           0.798, 0.220, 0.052, 0.010, 0.00485,0.00294,0.000178, \
           0.000108, 0.0000656, 0.0000511, 0.0000310, 0.0000146, \
           0.00000691, 0.00000419, 0.00000327, 0.00000254 /)


  if ( any(z.lt.min(zsa)) .or. any(z.gt.max(zsa))) then
    print("hgt_to_pres: Fatal: The height values do not fall between")
    print(min(zsa) + " km and " + max(zsa) + " km.")
    print("Execution halted.")
    exit
  end if
 
  nz   = dimsizes(z)
  pres = new(nz,float)

  do i = 0,nz-1
    found = False
    j = 0

    do while(.not.found.and.j.le.nsa-2)
     if ( ( z(i) .ge. zsa(j) ) .and. ( z(i) .le. zsa(j+1) ) ) then
       pres(i) = psa(j)* \
                ( psa(j+1)/psa(j) )^( (z(i)-zsa(j))/(zsa(j+1)-zsa(j)) )
       found = True 
     end if
     j = j + 1
    end do
  end do
 
  delete(psa)
  delete(zsa)
 
  return(pres)
 
end

