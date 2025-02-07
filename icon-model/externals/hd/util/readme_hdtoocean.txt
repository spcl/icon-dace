! readme_hdtoocean.txt - Readme on how to convert HD discharges to an ocean grid 
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: CC-BY-4.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________


1. Run frun_hdtoocean.com to generate coupling file
2. Run convert_hdtoocean.com to convert the discharge output


Conversion of HD mouth points into NEMO mouth points
  Method INEMOU = 3    Use both an existing mask and the coastal ocean points derived from ocean model sea mask

  For a specific HD mouth point, first, the closest point in a prescribed mouth mask is searched within a search radius DISMAX1.
  If no point is found, closest coastal point is searched within a search radius DISMAX2, whereat the coastal points where derived from the NEMO sea mask.

0.5 degree
        DISMAX1=25000.    ! 25 km for primary mask
        DISMAX2=100000.   ! 100 km for secondary mask (default for only 1 mask)
5 Min.
        DISMAX1=4000.     ! 4 km for primary mask
        DISMAX2=200000.   ! normally, 16 km for secondary mask (default for only 1 mask)
                            but NEMO ocean coast too smooth --> 200 km 

