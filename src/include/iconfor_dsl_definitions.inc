! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! DSL definitions 

#define onCells_3D_tracers  REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define onCells_2D_tracers  REAL(wp), POINTER, DIMENSION(:,:,:)
#define onCells_HalfLevels_tracers REAL(wp), POINTER, DIMENSION(:,:,:,:)

#define onCells_3D_Type(tp) TYPE(tp), POINTER, DIMENSION(:,:,:)
#define onCells_3D_RealPrecision(pr) REAL(pr), POINTER, DIMENSION(:,:,:)
#define onCells_3D_Real REAL(wp),  POINTER, DIMENSION(:,:,:)
#define onCells_3D_sp REAL(sp),  POINTER, DIMENSION(:,:,:)
#define onCells_HalfLevels_tracers REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define onCells_HalfLevels         REAL(wp), POINTER, DIMENSION(:,:,:)


#define onCells_2D_Type(tp) TYPE(tp), POINTER, DIMENSION(:,:)
#define onCells_2D_Type(tp) TYPE(tp), POINTER, DIMENSION(:,:)
#define onCells_2D_Real REAL(wp), POINTER, DIMENSION(:,:)
#define onCells_2D_RealPrecision(pr) REAL(pr), POINTER,  DIMENSION(:,:)
#define onCells_3D_Int INTEGER, POINTER, DIMENSION(:,:,:)
#define onCells_2D_Int INTEGER, POINTER, DIMENSION(:,:)

#define onCells_1D_Real REAL(wp), POINTER, DIMENSION(:)
#define onCells_HalfLevels1D_Real REAL(wp), POINTER, DIMENSION(:)

#define onEdges_3D_tracers  REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define onEdges_3D_RealPrecision(pr) REAL(pr), POINTER, DIMENSION(:,:,:)
#define onEdges_3D_Real REAL(wp), POINTER, DIMENSION(:,:,:)
#define onEdges_3D_Int INTEGER, POINTER, DIMENSION(:,:,:)
#define onEdges_3D_Int INTEGER, POINTER, DIMENSION(:,:,:)
#define onEdges_3D_Type(tp) TYPE(tp), POINTER, DIMENSION(:,:,:)

#define onEdges_HalfLevels REAL(wp), POINTER, DIMENSION(:,:,:)
#define onEdges_2D_Real REAL(wp), POINTER, DIMENSION(:,:)
#define onEdges_2D_RealPrecision(pr) REAL(pr), POINTER,  DIMENSION(:,:)
#define onEdges_2D_Int INTEGER, POINTER, DIMENSION(:,:)

#define onVertices_3D_Real REAL(wp), POINTER, DIMENSION(:,:,:)
#define onVertices_3D_Int INTEGER, POINTER, DIMENSION(:,:,:)
#define onVertices_3D_Type(tp) TYPE(tp), POINTER, DIMENSION(:,:,:)

#define onVertices_2D_Real REAL(wp), POINTER, DIMENSION(:,:)

!---------------------
! block definitions
#define onCellsBlock_3D_Real REAL(wp),  POINTER, DIMENSION(:,:)
#define onCellsBlock_2D_Real REAL(wp),  POINTER, DIMENSION(:)
#define blockList_Type(tp)  TYPE(tp), POINTER, DIMENSION(:)

!---------------------
! mappings
#define mapEdgesToCells_3D    REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define mapEdgesToVertices_3D REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define mapEdgesToEdges_3D    REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define mapEdgesToEdges_2D    REAL(wp), POINTER, DIMENSION(:,:,:)

#define mapCellsToEdges_3D    REAL(wp), POINTER, DIMENSION(:,:,:,:)
#define mapCellsToEdges_3D_RealPrecision(pr)  REAL(pr), POINTER, DIMENSION(:,:,:,:)
#define mapCellsToEdges_2D    REAL(wp), POINTER, DIMENSION(:,:,:)
#define mapCellsToEdges_2D_RealPrecision(pr)  REAL(pr), POINTER, DIMENSION(:,:,:)
#define mapCellsToCells_2D    REAL(wp), POINTER, DIMENSION(:,:,:)

#define mapEdgesToCells_3D_RealPrecision(pr) REAL(pr), POINTER, DIMENSION(:,:,:,:)
#define mapEdgesToCells_2D_RealPrecision(pr) REAL(pr), POINTER, DIMENSION(:,:,:)
#define mapEdgesToEdges_2D_RealPrecision(pr) REAL(pr), POINTER, DIMENSION(:,:,:)

#define mapEdgesToCells_3D_Type(tp)     TYPE(tp), POINTER, DIMENSION(:,:,:,:)
#define mapEdgesToVertices_3D_Type(tp)  TYPE(tp), POINTER, DIMENSION(:,:,:,:)


!---------------------
! connectivity
#define onVertices_3D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:,:)
#define onCells_3D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:,:)
#define onEdges_3D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:,:)

#define onVertices_2D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:)
#define onCells_2D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:)
#define onEdges_2D_Connectivity INTEGER, POINTER, DIMENSION(:,:,:)


!---------------------
! generic types
#define onGrid_1D_Real REAL(wp), POINTER, DIMENSION(:)
#define onGrid_1D_Int   INTEGER, POINTER, DIMENSION(:)
#define onGrid_HalfLevels1D_Real REAL(wp), POINTER, DIMENSION(:)

!---------------------
! shortcuts

#define onCells_3D onCells_3D_Real 
#define onEdges_3D onEdges_3D_Real
#define onVertices_3D onVertices_3D_Real

#define onCells_2D onCells_2D_Real
#define onEdges_2D onEdges_2D_Real
#define onVertices_2D onVertices_2D_Real

#define onCells_1D onCells_1D_Real
#define onCells_HalfLevels1D onCells_HalfLevels1D_Real

#define onCells_Type(tp) onCells_3D_Type(tp)
#define onEdges_Type(tp) onEdges_3D_Type(tp)
#define onVertices_Type(tp) onVertices_3D_Type(tp)

#define onCells onCells_3D
#define onEdges onEdges_3D
#define onVertices onVertices_3D

#define onCells_tracers onCells_3D_tracers
#define onEdges_tracers onEdges_3D_tracers
#define onEdges_Tracers onEdges_3D_tracers

#define onCellsBlock onCellsBlock_3D_Real

#define mapEdgesToCells      mapEdgesToCells_3D
#define mapEdgesToVertices   mapEdgesToVertices_3D
#define mapEdgesToEdges      mapEdgesToEdges_3D
#define mapCellsToEdges      mapCellsToEdges_3D

#define onVertices_Connectivity onVertices_2D_Connectivity
#define onCells_Connectivity onCells_2D_Connectivity
#define onEdges_Connectivity onEdges_2D_Connectivity

#define onGrid_1D onGrid_1D_Real
#define onGrid_HalfLevels1D onGrid_HalfLevels1D_Real
#define mapCellsToEdges mapCellsToEdges_3D

!---------------------
! Upper-lower case
#define OnCells_2D_Int onCells_2D_Int
#define OnEdges_2D_Int onEdges_2D_Int
#define OnCells onCells
#define OnEdges onEdges
#define OnVertices onVertices
#define onCells_3D_Tracers onCells_3D_tracers
#define onCells_Tracers onCells_3D_tracers
