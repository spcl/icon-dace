! ICON-Land
!
! ---------------------------------------
! Copyright (C) 2013-2024, MPI-M, MPI-BGC
!
! Contact: icon-model.org
! Authors: AUTHORS.md
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------

PROGRAM example

  dsl4jsb_Use_processes SOIL_, &
    & VEG_
  dsl4jsb_Use_memory(ATM2LAND_)
  dsl4jsb_Def_memory(SOIL_)
  dsl4jsb_Def_memory_tile(SOIL_, thisTile)
  dsl4jsb_get_memory(SOIL_) ! comment
  dsl4jsb_get_memory_tile(SOIL_, currentChild) ! get soil memory from tile source
  var3(:) = a1 * dsl4jsb_memory(SOIL_)%var1(:) + a2 * dsl4jsb_memory(SOIL_)%var2(:)
  dsl4jsb_Use_config(SOIL_)
  dsl4jsb_Def_config(SOIL_)
  dsl4jsb_Get_config(SOIL_)
  x = dsl4jsb_Config(SOIL_)%par1 * dsl4jsb_Config(SOIL_)%par2 * ws(:) ! comment
  CALL sub(dssl4jsb_Lctlib , arg2, dssl4jsb_Lctlib_param(abc)) ! This doesn't work!!!
  CALL sub(dsl4jsb_Lctlib, arg2)
  CALL sub(dsl4jsb_Lctlib)
  CALL sub(dsl4jsb_Lctlib_param(abc), dsl4jsb_Lctlib_param(def))
  dsl4jsb_Real2D_onChunk :: &
    & ws1, ws2, &
    & ws3
  dsl4jsb_Get_var2D_onChunk(SOIL_,ws)
  dsl4jsb_Get_var2D_onChunk(SOIL_, ws)
  dsl4jsb_Get_var2d_onChunk_tile(SOIL_, ws, sourceT)
  dsl4jsb_Get_var2d_onChunk_tile_name(SOIL_, ws, sourceT)
  dsl4jsb_var2d_onChunk(SOIL_,ws)
  CALL subroutine( &
    & dsl4jsb_var2D_onChunk (SOIL_, var1   ), &
    & dsl4jsb_var2D_onChunk(SOIL_, var2   )  &
    & )
  dsl4jsb_Real2D_onDomain :: ws
  dsl4jsb_Get_var2D_onDomain(SOIL_,ws)
  dsl4jsb_var2D_onDomain(SOIL_,ws)
  dsl4jsb_var_ptr(SOIL_,ws)
  x = dsl4jsb_Config(SOIL_)%par1 * dsl4jsb_var_ptr(SOIL_, ws)(:) ! comment
  dsl4jsb_Real3D_onChunk :: ws_l
  dsl4jsb_Get_var3D_onChunk(SOIL_,ws_l)
  dsl4jsb_var3D_onChunk(SOIL_,ws_l)
  dsl4jsb_Real3D_onDomain :: ws_l
  dsl4jsb_Get_var3D_onDomain(SOIL_,ws_l)
  dsl4jsb_var3D_onDomain(SOIL_,ws_l)
  dsl4jsb_Aggregate_onChunk(SOIL_,ws, aggregator)

  ! Test abort on error
  dsl4jsb_Def_menory(SOIL_)

END PROGRAM example
