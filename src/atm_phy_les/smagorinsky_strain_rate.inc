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

! see mo_sgs_turbmetric for an explanation
vn_vert1  = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
            p_patch%edges%primal_normal_vert(je,jb,1)%v1 + &
            v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
            p_patch%edges%primal_normal_vert(je,jb,1)%v2

vn_vert2  = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
            p_patch%edges%primal_normal_vert(je,jb,2)%v1 + &
            v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
            p_patch%edges%primal_normal_vert(je,jb,2)%v2

vt_vert2  = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
            p_patch%edges%dual_normal_vert(je,jb,2)%v1 + &
            v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * &
            p_patch%edges%dual_normal_vert(je,jb,2)%v2

vt_vert1  = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
            p_patch%edges%dual_normal_vert(je,jb,1)%v1 + &
            v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * &
            p_patch%edges%dual_normal_vert(je,jb,1)%v2

vn_vert3  = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
            p_patch%edges%primal_normal_vert(je,jb,3)%v1 + &
            v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
            p_patch%edges%primal_normal_vert(je,jb,3)%v2

vn_vert4  = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
            p_patch%edges%primal_normal_vert(je,jb,4)%v1 + &
            v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
            p_patch%edges%primal_normal_vert(je,jb,4)%v2

vt_vert4  = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
            p_patch%edges%dual_normal_vert(je,jb,4)%v1 + &
            v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * &
            p_patch%edges%dual_normal_vert(je,jb,4)%v2

vt_vert3  = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
            p_patch%edges%dual_normal_vert(je,jb,3)%v1 + &
            v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * &
            p_patch%edges%dual_normal_vert(je,jb,3)%v2

!W at full levels
w_full_c1 = 0.5_wp * (                                        &
           p_nh_prog%w(iecidx(je,jb,1), jk, iecblk(je,jb,1)) + &
           p_nh_prog%w(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )

w_full_c2 = 0.5_wp * (                                        &
           p_nh_prog%w(iecidx(je,jb,2),jk,  iecblk(je,jb,2)) + &
           p_nh_prog%w(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) )

!W at full levels vertices from w at vertices at interface levels
w_full_v1 = 0.5_wp * (                                    &
           w_vert(ividx(je,jb,1), jk, ivblk(je,jb,1)) +    &
           w_vert(ividx(je,jb,1),jk+1,ivblk(je,jb,1)) )

w_full_v2 = 0.5_wp * (                                    &
           w_vert(ividx(je,jb,2), jk, ivblk(je,jb,2)) +    &
           w_vert(ividx(je,jb,2),jk+1,ivblk(je,jb,2)) )

!Strain rates at edge center
D_11      =  2._wp * (vn_vert4-vn_vert3) * &
           p_patch%edges%inv_vert_vert_length(je,jb) &
              - 2._wp*(vn_ie(je,jk,jb)-vn_ie(je,jk+1,jb))*     &
                p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
                p_nh_metrics%ddxn_z_full(je,jk,jb)

D_12       =  p_patch%edges%tangent_orientation(je,jb) *  &
              (vn_vert2-vn_vert1) * p_patch%edges%inv_primal_edge_length(je,jb)&
            + (vt_vert4-vt_vert3) * p_patch%edges%inv_vert_vert_length(je,jb) &
            - (vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb)) * &
              p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
              p_nh_metrics%ddxn_z_full(je,jk,jb)         &
            - (vn_ie(je,jk,jb)-vn_ie(je,jk+1,jb))*     &
              p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
              p_nh_metrics%ddxt_z_full(je,jk,jb)

! todo: vertical metric term multiplicator is missing
D_13       =   (vn_ie(je,jk,jb)-vn_ie(je,jk+1,jb)) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) &
             + (w_full_c2 - w_full_c1) * p_patch%edges%inv_dual_edge_length(je,jb)     &
             - (w_ie(je,jk,jb) - w_ie(je,jk+1,jb)) *         &
               p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * &
               p_nh_metrics%ddxn_z_full(je,jk,jb)

D_22       =   2._wp*(vt_vert2-vt_vert1)*p_patch%edges%tangent_orientation(je,jb) * &
               p_patch%edges%inv_primal_edge_length(je,jb)    &
             + 2._wp*(vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb)) *      &
               p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * p_nh_metrics%ddxt_z_full(je,jk,jb)

! todo: vertical metric term multiplicator is missing
D_23       =   (vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb)) &
  &            * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) &
             + p_patch%edges%tangent_orientation(je,jb) * (w_full_v2 - w_full_v1) &
               * p_patch%edges%inv_primal_edge_length(je,jb)    &
             + (w_ie(je,jk,jb) - w_ie(je,jk+1,jb)) &
               * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) &
               * p_nh_metrics%ddxt_z_full(je,jk,jb)

! todo: vertical metric term multiplicator is missing
D_33       =   2._wp * (w_ie(je,jk,jb) - w_ie(je,jk+1,jb)) &
               * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

!Mechanical prod is half of this value divided by km
shear(je,jk,jb) = D_11**2 + D_22**2 + D_33**2 &
  &         + 2._wp * ( D_12**2 + D_13**2 + D_23**2 )

!calculate divergence to get the deviatoric part of stress tensor in
!diffusion: D_11-1/3*(D_11+D_22+D_33)
div_of_stress(je,jk,jb)  = 0.5_wp * ( D_11 + D_22 + D_33 )
