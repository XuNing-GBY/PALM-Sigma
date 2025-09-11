!> @file diffusion_w.f90
!--------------------------------------------------------------------------------------------------!
! This file was modified to work as part of PALM-Sigma by Xu Ning in 2025.
!
! PALM-Sigma consists of a series of user-defined modules and modified PALM modules, developed for
! wave-phase-resolved LES. PALM-Sigma should be used as USER_CODE based on the framework of PALM.
!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the terms of the GNU General
! Public License as published by the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
! Public License for more details.
!
! You should have received a copy of the GNU General Public License along with PALM. If not, see
! <http://www.gnu.org/licenses/>.
!
! Copyright 2025 Xu Ning
! Copyright 1997-2021 Leibniz Universitaet Hannover
!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Diffusion term of the w-component
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_w_mod


    PRIVATE
    PUBLIC diffusion_w

    INTERFACE diffusion_w
       MODULE PROCEDURE diffusion_w
       MODULE PROCEDURE diffusion_w_ij
    END INTERFACE diffusion_w

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w

       USE arrays_3d,                                                                              &
           ONLY :  ddzu, ddzw, drho_air_zw, km, rho_air, tend, u, v, w

       USE grid_variables,                                                                         &
           ONLY :  ddx, ddy

       USE indices,                                                                                &
           ONLY :  nxl, nxr, nyn, nys, nzb, nzt, topo_flags

       USE kinds

       USE surface_mod,                                                                            &
           ONLY :  surf_def_v, surf_lsm_v, surf_usm_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  kmxm              !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp              !<diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym              !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp              !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_east         !< flag to mask vertical wall east of the grid point
       REAL(wp) ::  mask_north        !< flag to mask vertical wall north of the grid point
       REAL(wp) ::  mask_south        !< flag to mask vertical wall south of the grid point
       REAL(wp) ::  mask_west         !< flag to mask vertical wall west of the grid point



       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, l, m) &
       !$ACC PRIVATE(surf_e, surf_s, flag, kmxm, kmxp, kmym, kmyp) &
       !$ACC PRIVATE(mask_west, mask_east, mask_south, mask_north) &
       !$ACC PRESENT(topo_flags, km) &
       !$ACC PRESENT(u, v, w) &
       !$ACC PRESENT(ddzu, ddzw, rho_air, drho_air_zw) &
       !$ACC PRESENT(surf_def_v(0:3)) &
       !$ACC PRESENT(surf_lsm_v(0:3)) &
       !$ACC PRESENT(surf_usm_v(0:3)) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt-1
!
!--             Predetermine flag to mask topography and wall-bounded grid points.
                flag       = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i),   3 ) )
                mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 3 ) )
                mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 3 ) )
                mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 3 ) )
                mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 3 ) )
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * ( km(k,j,i)   + km(k,j,i+1)   +                                   &
                                   km(k+1,j,i) + km(k+1,j,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i)   + km(k,j,i-1)   +                                   &
                                   km(k+1,j,i) + km(k+1,j,i-1) )
                kmyp = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +                                   &
                                   km(k,j+1,i) + km(k+1,j+1,i) )
                kmym = 0.25_wp * ( km(k,j,i)   + km(k+1,j,i)   +                                   &
                                   km(k,j-1,i) + km(k+1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                                          &
                              + ( mask_east *  kmxp * (                                            &
                                          ( w(k,j,i+1)   - w(k,j,i)   ) * ddx                      &
                                        + ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzu(k+1)                &
                                                      )                                            &
                                - mask_west * kmxm *  (                                            &
                                          ( w(k,j,i)     - w(k,j,i-1) ) * ddx                      &
                                        + ( u(k+1,j,i)   - u(k,j,i)   ) * ddzu(k+1)                &
                                                      )                                            &
                                ) * ddx                                 * flag                     &
                              + ( mask_north * kmyp * (                                            &
                                          ( w(k,j+1,i)   - w(k,j,i)   ) * ddy                      &
                                        + ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzu(k+1)                &
                                                      )                                            &
                                - mask_south * kmym * (                                            &
                                          ( w(k,j,i)     - w(k,j-1,i) ) * ddy                      &
                                        + ( v(k+1,j,i)   - v(k,j,i)   ) * ddzu(k+1)                &
                                                      )                                            &
                                ) * ddy                                 * flag                     &
                              + 2.0_wp * (                                                         &
                                km(k+1,j,i) * ( w(k+1,j,i) - w(k,j,i) )   * ddzw(k+1)              &
                                            * rho_air(k+1)                                         &
                              - km(k,j,i)   * ( w(k,j,i)   - w(k-1,j,i) ) * ddzw(k)                &
                                            * rho_air(k)                                           &
                                         ) * ddzu(k+1) * drho_air_zw(k) * flag
             ENDDO

!
!--          Add horizontal momentum flux v'w' at north- (l=0) and south-facing (l=1) surfaces.
!--          Note, in the the flat case, loops won't be entered as start_index > end_index.
!--          Furtermore, note, no vertical natural surfaces so far.
!--          Default-type surfaces
             DO  l = 0, 1
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_w(m) * ddy
                ENDDO
             ENDDO
!
!--          Natural-type surfaces
             DO  l = 0, 1
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_w(m) * ddy
                ENDDO
             ENDDO
!
!--          Urban-type surfaces
             DO  l = 0, 1
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_w(m) * ddy
                ENDDO
             ENDDO
!
!--          Add horizontal momentum flux u'w' at east- (l=2) and west-facing (l=3) surface.
!--          Default-type surfaces
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_w(m) * ddx
                ENDDO
             ENDDO
!
!--          Natural-type surfaces
             DO  l = 2, 3
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_w(m) * ddx
                ENDDO
             ENDDO
!
!--          Urban-type surfaces
             DO  l = 2, 3
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_w(m) * ddx
                ENDDO
             ENDDO

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_w


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_w_ij( i, j )

       ! BOC - Xu
       USE arrays_3d,                                                                              &
           ONLY :  ddzu, ddzw, drho_air_zw, km, rho_air, tend, u, v, w, ddzzu, ddzzw
       ! EOC - Xu

       USE grid_variables,                                                                         &
           ONLY :  ddx, ddy

       USE indices,                                                                                &
           ONLY :  nzb, nzt, topo_flags

       USE kinds

       USE surface_mod,                                                                            &
           ONLY :  surf_def_v, surf_lsm_v, surf_usm_v

       ! BOC - Xu   
       USE wind_wave_inta_mod,                                                                     &
           ONLY:  dzeta_dt_s, dzeta_dt_u, dzeta_dt_v, dzeta_dt_w,                                  &
                  dzeta_dx_s, dzeta_dx_u, dzeta_dx_v, dzeta_dx_w,                                  &
                  dzeta_dy_s, dzeta_dy_u, dzeta_dy_v, dzeta_dy_w,                                  &
                  dzeta_dz_u, dzeta_dz_v, dzeta_dz_s
       ! EOC - Xu

       IMPLICIT NONE


       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  kmxm              !< diffusion coefficient on leftward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmxp              !< diffusion coefficient on rightward side of the w-gridbox - interpolated onto xu-y grid
       REAL(wp) ::  kmym              !< diffusion coefficient on southward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  kmyp              !< diffusion coefficient on northward side of the w-gridbox - interpolated onto x-yv grid
       REAL(wp) ::  mask_east         !< flag to mask vertical wall east of the grid point
       REAL(wp) ::  mask_north        !< flag to mask vertical wall north of the grid point
       REAL(wp) ::  mask_south        !< flag to mask vertical wall south of the grid point
       REAL(wp) ::  mask_west         !< flag to mask vertical wall west of the grid point


       ! BOC - Xu
       DO  k = nzb+1, nzt-1
!
!--       Predetermine flag to mask topography and wall-bounded grid points.
          flag       = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i),   3 ) )
          mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 3 ) )
          mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 3 ) )
          mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 3 ) )
          mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 3 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k+1,j,i)+km(k+1,j,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k+1,j,i)+km(k+1,j,i-1) )
          kmyp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j+1,i)+km(k+1,j+1,i) )
          kmym = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )

          tend(k,j,i) = tend(k,j,i)                                                                &
                      + ( kmxp        * ( ( u(k+1,j,i+1) - u(k,j,i+1) ) * ddzzu(k+1)               &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i+1) ) * 0.5_wp     &
                                        + ( w(k,j,i+1) - w(k,j,i) ) * ddx                          &
                                        + ( ( w(k+1,j,i+1) + w(k+1,j,i) + w(k,j,i+1) + w(k,j,i) )  &
                                          - ( w(k-1,j,i+1) + w(k-1,j,i) + w(k,j,i+1) + w(k,j,i) )  &
                                          ) * 0.25_wp * ddzzu(k+1)                                 &
                                          * ( dzeta_dx_w(k,j,i) + dzeta_dx_w(k,j,i+1) ) * 0.5_wp   &
                                        ) * mask_east                                              &
                        - kmxm        * ( ( u(k+1,j,i) - u(k,j,i) ) * ddzzu(k+1)                   &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i-1) ) * 0.5_wp     &
                                        + ( w(k,j,i) - w(k,j,i-1) ) * ddx                          &
                                        + ( ( w(k+1,j,i-1) + w(k+1,j,i) + w(k,j,i-1) + w(k,j,i) )  &
                                          - ( w(k-1,j,i-1) + w(k-1,j,i) + w(k,j,i-1) + w(k,j,i) )  &
                                          ) * 0.25_wp * ddzzu(k+1)                                 &
                                          * ( dzeta_dx_w(k,j,i) + dzeta_dx_w(k,j,i-1) ) * 0.5_wp   &
                                        ) * mask_west                                              &
                        ) * ddx * flag                                                             &
                      + ( km(k+1,j,i) * ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k+2,j,i+1) + u(k+2,j,i)  &
                                          - u(k+1,j,i+1) - u(k+1,j,i) - u(k,j,i+1)   - u(k,j,i)    &
                                          ) * 0.25_wp * ddzzw(k+1)                                 &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        + ( w(k+1,j,i+1) + w(k+1,j,i) + w(k,j,i+1) + w(k,j,i)      &
                                          - w(k+1,j,i-1) - w(k+1,j,i) - w(k,j,i-1) - w(k,j,i)      &
                                          ) * 0.25_wp * ddx                                        &
                                        + ( w(k+1,j,i) - w(k,j,i) ) * ddzzw(k+1)                   &
                                          * ( dzeta_dx_w(k+1,j,i) + dzeta_dx_w(k,j,i) ) * 0.5_wp   &
                                        ) * rho_air(k+1)                                           &
                        - km(k,j,i)   * ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1)   + u(k,j,i)    &
                                          - u(k,j,i+1)   - u(k,j,i)   - u(k-1,j,i+1) - u(k-1,j,i)  &
                                          ) * 0.25_wp * ddzzw(k)                                   &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        + ( w(k-1,j,i+1) + w(k-1,j,i) + w(k,j,i+1) + w(k,j,i)      &
                                          - w(k-1,j,i-1) - w(k-1,j,i) - w(k,j,i-1) - w(k,j,i)      &
                                          ) * 0.25_wp * ddx                                        &
                                        + ( w(k,j,i) - w(k-1,j,i) ) * ddzzw(k)                     &
                                          * ( dzeta_dx_w(k-1,j,i) + dzeta_dx_w(k,j,i) ) * 0.5_wp   &
                                        ) * rho_air(k)                                             &
                        ) * ddzzu(k+1) * drho_air_zw(k) * dzeta_dx_w(k,j,i) * flag                 &
                      + ( kmyp        * ( ( v(k+1,j+1,i) - v(k,j+1,i) ) * ddzzu(k+1)               &
                                          * ( dzeta_dz_s(j+1,i) + dzeta_dz_s(j,i) ) * 0.5_wp     &
                                        + ( w(k,j+1,i) - w(k,j,i) ) * ddy                          &
                                        + ( w(k+1,j+1,i) + w(k+1,j,i) + w(k,j+1,i) + w(k,j,i)      &
                                          - w(k-1,j+1,i) - w(k-1,j,i) - w(k,j+1,i) - w(k,j,i)      &
                                          ) * 0.25_wp * ddzzu(k+1)                                 &
                                          * ( dzeta_dy_w(k,j+1,i) + dzeta_dy_w(k,j,i) ) * 0.5_wp   &
                                        ) * mask_north                                             &
                        - kmym        * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzzu(k+1)                   &
                                          * ( dzeta_dz_s(j-1,i) + dzeta_dz_s(j,i) ) * 0.5_wp     &
                                        + ( w(k,j,i) - w(k,j-1,i) ) * ddy                          &
                                        + ( w(k+1,j-1,i) + w(k+1,j,i) + w(k,j-1,i) + w(k,j,i)      &
                                          - w(k-1,j-1,i) - w(k-1,j,i) - w(k,j-1,i) - w(k,j,i)      &
                                          ) * 0.25_wp * ddzzu(k+1)                                 &
                                          * ( dzeta_dy_w(k,j-1,i) + dzeta_dy_w(k,j,i) ) * 0.5_wp   &
                                        ) * mask_south                                             &
                        ) * ddy * flag                                                             &
                      + ( km(k+1,j,i) * ( ( v(k+2,j+1,i) + v(k+2,j,i) + v(k+1,j+1,i) + v(k+1,j,i)  &
                                          - v(k+1,j+1,i) - v(k+1,j,i) - v(k,j+1,i) - v(k,j,i)      &
                                          ) * 0.25_wp * ddzzw(k+1)                                 &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        + ( w(k+1,j+1,i) + w(k+1,j,i) + w(k,j+1,i) + w(k,j,i)      &
                                          - w(k+1,j-1,i) - w(k+1,j,i) - w(k,j-1,i) - w(k,j,i)      &
                                          ) * 0.25_wp * ddy                                        &
                                        + ( w(k+1,j,i) - w(k,j,i) ) * ddzzw(k+1)                   &
                                          * ( dzeta_dy_w(k+1,j,i) + dzeta_dy_w(k,j,i) ) * 0.5_wp   &
                                        ) * rho_air(k+1)                                           &
                        - km(k,j,i)   * ( ( v(k+1,j+1,i) + v(k+1,j,i) + v(k,j+1,i) + v(k,j,i)      &
                                          - v(k-1,j+1,i) - v(k-1,j,i) - v(k,j+1,i) - v(k,j,i)      &
                                          ) * 0.25_wp * ddzzw(k)                                   &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        + ( w(k-1,j+1,i) + w(k-1,j,i) + w(k,j+1,i) + w(k,j,i)      &
                                          - w(k-1,j-1,i) - w(k-1,j,i) - w(k,j-1,i) - w(k,j,i)      &
                                          ) * 0.25_wp * ddy                                        &
                                        + ( w(k,j,i) - w(k-1,j,i) ) * ddzzw(k)                     &
                                          * ( dzeta_dy_w(k-1,j,i) + dzeta_dy_w(k,j,i) ) * 0.5_wp   &
                                        ) * rho_air(k)                                             &
                        ) * ddzzu(k+1) * drho_air_zw(k) * dzeta_dy_w(k,j,i) * flag                 &
                      + ( km(k+1,j,i) * ( ( w(k+1,j,i) - w(k,j,i) ) * ddzzw(k+1)                   &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        ) * rho_air(k+1)                                           &
                        - km(k,j,i)   * ( ( w(k,j,i) - w(k-1,j,i) ) * ddzzw(k)                     &
                                          * ( dzeta_dz_s(j,i) + dzeta_dz_s(j,i) ) * 0.5_wp       &
                                        ) * rho_air(k)                                             &
                        ) * 2.0_wp * ddzzu(k+1) * drho_air_zw(k) * dzeta_dz_s(j,i) * flag
       ENDDO
!
!--    Add horizontal momentum flux v'w' at north- (l=0) and south-facing (l=1) surfaces. Note, in
!--    the the flat case, loops won't be entered as start_index > end_index. Furtermore, note, no
!--    vertical natural surfaces so far.
!--    Default-type surfaces
       DO  l = 0, 1
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_w(m) * ddy
          ENDDO
       ENDDO
!
!--    Natural-type surfaces
       DO  l = 0, 1
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_w(m) * ddy
          ENDDO
       ENDDO
!
!--    Urban-type surfaces
       DO  l = 0, 1
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_w(m) * ddy
          ENDDO
       ENDDO
!
!--    Add horizontal momentum flux u'w' at east- (l=2) and west-facing (l=3) surfaces.
!--    Default-type surfaces
       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_w(m) * ddx
          ENDDO
       ENDDO
!
!--    Natural-type surfaces
       DO  l = 2, 3
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_w(m) * ddx
          ENDDO
       ENDDO
!
!--    Urban-type surfaces
       DO  l = 2, 3
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_w(m) * ddx
          ENDDO
       ENDDO
       ! EOC - Xu

    END SUBROUTINE diffusion_w_ij

 END MODULE diffusion_w_mod
