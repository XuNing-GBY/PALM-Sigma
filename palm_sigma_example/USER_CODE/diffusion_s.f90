!> @file diffusion_s.f90
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
!> Diffusion term of scalar quantities (temperature and water content)
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_s_mod


    PRIVATE
    PUBLIC diffusion_s

    INTERFACE diffusion_s
       MODULE PROCEDURE diffusion_s
       MODULE PROCEDURE diffusion_s_ij
    END INTERFACE diffusion_s

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s( s, s_flux_def_h_up,    s_flux_def_h_down,                              &
                               s_flux_t,                                                           &
                               s_flux_lsm_h_up,    s_flux_lsm_h_down,                              &
                               s_flux_usm_h_up,    s_flux_usm_h_down,                              &
                               s_flux_def_v_north, s_flux_def_v_south,                             &
                               s_flux_def_v_east,  s_flux_def_v_west,                              &
                               s_flux_lsm_v_north, s_flux_lsm_v_south,                             &
                               s_flux_lsm_v_east,  s_flux_lsm_v_west,                              &
                               s_flux_usm_v_north, s_flux_usm_v_south,                             &
                               s_flux_usm_v_east,  s_flux_usm_v_west )

       USE arrays_3d,                                                                              &
           ONLY:  ddzu, ddzw, kh, tend, drho_air, rho_air_zw

       USE control_parameters,                                                                     &
           ONLY: use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                                         &
           ONLY:  ddx, ddx2, ddy, ddy2

       USE indices,                                                                                &
           ONLY:  nxl, nxlg, nxr, nxrg, nyn, nyng, nys, nysg, nzb, nzt, topo_flags

       USE kinds

       USE surface_mod,                                                                            &
           ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  mask_bottom       !< flag to mask vertical upward-facing surface
       REAL(wp) ::  mask_east         !< flag to mask vertical surface east of the grid point
       REAL(wp) ::  mask_north        !< flag to mask vertical surface north of the grid point
       REAL(wp) ::  mask_south        !< flag to mask vertical surface south of the grid point
       REAL(wp) ::  mask_top          !< flag to mask vertical downward-facing surface
       REAL(wp) ::  mask_west         !< flag to mask vertical surface west of the grid point

       REAL(wp), DIMENSION(1:surf_def_h(0)%ns) ::  s_flux_def_h_up    !< flux at horizontal upward-facing default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_h(1)%ns) ::  s_flux_def_h_down  !< flux at horizontal donwward-facing default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(2)%ns) ::  s_flux_def_v_east  !< flux at east-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(0)%ns) ::  s_flux_def_v_north !< flux at north-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(1)%ns) ::  s_flux_def_v_south !< flux at south-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(3)%ns) ::  s_flux_def_v_west  !< flux at west-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_h(0)%ns) ::  s_flux_lsm_h_up    !< flux at horizontal upward-facing natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_h(1)%ns) ::  s_flux_lsm_h_down  !< flux at horizontal downward-facing natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(2)%ns) ::  s_flux_lsm_v_east  !< flux at east-facing vertical natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(0)%ns) ::  s_flux_lsm_v_north !< flux at north-facing vertical natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(1)%ns) ::  s_flux_lsm_v_south !< flux at south-facing vertical natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(3)%ns) ::  s_flux_lsm_v_west  !< flux at west-facing vertical natural-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_h(0)%ns) ::  s_flux_usm_h_up    !< flux at horizontal upward-facing urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_h(1)%ns) ::  s_flux_usm_h_down  !< flux at horizontal downward-facing urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(2)%ns) ::  s_flux_usm_v_east  !< flux at east-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(0)%ns) ::  s_flux_usm_v_north !< flux at north-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(1)%ns) ::  s_flux_usm_v_south !< flux at south-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(3)%ns) ::  s_flux_usm_v_west  !< flux at west-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_def_h(2)%ns) ::  s_flux_t           !< flux at model top

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< treated scalar


       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, m) &
       !$ACC PRIVATE(surf_e, surf_s, flag, mask_top, mask_bottom) &
       !$ACC PRIVATE(mask_north, mask_south, mask_west, mask_east) &
       !$ACC PRESENT(topo_flags, kh) &
       !$ACC PRESENT(s) &
       !$ACC PRESENT(ddzu, ddzw, drho_air, rho_air_zw) &
       !$ACC PRESENT(surf_def_h(0:2), surf_def_v(0:3)) &
       !$ACC PRESENT(surf_lsm_h(0:1), surf_lsm_v(0:3)) &
       !$ACC PRESENT(surf_usm_h(0:1), surf_usm_v(0:3)) &
       !$ACC PRESENT(s_flux_def_h_up, s_flux_def_h_down) &
       !$ACC PRESENT(s_flux_t) &
       !$ACC PRESENT(s_flux_def_v_north, s_flux_def_v_south) &
       !$ACC PRESENT(s_flux_def_v_east, s_flux_def_v_west) &
       !$ACC PRESENT(s_flux_lsm_h_up, s_flux_lsm_h_down) &
       !$ACC PRESENT(s_flux_lsm_v_north, s_flux_lsm_v_south) &
       !$ACC PRESENT(s_flux_lsm_v_east, s_flux_lsm_v_west) &
       !$ACC PRESENT(s_flux_usm_h_up, s_flux_usm_h_down) &
       !$ACC PRESENT(s_flux_usm_v_north, s_flux_usm_v_south) &
       !$ACC PRESENT(s_flux_usm_v_east, s_flux_usm_v_west) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nys,nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb+1, nzt
!
!--             Predetermine flag to mask topography and wall-bounded grid points
                flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--             Predetermine flag to mask wall-bounded grid points, equivalent to former s_outer
!--             array
                mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 0 ) )
                mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 0 ) )
                mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 0 ) )
                mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 0 ) )

                tend(k,j,i) = tend(k,j,i)                                                          &
                                                  + 0.5_wp * (                                     &
                                mask_east  * ( kh(k,j,i) + kh(k,j,i+1) )                           &
                                           * ( s(k,j,i+1) - s(k,j,i)   )                           &
                              - mask_west  * ( kh(k,j,i) + kh(k,j,i-1) )                           &
                                           * ( s(k,j,i)   - s(k,j,i-1) )                           &
                                                             ) * ddx2 * flag                       &
                                                  + 0.5_wp * (                                     &
                                mask_north * ( kh(k,j,i) + kh(k,j+1,i) )                           &
                                           * ( s(k,j+1,i) - s(k,j,i)   )                           &
                              - mask_south * ( kh(k,j,i) + kh(k,j-1,i) )                           &
                                           * ( s(k,j,i)   - s(k,j-1,i) )                           &
                                                             ) * ddy2 * flag
             ENDDO

!
!--          Apply prescribed horizontal wall heatflux where necessary. First, determine start and
!--          end index for respective (j,i)-index. Please note, in the flat case following loop will
!--          not be entered, as surf_s=1 and surf_e=0. Furtermore, note, no vertical natural
!--          surfaces so far.
!--          First, for default-type surfaces.
!--          North-facing vertical default-type surfaces
             surf_s = surf_def_v(0)%start_index(j,i)
             surf_e = surf_def_v(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_def_v(0)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_def_v_north(m) * ddy
             ENDDO
!
!--          South-facing vertical default-type surfaces
             surf_s = surf_def_v(1)%start_index(j,i)
             surf_e = surf_def_v(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_def_v(1)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_def_v_south(m) * ddy
             ENDDO
!
!--          East-facing vertical default-type surfaces
             surf_s = surf_def_v(2)%start_index(j,i)
             surf_e = surf_def_v(2)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_def_v(2)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_def_v_east(m) * ddx
             ENDDO
!
!--          West-facing vertical default-type surfaces
             surf_s = surf_def_v(3)%start_index(j,i)
             surf_e = surf_def_v(3)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_def_v(3)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_def_v_west(m) * ddx
             ENDDO
!
!--          Now, for natural-type surfaces.
!--          North-facing
             surf_s = surf_lsm_v(0)%start_index(j,i)
             surf_e = surf_lsm_v(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_lsm_v(0)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_north(m) * ddy
             ENDDO
!
!--          South-facing
             surf_s = surf_lsm_v(1)%start_index(j,i)
             surf_e = surf_lsm_v(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_lsm_v(1)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_south(m) * ddy
             ENDDO
!
!--          East-facing
             surf_s = surf_lsm_v(2)%start_index(j,i)
             surf_e = surf_lsm_v(2)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_lsm_v(2)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_east(m) * ddx
             ENDDO
!
!--          West-facing
             surf_s = surf_lsm_v(3)%start_index(j,i)
             surf_e = surf_lsm_v(3)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_lsm_v(3)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_west(m) * ddx
             ENDDO
!
!--          Now, for urban-type surfaces.
!--          North-facing
             surf_s = surf_usm_v(0)%start_index(j,i)
             surf_e = surf_usm_v(0)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_usm_v(0)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_north(m) * ddy
             ENDDO
!
!--          South-facing
             surf_s = surf_usm_v(1)%start_index(j,i)
             surf_e = surf_usm_v(1)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_usm_v(1)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_south(m) * ddy
             ENDDO
!
!--          East-facing
             surf_s = surf_usm_v(2)%start_index(j,i)
             surf_e = surf_usm_v(2)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_usm_v(2)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_east(m) * ddx
             ENDDO
!
!--          West-facing
             surf_s = surf_usm_v(3)%start_index(j,i)
             surf_e = surf_usm_v(3)%end_index(j,i)
             DO  m = surf_s, surf_e
                k           = surf_usm_v(3)%k(m)
                tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_west(m) * ddx
             ENDDO

!
!--          Compute vertical diffusion. In case that surface fluxes have been prescribed or
!--          computed at bottom and/or top, index k starts/ends at nzb+2 or nzt-1, respectively.
!--          Model top is also mask if top flux is given.
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 0 is used to mask
!--             topography in general, and flag 8 implies information about use_surface_fluxes.
!--             Flag 9 is used to control flux at model top.
                mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
                mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) ) *           &
                              MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )
                flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

                tend(k,j,i) = tend(k,j,i)                                                          &
                                       + 0.5_wp * (                                                &
                                      ( kh(k,j,i) + kh(k+1,j,i) ) *                                &
                                          ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)                      &
                                                            * rho_air_zw(k)                        &
                                                            * mask_top                             &
                                    - ( kh(k,j,i) + kh(k-1,j,i) ) *                                &
                                          ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)                        &
                                                            * rho_air_zw(k-1)                      &
                                                            * mask_bottom                          &
                                                  ) * ddzw(k) * drho_air(k)                        &
                                                              * flag
             ENDDO

!
!--          Vertical diffusion at horizontal walls.
             IF ( use_surface_fluxes )  THEN
!
!--             Default-type surfaces, upward-facing
                surf_s = surf_def_h(0)%start_index(j,i)
                surf_e = surf_def_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(0)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_def_h_up(m) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Default-type surfaces, downward-facing
                surf_s = surf_def_h(1)%start_index(j,i)
                surf_e = surf_def_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(1)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_def_h_down(m) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Natural-type surfaces, upward-facing
                surf_s = surf_lsm_h(0)%start_index(j,i)
                surf_e = surf_lsm_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_lsm_h(0)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_lsm_h_up(m) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Natural-type surfaces, downward-facing
                surf_s = surf_lsm_h(1)%start_index(j,i)
                surf_e = surf_lsm_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_lsm_h(1)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_lsm_h_down(m) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Urban-type surfaces, upward-facing
                surf_s = surf_usm_h(0)%start_index(j,i)
                surf_e = surf_usm_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_usm_h(0)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_usm_h_up(m) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Urban-type surfaces, downward-facing
                surf_s = surf_usm_h(1)%start_index(j,i)
                surf_e = surf_usm_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_usm_h(1)%k(m)
                   tend(k,j,i) = tend(k,j,i) + s_flux_usm_h_down(m) * ddzw(k) * drho_air(k)
                ENDDO

             ENDIF
!
!--          Vertical diffusion at the last computational gridpoint along z-direction
             IF ( use_top_fluxes )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e

                   k   = surf_def_h(2)%k(m)
                   tend(k,j,i) = tend(k,j,i) + ( - s_flux_t(m) ) * ddzw(k) * drho_air(k)
                ENDDO
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_s

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_s_ij( i, j, s,                                                            &
                               s_flux_def_h_up,    s_flux_def_h_down,                              &
                               s_flux_t,                                                           &
                               s_flux_lsm_h_up,    s_flux_lsm_h_down,                              &
                               s_flux_usm_h_up,    s_flux_usm_h_down,                              &
                               s_flux_def_v_north, s_flux_def_v_south,                             &
                               s_flux_def_v_east,  s_flux_def_v_west,                              &
                               s_flux_lsm_v_north, s_flux_lsm_v_south,                             &
                               s_flux_lsm_v_east,  s_flux_lsm_v_west,                              &
                               s_flux_usm_v_north, s_flux_usm_v_south,                             &
                               s_flux_usm_v_east,  s_flux_usm_v_west )

       ! BOC -Xu
       USE arrays_3d,                                                                              &
           ONLY:  ddzu, ddzw, kh, tend, drho_air, rho_air_zw, ddzzu, ddzzw
       ! EOC - Xu

       USE control_parameters,                                                                     &
           ONLY: use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                                         &
           ONLY:  ddx, ddx2, ddy, ddy2

       USE indices,                                                                                &
           ONLY:  nxlg, nxrg, nyng, nysg, nzb, nzt, topo_flags

       USE kinds

       USE surface_mod,                                                                            &
           ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

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
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp) ::  flag              !< flag to mask topography grid points
       REAL(wp) ::  mask_bottom       !< flag to mask vertical upward-facing surface
       REAL(wp) ::  mask_east         !< flag to mask vertical surface east of the grid point
       REAL(wp) ::  mask_north        !< flag to mask vertical surface north of the grid point
       REAL(wp) ::  mask_south        !< flag to mask vertical surface south of the grid point
       REAL(wp) ::  mask_top          !< flag to mask vertical downward-facing surface
       REAL(wp) ::  mask_west         !< flag to mask vertical surface west of the grid point

       REAL(wp), DIMENSION(1:surf_def_h(1)%ns) ::  s_flux_def_h_down  !< flux at horizontal donwward-facing default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_h(0)%ns) ::  s_flux_def_h_up    !< flux at horizontal upward-facing default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(2)%ns) ::  s_flux_def_v_east  !< flux at east-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(0)%ns) ::  s_flux_def_v_north !< flux at north-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(1)%ns) ::  s_flux_def_v_south !< flux at south-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_def_v(3)%ns) ::  s_flux_def_v_west  !< flux at west-facing vertical default-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_h(0)%ns) ::  s_flux_lsm_h_up    !< flux at horizontal upward-facing natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_h(1)%ns) ::  s_flux_lsm_h_down    !< flux at horizontal downward-facing natural-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(2)%ns) ::  s_flux_lsm_v_east  !< flux at east-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(0)%ns) ::  s_flux_lsm_v_north !< flux at north-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(1)%ns) ::  s_flux_lsm_v_south !< flux at south-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_lsm_v(3)%ns) ::  s_flux_lsm_v_west  !< flux at west-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_h(0)%ns) ::  s_flux_usm_h_up    !< flux at horizontal upward-facing urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_h(1)%ns) ::  s_flux_usm_h_down  !< flux at horizontal downward-facing urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(2)%ns) ::  s_flux_usm_v_east  !< flux at east-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(0)%ns) ::  s_flux_usm_v_north !< flux at north-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(1)%ns) ::  s_flux_usm_v_south !< flux at south-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_usm_v(3)%ns) ::  s_flux_usm_v_west  !< flux at west-facing vertical urban-type surfaces
       REAL(wp), DIMENSION(1:surf_def_h(2)%ns) ::  s_flux_t           !< flux at model top

       REAL(wp), DIMENSION(nzb:nzt+1,nysg:nyng,nxlg:nxrg) ::  s  !< treated scalar

       ! BOC - Xu
!
!--    Compute horizontal diffusion
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography and wall-bounded grid points
          flag = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )
!
!--       Predetermine flag to mask wall-bounded grid points, equivalent to former s_outer array
          mask_west  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 0 ) )
          mask_east  = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 0 ) )
          mask_south = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j-1,i), 0 ) )
          mask_north = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j+1,i), 0 ) )

          mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
          mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) )  *        &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )

!
!--       Finally, determine flag to mask both topography itself as well as wall-bounded grid
!--       points, which will be treated further below

          tend(k,j,i) = tend(k,j,i)                                                                &
                      + ( ( kh(k,j,i+1) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k,j,i+1) - s(k,j,i) ) * ddx                                        &
                          + ( s(k+1,j,i+1) + s(k+1,j,i) + s(k,j,i+1) + s(k,j,i)                    &
                            - s(k-1,j,i+1) - s(k-1,j,i) - s(k,j,i+1) - s(k,j,i)                    &
                            ) * 0.25_wp * ddzzw(k) * dzeta_dx_u(k,j,i+1)                           &
                          ) * mask_east                                                            &
                        - ( kh(k,j,i-1) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k,j,i) - s(k,j,i-1) ) * ddx                                        &
                          + ( s(k+1,j,i-1) + s(k+1,j,i) + s(k,j,i-1) + s(k,j,i)                    &
                            - s(k-1,j,i-1) - s(k-1,j,i) - s(k,j,i-1) - s(k,j,i)                    &
                            ) * 0.25_wp * ddzzw(k) * dzeta_dx_u(k,j,i)                             &
                          ) * mask_west                                                            &
                        ) * ddx * flag                                                             &
                      + ( ( kh(k+1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k+1,j,i+1) + s(k+1,j,i) + s(k,j,i+1) + s(k,j,i)                    &
                            - s(k+1,j,i-1) - s(k+1,j,i) - s(k,j,i-1) - s(k,j,i)                    &
                            ) * 0.25_wp * ddx                                                      &
                          + ( s(k+1,j,i) - s(k,j,i) ) * ddzzu(k+1) * dzeta_dx_w(k,j,i)             &
                          ) * rho_air_zw(k) * mask_top                                             &
                        - ( kh(k-1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k-1,j,i+1) + s(k-1,j,i) + s(k,j,i+1) + s(k,j,i)                    &
                            - s(k-1,j,i-1) - s(k-1,j,i) - s(k,j,i-1) - s(k,j,i)                    &
                            ) * 0.25_wp * ddx                                                      &
                          + ( s(k,j,i) - s(k-1,j,i) ) * ddzzu(k) * dzeta_dx_w(k-1,j,i)             &
                          ) * rho_air_zw(k-1) * mask_bottom                                        &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dx_s(k,j,i) * flag                      &
                      + ( ( kh(k,j+1,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k,j+1,i) - s(k,j,i) ) * ddy                                        &
                          + ( s(k+1,j+1,i) + s(k+1,j,i) + s(k,j+1,i) + s(k,j,i)                    &
                            - s(k-1,j+1,i) - s(k-1,j,i) - s(k,j+1,i) - s(k,j,i)                    &
                            ) * 0.25_wp * ddzzw(k) * dzeta_dy_v(k,j+1,i)                           &
                          ) * mask_north                                                           &
                        - ( kh(k,j-1,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k,j,i) - s(k,j-1,i) ) * ddy                                        &
                          + ( s(k+1,j-1,i) + s(k+1,j,i) + s(k,j-1,i) + s(k,j,i)                    &
                            - s(k-1,j-1,i) - s(k-1,j,i) - s(k,j-1,i) - s(k,j,i)                    &
                            ) * 0.25_wp * ddzzw(k) * dzeta_dy_v(k,j,i)                             &
                          ) * mask_south                                                           &
                        ) * ddy * flag                                                             &
                      + ( ( kh(k+1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k+1,j+1,i) + s(k+1,j,i) + s(k,j+1,i) + s(k,j,i)                    &
                            - s(k+1,j-1,i) - s(k+1,j,i) - s(k,j-1,i) - s(k,j,i)                    &
                            ) * 0.25_wp * ddy                                                      &
                          + ( s(k+1,j,i) - s(k,j,i) ) * ddzzu(k+1) * dzeta_dy_w(k,j,i)             &
                          ) * rho_air_zw(k) * mask_top                                             &
                        - ( kh(k-1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( ( s(k-1,j+1,i) + s(k-1,j,i) + s(k,j+1,i) + s(k,j,i)                    &
                            - s(k-1,j-1,i) - s(k-1,j,i) - s(k,j-1,i) - s(k,j,i)                    &
                            ) * 0.25_wp * ddy                                                      &
                          + ( s(k,j,i) - s(k-1,j,i) ) * ddzzu(k) * dzeta_dy_w(k-1,j,i)             &
                          ) * rho_air_zw(k-1) * mask_bottom                                        &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dy_s(k,j,i) * flag                      &
                      + ( ( kh(k+1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( s(k+1,j,i) - s(k,j,i) ) * ddzzu(k+1) * dzeta_dz_s(j,i)                &
                          * rho_air_zw(k) * mask_top                                               &
                        - ( kh(k-1,j,i) + kh(k,j,i) ) * 0.5_wp *                                   &
                          ( s(k,j,i) - s(k-1,j,i) ) * ddzzu(k) * dzeta_dz_s(j,i)                  &
                          * rho_air_zw(k-1) * mask_bottom                                          &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dz_s(j,i) * flag
          ! EOC - Xu
       ENDDO

!
!--    Apply prescribed horizontal wall heatflux where necessary. First, determine start and end
!--    index for respective (j,i)-index. Please note, in the flat case following loops will not be
!--    entered, as surf_s=1 and surf_e=0. Furtermore, note, no vertical natural surfaces so far.
!--    First, for default-type surfaces.
!--    North-facing vertical default-type surfaces
       surf_s = surf_def_v(0)%start_index(j,i)
       surf_e = surf_def_v(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_def_v(0)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_def_v_north(m) * ddy
       ENDDO
!
!--    South-facing vertical default-type surfaces
       surf_s = surf_def_v(1)%start_index(j,i)
       surf_e = surf_def_v(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_def_v(1)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_def_v_south(m) * ddy
       ENDDO
!
!--    East-facing vertical default-type surfaces
       surf_s = surf_def_v(2)%start_index(j,i)
       surf_e = surf_def_v(2)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_def_v(2)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_def_v_east(m) * ddx
       ENDDO
!
!--    West-facing vertical default-type surfaces
       surf_s = surf_def_v(3)%start_index(j,i)
       surf_e = surf_def_v(3)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_def_v(3)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_def_v_west(m) * ddx
       ENDDO
!
!--    Now, for natural-type surfaces
!--    North-facing
       surf_s = surf_lsm_v(0)%start_index(j,i)
       surf_e = surf_lsm_v(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_lsm_v(0)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_north(m) * ddy
       ENDDO
!
!--    South-facing
       surf_s = surf_lsm_v(1)%start_index(j,i)
       surf_e = surf_lsm_v(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_lsm_v(1)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_south(m) * ddy
       ENDDO
!
!--    East-facing
       surf_s = surf_lsm_v(2)%start_index(j,i)
       surf_e = surf_lsm_v(2)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_lsm_v(2)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_east(m) * ddx
       ENDDO
!
!--    West-facing
       surf_s = surf_lsm_v(3)%start_index(j,i)
       surf_e = surf_lsm_v(3)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_lsm_v(3)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_lsm_v_west(m) * ddx
       ENDDO
!
!--    Now, for urban-type surfaces
!--    North-facing
       surf_s = surf_usm_v(0)%start_index(j,i)
       surf_e = surf_usm_v(0)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_usm_v(0)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_north(m) * ddy
       ENDDO
!
!--    South-facing
       surf_s = surf_usm_v(1)%start_index(j,i)
       surf_e = surf_usm_v(1)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_usm_v(1)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_south(m) * ddy
       ENDDO
!
!--    East-facing
       surf_s = surf_usm_v(2)%start_index(j,i)
       surf_e = surf_usm_v(2)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_usm_v(2)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_east(m) * ddx
       ENDDO
!
!--    West-facing
       surf_s = surf_usm_v(3)%start_index(j,i)
       surf_e = surf_usm_v(3)%end_index(j,i)
       DO  m = surf_s, surf_e
          k           = surf_usm_v(3)%k(m)
          tend(k,j,i) = tend(k,j,i) + s_flux_usm_v_west(m) * ddx
       ENDDO

       ! BOC - Xu
!
!--    Compute vertical diffusion. In case that surface fluxes have been prescribed or computed at
!--    bottom and/or top, index k starts/ends at nzb+2 or nzt-1, respectively. Model top is also
!--    mask if top flux is given.
       !DO  k = nzb+1, nzt
!
!--       Determine flags to mask topography below and above. Flag 0 is used to mask topography in
!--       general, and flag 8 implies information about use_surface_fluxes. Flag 9 is used to
!--       control flux at model top.
       !   mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
       !   mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) )  *                &
       !                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )
       !   flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 0 ) )

       !   tend(k,j,i) = tend(k,j,i)                                                                &
       !                                + 0.5_wp * (                                                &
       !                               ( kh(k,j,i) + kh(k+1,j,i) ) *                                &
       !                                   ( s(k+1,j,i)-s(k,j,i) ) * ddzu(k+1)                      &
       !                                                     * rho_air_zw(k)                        &
       !                                                     * mask_top                             &
       !                             - ( kh(k,j,i) + kh(k-1,j,i) ) *                                &
       !                                   ( s(k,j,i)-s(k-1,j,i) ) * ddzu(k)                        &
       !                                                     * rho_air_zw(k-1)                      &
       !                                                     * mask_bottom                          &
       !                                           ) * ddzw(k) * drho_air(k)                        &
       !                                                       * flag
       !ENDDO
       ! EOC - Xu

       ! BOC - Xu
!
!--    Vertical diffusion at horizontal walls.
!--    TO DO: Adjust for downward facing walls and mask already in main loop
       IF ( use_surface_fluxes )  THEN
!
!--       Default-type surfaces, upward-facing
          surf_s = surf_def_h(0)%start_index(j,i)
          surf_e = surf_def_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(0)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_def_h_up(m) * ddzzw(k) * drho_air(k)               &
                                                            * dzeta_dz_s(j,i)
          ENDDO
!
!--       Default-type surfaces, downward-facing
          surf_s = surf_def_h(1)%start_index(j,i)
          surf_e = surf_def_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(1)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_def_h_down(m) * ddzzw(k) * drho_air(k)             &
                                                              * dzeta_dz_s(j,i)
          ENDDO
!
!--       Natural-type surfaces, upward-facing
          surf_s = surf_lsm_h(0)%start_index(j,i)
          surf_e = surf_lsm_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_lsm_h(0)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_lsm_h_up(m) * ddzzw(k) * drho_air(k)               &
                                                            * dzeta_dz_s(j,i)
          ENDDO
!
!--       Natural-type surfaces, downward-facing
          surf_s = surf_lsm_h(1)%start_index(j,i)
          surf_e = surf_lsm_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_lsm_h(1)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_lsm_h_down(m) * ddzzw(k) * drho_air(k)             &
                                                              * dzeta_dz_s(j,i)
          ENDDO
!
!--       Urban-type surfaces, upward-facing
          surf_s = surf_usm_h(0)%start_index(j,i)
          surf_e = surf_usm_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_usm_h(0)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_usm_h_up(m) * ddzzw(k) * drho_air(k)               &
                                                            * dzeta_dz_s(j,i)
          ENDDO
!
!--       Urban-type surfaces, upward-facing
          surf_s = surf_usm_h(1)%start_index(j,i)
          surf_e = surf_usm_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_usm_h(1)%k(m)
             tend(k,j,i) = tend(k,j,i) + s_flux_usm_h_down(m) * ddzzw(k) * drho_air(k)             &
                                                              * dzeta_dz_s(j,i)
          ENDDO
       ENDIF
!
!--    Vertical diffusion at the last computational gridpoint along z-direction
       IF ( use_top_fluxes )  THEN
          surf_s = surf_def_h(2)%start_index(j,i)
          surf_e = surf_def_h(2)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(2)%k(m)
             tend(k,j,i) = tend(k,j,i) + ( - s_flux_t(m) ) * ddzzw(k) * drho_air(k)                &
                                                           * dzeta_dz_s(j,i)
          ENDDO
       ENDIF
       ! EOC - Xu

    END SUBROUTINE diffusion_s_ij

 END MODULE diffusion_s_mod
