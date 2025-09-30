!> @file diffusion_v.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of PALM-Sigma.
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
!> Diffusion term of the v-component
!--------------------------------------------------------------------------------------------------!
 MODULE diffusion_v_mod


    PRIVATE
    PUBLIC diffusion_v

    INTERFACE diffusion_v
       MODULE PROCEDURE diffusion_v
       MODULE PROCEDURE diffusion_v_ij
    END INTERFACE diffusion_v

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v

       USE arrays_3d,                                                                              &
           ONLY:  ddzu, ddzw, drho_air, km, rho_air_zw, tend, u, v, w 

       USE control_parameters,                                                                     &
           ONLY:  constant_top_momentumflux, use_surface_fluxes,  use_top_fluxes

       USE grid_variables,                                                                         &
           ONLY:  ddx, ddy, ddy2

       USE indices,                                                                                &
           ONLY:  nxl, nxr, nyn, nysv, nzb, nzt, topo_flags

       USE kinds

       USE surface_mod,                                                                            &
           ONLY :  surf_def_h, surf_def_v, surf_lsm_h, surf_lsm_v, surf_usm_h, surf_usm_v

       IMPLICIT NONE

       INTEGER(iwp) ::  i             !< running index x direction
       INTEGER(iwp) ::  j             !< running index y direction
       INTEGER(iwp) ::  k             !< running index z direction
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  flag          !< flag to mask topography grid points
       REAL(wp)     ::  kmxm          !< diffusion coefficient on leftward side of the v-gridbox - interpolated onto xu-yv grid
       REAL(wp)     ::  kmxp          !< diffusion coefficient on rightward side of the v-gridbox - interpolated onto xu-yv grid
       REAL(wp)     ::  kmzm          !< diffusion coefficient on bottom of the gridbox - interpolated onto yv-zw grid
       REAL(wp)     ::  kmzp          !< diffusion coefficient on top of the gridbox - interpolated onto yv-zw grid
       REAL(wp)     ::  mask_bottom   !< flag to mask vertical upward-facing surface
       REAL(wp)     ::  mask_east     !< flag to mask vertical surface south of the grid point
       REAL(wp)     ::  mask_west     !< flag to mask vertical surface north of the grid point
       REAL(wp)     ::  mask_top      !< flag to mask vertical downward-facing surface

       !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(i, j, k, l, m) &
       !$ACC PRIVATE(surf_e, surf_s, flag, kmxm, kmxp, kmzm, kmzp) &
       !$ACC PRIVATE(mask_bottom, mask_east, mask_west, mask_top) &
       !$ACC PRESENT(topo_flags, km) &
       !$ACC PRESENT(u, v, w) &
       !$ACC PRESENT(ddzu, ddzw, drho_air, rho_air_zw) &
       !$ACC PRESENT(surf_def_h(0:2), surf_def_v(2:3)) &
       !$ACC PRESENT(surf_lsm_h(0:1), surf_lsm_v(2:3)) &
       !$ACC PRESENT(surf_usm_h(0:1), surf_usm_v(2:3)) &
       !$ACC PRESENT(tend)
       DO  i = nxl, nxr
          DO  j = nysv, nyn
!
!--          Compute horizontal diffusion
             DO  k = nzb+1, nzt

!
!--             Predetermine flag to mask topography and wall-bounded grid points.
!--             It is sufficient to masked only east- and west-facing surfaces, which need special
!--             treatment for the v-component.
                flag      = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i),   2 ) )
                mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 2 ) )
                mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 2 ) )
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
                kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )

                tend(k,j,i) = tend(k,j,i) +          (                                             &
                                mask_east * kmxp * (                                               &
                                       ( v(k,j,i+1) - v(k,j,i)     ) * ddx                         &
                                     + ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy                         &
                                                   )                                               &
                              - mask_west * kmxm * (                                               &
                                       ( v(k,j,i) - v(k,j,i-1) ) * ddx                             &
                                     + ( u(k,j,i) - u(k,j-1,i) ) * ddy                             &
                                                   )                                               &
                                                     ) * ddx  * flag                               &
                                          + 2.0_wp * (                                             &
                                        km(k,j,i)   * ( v(k,j+1,i) - v(k,j,i)   )                  &
                                      - km(k,j-1,i) * ( v(k,j,i)   - v(k,j-1,i) )                  &
                                                     ) * ddy2 * flag

             ENDDO

!
!--          Add horizontal momentum flux v'u' at east- (l=2) and west-facing (l=3) surfaces. Note,
!--          in the the flat case, loops won't be entered as start_index > end_index. Furtermore,
!--          note, no vertical natural surfaces so far.
!--          Default-type surfaces
             DO  l = 2, 3
                surf_s = surf_def_v(l)%start_index(j,i)
                surf_e = surf_def_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_def_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_uv(m) * ddx
                ENDDO
             ENDDO
!
!--          Natural-type surfaces
             DO  l = 2, 3
                surf_s = surf_lsm_v(l)%start_index(j,i)
                surf_e = surf_lsm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_lsm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_uv(m) * ddx
                ENDDO
             ENDDO
!
!--          Urban-type surfaces
             DO  l = 2, 3
                surf_s = surf_usm_v(l)%start_index(j,i)
                surf_e = surf_usm_v(l)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k           = surf_usm_v(l)%k(m)
                   tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_uv(m) * ddx
                ENDDO
             ENDDO
!
!--          Compute vertical diffusion. In case of simulating a surface layer, respective grid
!--          diffusive fluxes are masked (flag 10) within this loop, and added further below, else,
!--          simple gradient approach is applied. Model top is also mask if top-momentum flux is
!--          given.
             DO  k = nzb+1, nzt
!
!--             Determine flags to mask topography below and above. Flag 2 is used to mask
!--             topography in general, while flag 8 implies also information about
!--             use_surface_fluxes. Flag 9 is used to control momentum flux at model top.
                mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
                mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) ) *           &
                              MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )
                flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
!
!--             Interpolate eddy diffusivities on staggered gridpoints
                kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
                kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

                tend(k,j,i) = tend(k,j,i)                                                          &
                              & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)                 &
                              &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy                       &
                              &            ) * rho_air_zw(k)   * mask_top                          &
                              &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)               &
                              &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy                   &
                              &            ) * rho_air_zw(k-1) * mask_bottom                       &
                              &   ) * ddzw(k) * drho_air(k) * flag
             ENDDO

!
!--          Vertical diffusion at the first grid point above the surface, if the momentum flux at
!--          the bottom is given by the Prandtl law or if it is prescribed by the user.
!--          Difference quotient of the momentum flux is not formed over half of the grid spacing
!--          (2.0*ddzw(k)) any more, since the comparison with other (LES) models showed that the
!--          values of the momentum flux becomes too large in this case.
             IF ( use_surface_fluxes )  THEN
!
!--             Default-type surfaces, upward-facing
                surf_s = surf_def_h(0)%start_index(j,i)
                surf_e = surf_def_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(0)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - ( - surf_def_h(0)%vsws(m) ) ) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Default-type surfaces, dowward-facing
                surf_s = surf_def_h(1)%start_index(j,i)
                surf_e = surf_def_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_def_h(1)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - surf_def_h(1)%vsws(m) ) * ddzw(k) * drho_air(k)
                ENDDO
!
!--             Natural-type surfaces, upward-facing
                surf_s = surf_lsm_h(0)%start_index(j,i)
                surf_e = surf_lsm_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_lsm_h(0)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - ( - surf_lsm_h(0)%vsws(m) ) ) * ddzw(k) * drho_air(k)

                ENDDO
!
!--             Natural-type surfaces, downward-facing
                surf_s = surf_lsm_h(1)%start_index(j,i)
                surf_e = surf_lsm_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_lsm_h(1)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 +  ( - surf_lsm_h(1)%vsws(m) ) * ddzw(k) * drho_air(k)

                ENDDO
!
!--             Urban-type surfaces, upward-facing
                surf_s = surf_usm_h(0)%start_index(j,i)
                surf_e = surf_usm_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_usm_h(0)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - ( - surf_usm_h(0)%vsws(m) ) ) * ddzw(k) * drho_air(k)

                ENDDO
!
!--             Urban-type surfaces, downward-facing
                surf_s = surf_usm_h(1)%start_index(j,i)
                surf_e = surf_usm_h(1)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k   = surf_usm_h(1)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - surf_usm_h(1)%vsws(m) ) * ddzw(k) * drho_air(k)

                ENDDO
             ENDIF
!
!--          Add momentum flux at model top
             IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
                surf_s = surf_def_h(2)%start_index(j,i)
                surf_e = surf_def_h(2)%end_index(j,i)
                DO  m = surf_s, surf_e

                   k   = surf_def_h(2)%k(m)

                   tend(k,j,i) = tend(k,j,i)                                                       &
                                 + ( - surf_def_h(2)%vsws(m) ) * ddzw(k) * drho_air(k)
                ENDDO
             ENDIF

          ENDDO
       ENDDO

    END SUBROUTINE diffusion_v


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
    SUBROUTINE diffusion_v_ij( i, j )

       ! BOC - Xu
       USE arrays_3d,                                                                              &
           ONLY:  ddzu, ddzw, drho_air, km, tend, u, v, w, rho_air_zw, ddzzu, ddzzw
       ! EOC - Xu

       USE control_parameters,                                                                     &
           ONLY:  constant_top_momentumflux, use_surface_fluxes, use_top_fluxes

       USE grid_variables,                                                                         &
           ONLY:  ddx, ddy, ddy2

       USE indices,                                                                                &
           ONLY:  nzb, nzt, topo_flags

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
       INTEGER(iwp) ::  l             !< running index of surface type, south- or north-facing wall
       INTEGER(iwp) ::  m             !< running index surface elements
       INTEGER(iwp) ::  surf_e        !< End index of surface elements at (j,i)-gridpoint
       INTEGER(iwp) ::  surf_s        !< Start index of surface elements at (j,i)-gridpoint

       REAL(wp)     ::  flag          !< flag to mask topography grid points
       REAL(wp)     ::  kmxm          !< diffusion coefficient on leftward side of the v-gridbox - interpolated onto xu-yv grid
       REAL(wp)     ::  kmxp          !< diffusion coefficient on rightward side of the v-gridbox - interpolated onto xu-yv grid
       REAL(wp)     ::  kmzm          !< diffusion coefficient on bottom of the gridbox - interpolated onto xu-zw grid
       REAL(wp)     ::  kmzp          !< diffusion coefficient on top of the gridbox - interpolated onto xu-zw grid
       REAL(wp)     ::  mask_bottom   !< flag to mask vertical upward-facing surface
       REAL(wp)     ::  mask_east     !< flag to mask vertical surface south of the grid point
       REAL(wp)     ::  mask_west     !< flag to mask vertical surface north of the grid point
       REAL(wp)     ::  mask_top      !< flag to mask vertical downward-facing surface

       ! BOC - Xu
!
!--    Compute horizontal diffusion
       DO  k = nzb+1, nzt
!
!--       Predetermine flag to mask topography and wall-bounded grid points.
!--       It is sufficient to masked only east- and west-facing surfaces, which need special
!--       treatment for the v-component.
          flag      = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i),   2 ) )
          mask_east = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i+1), 2 ) )
          mask_west = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i-1), 2 ) )
          mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
          mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) ) *         &
                        MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
          kmxp = 0.25_wp * ( km(k,j,i)+km(k,j,i+1)+km(k,j-1,i)+km(k,j-1,i+1) )
          kmxm = 0.25_wp * ( km(k,j,i)+km(k,j,i-1)+km(k,j-1,i)+km(k,j-1,i-1) )
          kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
          kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

          tend(k,j,i) = tend(k,j,i) +                                                              &
                      + ( kmxp * ( ( u(k,j,i+1) - u(k,j-1,i+1) ) * ddy                             &
                                 + ( ( u(k,j,i+1) + u(k,j-1,i+1) + u(k+1,j,i+1) + u(k+1,j-1,i+1) ) &
                                   - ( u(k,j,i+1) + u(k,j-1,i+1) + u(k-1,j,i+1) + u(k-1,j-1,i+1) ) &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dy_v(k,j,i) + dzeta_dy_v(k,j,i+1) ) * 0.5_wp          &
                                 + ( v(k,j,i+1) - v(k,j,i) ) * ddx                                 &
                                 + ( ( v(k+1,j,i+1) + v(k+1,j,i) + v(k,j,i+1) + v(k,j,i) )         &
                                   - ( v(k-1,j,i+1) + v(k-1,j,i) + v(k,j,i+1) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dx_v(k,j,i) + dzeta_dx_v(k,j,i+1) ) * 0.5_wp          &
                                 ) * mask_east                                                     &
                        - kmxm * ( ( u(k,j,i) - u(k,j-1,i) ) * ddy                                 &
                                 + ( ( u(k,j,i) + u(k,j-1,i) + u(k+1,j,i) + u(k+1,j-1,i) )         &
                                   - ( u(k,j,i) + u(k,j-1,i) + u(k-1,j,i) + u(k-1,j-1,i) )         &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dy_v(k,j,i) + dzeta_dy_v(k,j,i-1) ) * 0.5_wp          &
                                 + ( v(k,j,i) - v(k,j,i-1) ) * ddx                                 &
                                 + ( ( v(k+1,j,i-1) + v(k+1,j,i) + v(k,j,i-1) + v(k,j,i) )         &
                                   - ( v(k-1,j,i-1) + v(k-1,j,i) + v(k,j,i-1) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dx_v(k,j,i) + dzeta_dx_v(k,j,i-1) ) * 0.5_wp          &
                                 ) * mask_west                                                     &
                        ) * ddx * flag                                                             &
                      + ( kmzp * ( ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1) + u(k,j,i) )         &
                                   - ( u(k+1,j-1,i+1) + u(k+1,j-1,i) + u(k,j-1,i+1) + u(k,j-1,i) ) &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k+1,j-1,i+1) + u(k+1,j-1,i) ) &
                                   - ( u(k,j,i+1) + u(k,j,i) + u(k,j-1,i+1) + u(k,j-1,i) )         &
                                   ) * 0.25_wp * ddzzu(k+1)                                        &
                                   * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 + ( ( v(k+1,j,i+1) + v(k+1,j,i) + v(k,j,i+1) + v(k,j,i) )         &
                                   - ( v(k+1,j,i-1) + v(k+1,j,i) + v(k,j,i-1) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( v(k+1,j,i) - v(k,j,i) ) * ddzzu(k+1)                          &
                                   * ( dzeta_dx_v(k,j,i) + dzeta_dx_v(k+1,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k) * mask_top                                      &
                        - kmzm * ( ( ( u(k-1,j,i+1) + u(k-1,j,i) + u(k,j,i+1) + u(k,j,i) )         &
                                   - ( u(k-1,j-1,i+1) + u(k-1,j-1,i) + u(k,j-1,i+1) + u(k,j-1,i) ) &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( ( u(k,j,i+1) + u(k,j,i) + u(k,j-1,i+1) + u(k,j-1,i) )         &
                                   - ( u(k-1,j,i+1) + u(k-1,j,i) + u(k-1,j-1,i+1) + u(k-1,j-1,i) ) &
                                   ) * 0.25_wp * ddzzu(k)                                          &
                                   * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 + ( ( v(k-1,j,i+1) + v(k-1,j,i) + v(k,j,i+1) + v(k,j,i) )         &
                                   - ( v(k-1,j,i-1) + v(k-1,j,i) + v(k,j,i-1) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( v(k,j,i) - v(k-1,j,i) ) * ddzzu(k)                            &
                                   * ( dzeta_dx_v(k,j,i) + dzeta_dx_v(k-1,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1) * mask_bottom                                 &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dx_v(k,j,i) * flag                      &
                      + ( km(k,j,i)   * ( ( v(k,j+1,i) - v(k,j,i) ) * ddy                          &
                                        + ( ( v(k+1,j+1,i) + v(k+1,j,i) + v(k,j+1,i) + v(k,j,i) )  &
                                          - ( v(k-1,j+1,i) + v(k-1,j,i) + v(k,j+1,i) + v(k,j,i) )  &
                                          ) * 0.25_wp * ddzzw(k)                                   &
                                          * ( dzeta_dy_v(k,j+1,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp   &
                                        )                                                          &
                        - km(k,j-1,i) * ( ( v(k,j,i) - v(k,j-1,i) ) * ddy                          &
                                        + ( ( v(k+1,j-1,i) + v(k+1,j,i) + v(k,j-1,i) + v(k,j,i) )  &
                                          - ( v(k-1,j-1,i) + v(k-1,j,i) + v(k,j-1,i) + v(k,j,i) )  &
                                          ) * 0.25_wp * ddzzw(k)                                   &
                                          * ( dzeta_dy_v(k,j-1,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp   &
                                        )                                                          &
                        ) * 2.0_wp * ddy * flag                                                    &
                      + ( kmzp * ( ( ( v(k+1,j+1,i) + v(k+1,j,i) + v(k,j+1,i) + v(k,j,i) )         &
                                   - ( v(k+1,j-1,i) + v(k+1,j,i) + v(k,j-1,i) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( v(k+1,j,i) - v(k,j,i) ) * ddzzu(k+1)                          &
                                   * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k) * mask_top                                      &
                        - kmzm * ( ( ( v(k-1,j+1,i) + v(k-1,j,i) + v(k,j+1,i) + v(k,j,i) )         &
                                   - ( v(k-1,j-1,i) + v(k-1,j,i) + v(k,j-1,i) + v(k,j,i) )         &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( v(k,j,i) - v(k-1,j,i) ) * ddzzu(k)                            &
                                   * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1) * mask_bottom                                 &
                        ) * 2.0_wp * ddzzw(k) * drho_air(k) * dzeta_dy_v(k,j,i) * flag             &
                      + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzzu(k+1)                          &
                                   * ( dzeta_dz_v(j,i) + dzeta_dz_v(j,i) ) * 0.5_wp                &
                                 + ( w(k,j,i) - w(k,j-1,i) ) * ddy                                 &
                                 + ( ( w(k+1,j,i) + w(k+1,j-1,i) + w(k,j,i) + w(k,j-1,i) )         &
                                   - ( w(k-1,j,i) + w(k-1,j-1,i) + w(k,j,i) + w(k,j-1,i) )         &
                                   ) * 0.25_wp * ddzzu(k+1)                                        &
                                   * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k) * mask_top                                      &
                        - kmzm * ( ( v(k,j,i) - v(k-1,j,i) ) * ddzzu(k)                            &
                                   * ( dzeta_dz_v(j,i) + dzeta_dz_v(j,i) ) * 0.5_wp                &
                                 + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy                             &
                                 + ( ( w(k,j,i) + w(k,j-1,i) + w(k-1,j,i) + w(k-1,j-1,i) )         &
                                   - ( w(k-1,j,i) + w(k-1,j-1,i) + w(k-2,j,i) + w(k-2,j-1,i) )     &
                                   ) * 0.25_wp * ddzzu(k)                                          &
                                   * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1) * mask_bottom                                 &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dz_v(j,i) * flag
       ENDDO

!
!--    Add horizontal momentum flux v'u' at east- (l=2) and west-facing (l=3) surfaces. Note, in the
!--    the flat case, loops won't be entered as start_index > end_index. Furtermore, note, no
!--    vertical natural surfaces so far.
!--    Default-type surfaces
       DO  l = 2, 3
          surf_s = surf_def_v(l)%start_index(j,i)
          surf_e = surf_def_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_def_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_def_v(l)%mom_flux_uv(m) * ddx
          ENDDO
       ENDDO
!
!--    Natural-type surfaces
       DO  l = 2, 3
          surf_s = surf_lsm_v(l)%start_index(j,i)
          surf_e = surf_lsm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_lsm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_lsm_v(l)%mom_flux_uv(m) * ddx
          ENDDO
       ENDDO
!
!--    Urban-type surfaces
       DO  l = 2, 3
          surf_s = surf_usm_v(l)%start_index(j,i)
          surf_e = surf_usm_v(l)%end_index(j,i)
          DO  m = surf_s, surf_e
             k           = surf_usm_v(l)%k(m)
             tend(k,j,i) = tend(k,j,i) + surf_usm_v(l)%mom_flux_uv(m) * ddx
          ENDDO
       ENDDO
!
!--    Compute vertical diffusion. In case of simulating a surface layer, respective grid diffusive
!--    fluxes are masked (flag 8) within this loop, and added further below, else, simple gradient
!--    approach is applied. Model top is also mask if top-momentum flux is given.
       !DO  k = nzb+1, nzt
!
!--       Determine flags to mask topography below and above. Flag 2 is used to mask topography in
!--       general, while flag 10 implies also information about use_surface_fluxes. Flag 9 is used
!--       to control momentum flux at model top.
       !   mask_bottom = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )
       !   mask_top    = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 8 ) ) *                 &
       !                 MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k+1,j,i), 9 ) )
       !   flag        = MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 2 ) )
!
!--       Interpolate eddy diffusivities on staggered gridpoints
       !   kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j-1,i)+km(k+1,j-1,i) )
       !   kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j-1,i)+km(k-1,j-1,i) )

       !   tend(k,j,i) = tend(k,j,i)                                                                &
       !                 & + ( kmzp * ( ( v(k+1,j,i) - v(k,j,i) ) * ddzu(k+1)                       &
       !                 &            + ( w(k,j,i) - w(k,j-1,i) ) * ddy                             &
       !                 &            ) * rho_air_zw(k)   * mask_top                                &
       !                 &   - kmzm * ( ( v(k,j,i)   - v(k-1,j,i)   ) * ddzu(k)                     &
       !                 &            + ( w(k-1,j,i) - w(k-1,j-1,i) ) * ddy                         &
       !                 &            ) * rho_air_zw(k-1) * mask_bottom                             &
       !                 &   ) * ddzw(k) * drho_air(k) * flag
       !ENDDO

!
!--    Vertical diffusion at the first grid point above the surface, if the momentum flux at the
!--    bottom is given by the Prandtl law or if it is prescribed by the user.
!--    Difference quotient of the momentum flux is not formed over half of the grid spacing
!--    (2.0*ddzw(k)) any more, since the comparison with other (LES) models showed that the values
!--    of the momentum flux becomes too large in this case.
       IF ( use_surface_fluxes )  THEN
!
!--       Default-type surfaces, upward-facing
          surf_s = surf_def_h(0)%start_index(j,i)
          surf_e = surf_def_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(0)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - ( - surf_def_h(0)%vsws(m) ) ) * ddzzw(k) * drho_air(k) &
                                                                           * dzeta_dz_v(j,i)
          ENDDO
!
!--       Default-type surfaces, dowward-facing
          surf_s = surf_def_h(1)%start_index(j,i)
          surf_e = surf_def_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_def_h(1)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - surf_def_h(1)%vsws(m) ) * ddzzw(k) * drho_air(k)       &
                                                                     * dzeta_dz_v(j,i)
          ENDDO
!
!--       Natural-type surfaces, upward-facing
          surf_s = surf_lsm_h(0)%start_index(j,i)
          surf_e = surf_lsm_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_lsm_h(0)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - ( - surf_lsm_h(0)%vsws(m) ) ) * ddzzw(k) * drho_air(k) &
                                                                           * dzeta_dz_v(j,i)
          ENDDO
!
!--       Natural-type surfaces, downward-facing
          surf_s = surf_lsm_h(1)%start_index(j,i)
          surf_e = surf_lsm_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_lsm_h(1)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - surf_lsm_h(1)%vsws(m) ) * ddzzw(k) * drho_air(k)       &
                                                                     * dzeta_dz_v(j,i)
          ENDDO
!
!--       Urban-type surfaces, upward-facing
          surf_s = surf_usm_h(0)%start_index(j,i)
          surf_e = surf_usm_h(0)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_usm_h(0)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - ( - surf_usm_h(0)%vsws(m) ) ) * ddzzw(k) * drho_air(k) &
                                                                           * dzeta_dz_v(j,i)
          ENDDO
!
!--       Urban-type surfaces, downward-facing
          surf_s = surf_usm_h(1)%start_index(j,i)
          surf_e = surf_usm_h(1)%end_index(j,i)
          DO  m = surf_s, surf_e
             k   = surf_usm_h(1)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - surf_usm_h(1)%vsws(m) ) * ddzzw(k) * drho_air(k)       &
                                                                     * dzeta_dz_v(j,i)
          ENDDO
       ENDIF
!
!--    Add momentum flux at model top
       IF ( use_top_fluxes  .AND.  constant_top_momentumflux )  THEN
          surf_s = surf_def_h(2)%start_index(j,i)
          surf_e = surf_def_h(2)%end_index(j,i)
          DO  m = surf_s, surf_e

             k   = surf_def_h(2)%k(m)

             tend(k,j,i) = tend(k,j,i) + ( - surf_def_h(2)%vsws(m) ) * ddzzw(k) * drho_air(k)       &
                                                                     * dzeta_dz_v(j,i)
          ENDDO
       ENDIF
       ! EOC - Xu

    END SUBROUTINE diffusion_v_ij

 END MODULE diffusion_v_mod
