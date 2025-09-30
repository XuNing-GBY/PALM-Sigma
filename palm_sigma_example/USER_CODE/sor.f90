!> @file sor.f90
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
!
! Description:
! ------------
!> Solve the Poisson-equation with the SOR-Red/Black-scheme.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE sor( d, ddzzu, ddzzw, p )

    USE arrays_3d,                                                                                 &
        ONLY:  rho_air,                                                                            &
               rho_air_zw

    USE control_parameters,                                                                        &
        ONLY:  bc_dirichlet_l,                                                                     &
               bc_dirichlet_n,                                                                     &
               bc_dirichlet_r,                                                                     &
               bc_dirichlet_s,                                                                     &
               bc_lr_cyc,                                                                          &
               bc_ns_cyc,                                                                          &
               bc_radiation_l,                                                                     &
               bc_radiation_n,                                                                     &
               bc_radiation_r,                                                                     &
               bc_radiation_s,                                                                     &
               ibc_p_b,                                                                            &
               ibc_p_t,                                                                            &
               n_sor,                                                                              &
               omega_sor

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz

    USE grid_variables,                                                                            &
        ONLY:  ddx2,                                                                               &
               ddy2,                                                                               &
               ddx,                                                                                &
               ddy

    USE indices,                                                                                   &
        ONLY:  nbgp,                                                                               &
               nxl,                                                                                &
               nxlg,                                                                               &
               nxr,                                                                                &
               nxrg,                                                                               &
               nyn,                                                                                &
               nyng,                                                                               &
               nys,                                                                                &
               nysg,                                                                               &
               nz,                                                                                 &
               nzb,                                                                                &
               nzt

    USE kinds

    USE wind_wave_inta_mod,                                                                        &
        ONLY:  dzeta_dx_s, dzeta_dx_u,                                                             &
               dzeta_dy_s, dzeta_dy_v,                                                             &
               dzeta_dz_s

    IMPLICIT NONE

    INTEGER(iwp) ::  i     !<
    INTEGER(iwp) ::  j     !<
    INTEGER(iwp) ::  k     !<
    INTEGER(iwp) ::  n     !<
    INTEGER(iwp) ::  nxl1  !<
    INTEGER(iwp) ::  nxl2  !<
    INTEGER(iwp) ::  nys1  !<
    INTEGER(iwp) ::  nys2  !<

    REAL(wp) ::  resi

    REAL(wp) ::  ddzzu(1:nz+1)   !<
    REAL(wp) ::  ddzzw(1:nzt+1)  !<

    REAL(wp) ::  d(nzb+1:nzt,nys:nyn,nxl:nxr)      !<
    REAL(wp) ::  p(nzb:nzt+1,nysg:nyng,nxlg:nxrg)  !<

!
!-- Limits for RED- and BLACK-part.
    IF ( MOD( nxl , 2 ) == 0 )  THEN
       nxl1 = nxl
       nxl2 = nxl + 1
    ELSE
       nxl1 = nxl + 1
       nxl2 = nxl
    ENDIF
    IF ( MOD( nys , 2 ) == 0 )  THEN
       nys1 = nys
       nys2 = nys + 1
    ELSE
       nys1 = nys + 1
       nys2 = nys
    ENDIF

    DO  n = 1, n_sor

!
!--    RED-part
       DO  i = nxl1, nxr, 2
          DO  j = nys2, nyn, 2
             k = nzb+1 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                               &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                               &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k,j,i+1)   - p(k,j,i-1)   )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k,j+1,i)   - p(k,j-1,i)   )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k) * (                       &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                                                  ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                                                  ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) *                                       &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                             &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                             &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) ) ) * p(k,j,i)                      &
                        )
             k = nzb+2 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                  - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                  - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                  - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                  - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                  ( p(k,j,i-1) + p(k-1,j,i-1) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i-1) + p(k-2,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                  ( p(k,j,i+1) + p(k-1,j,i) + p(k-1,j,i+1)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                  ( p(k,j-1,i) + p(k-1,j-1,i) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j-1,i) + p(k-2,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                  ( p(k,j+1,i) + p(k-1,j,i) + p(k-1,j+1,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                          + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                              )                                     &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          ) * p(k,j,i)                                                              &
                        )
             DO  k = nzb+3, nzt-1
                 p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                          ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          ) *                                                                         &
                          ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                          + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                    - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                            - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                    - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                            + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                    + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                       )                              &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                    - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                            - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                    - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                            + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                    + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                       )                              &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                              ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                              ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                              ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                              ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                    ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                    ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                    ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                    ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                               )                      &
                          + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                              ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                            + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                                )                                     &
                          - d(k,j,i)                                                                  &
                          - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                            + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            ) * p(k,j,i)                                                              &
                          )
             ENDDO
             k = nzt
             p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                      ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                        + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                      + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      ) *                                                                         &
                      ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                      + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                      + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                        - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                        + dzeta_dx_s(k,j,i)   * ( p(k,j,i+1)   - p(k,j,i-1)                       &
                                                + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                   )                              &
                      + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                        - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                        + dzeta_dy_s(k,j,i)   * ( p(k,j+1,i)   - p(k,j-1,i)                       &
                                                + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                   )                              &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                          ( p(k+1,j,i-1) + p(k+1,j,i) - p(k-1,j,i-1) - p(k-1,j,i) )               &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                          ( p(k+1,j,i) + p(k+1,j,i+1) - p(k-1,j,i) - p(k-1,j,i+1) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                          ( p(k+1,j-1,i) + p(k+1,j,i) - p(k-1,j-1,i) - p(k-1,j,i) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                          ( p(k+1,j,i) + p(k+1,j+1,i) - p(k-1,j,i) - p(k-1,j+1,i) )               &
                                                                           )                      &
                      - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                           )                      &
                      + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                          ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                            )                                     &
                      - d(k,j,i)                                                                  &
                      - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        ) * p(k,j,i)                                                              &
                      )
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys1, nyn, 2
             k = nzb+1 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                               &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                               &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k,j,i+1)   - p(k,j,i-1)   )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k,j+1,i)   - p(k,j-1,i)   )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k) * (                       &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                                                  ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                                                  ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) *                                       &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                             &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                             &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) ) ) * p(k,j,i)                      &
                        )
             k = nzb+2 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                  - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                  - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                  - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                  - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                  ( p(k,j,i-1) + p(k-1,j,i-1) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i-1) + p(k-2,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                  ( p(k,j,i+1) + p(k-1,j,i) + p(k-1,j,i+1)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                  ( p(k,j-1,i) + p(k-1,j-1,i) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j-1,i) + p(k-2,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                  ( p(k,j+1,i) + p(k-1,j,i) + p(k-1,j+1,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                          + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                              )                                     &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          ) * p(k,j,i)                                                              &
                        )
             DO  k = nzb+3, nzt-1
                 p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                          ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          ) *                                                                         &
                          ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                          + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                    - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                            - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                    - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                            + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                    + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                       )                              &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                    - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                            - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                    - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                            + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                    + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                       )                              &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                              ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                              ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                              ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                              ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                    ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                    ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                    ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                    ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                               )                      &
                          + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                              ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                            + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                                )                                     &
                          - d(k,j,i)                                                                  &
                          - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                            + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            ) * p(k,j,i)                                                              &
                          )
             ENDDO
             k = nzt
             p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                      ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                        + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                      + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      ) *                                                                         &
                      ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                      + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                      + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                        - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                        + dzeta_dx_s(k,j,i)   * ( p(k,j,i+1)   - p(k,j,i-1)                       &
                                                + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                   )                              &
                      + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                        - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                        + dzeta_dy_s(k,j,i)   * ( p(k,j+1,i)   - p(k,j-1,i)                       &
                                                + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                   )                              &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                          ( p(k+1,j,i-1) + p(k+1,j,i) - p(k-1,j,i-1) - p(k-1,j,i) )               &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                          ( p(k+1,j,i) + p(k+1,j,i+1) - p(k-1,j,i) - p(k-1,j,i+1) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                          ( p(k+1,j-1,i) + p(k+1,j,i) - p(k-1,j-1,i) - p(k-1,j,i) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                          ( p(k+1,j,i) + p(k+1,j+1,i) - p(k-1,j,i) - p(k-1,j+1,i) )               &
                                                                           )                      &
                      - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                           )                      &
                      + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                          ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                            )                                     &
                      - d(k,j,i)                                                                  &
                      - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        ) * p(k,j,i)                                                              &
                      )
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF

!
!--    BLACK-part
       DO  i = nxl1, nxr, 2
          DO  j = nys1, nyn, 2
             k = nzb+1 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                               &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                               &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k,j,i+1)   - p(k,j,i-1)   )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k,j+1,i)   - p(k,j-1,i)   )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k) * (                       &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                                                  ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                                                  ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) *                                       &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                             &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                             &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) ) ) * p(k,j,i)                      &
                        )
             k = nzb+2 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                  - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                  - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                  - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                  - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                  ( p(k,j,i-1) + p(k-1,j,i-1) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i-1) + p(k-2,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                  ( p(k,j,i+1) + p(k-1,j,i) + p(k-1,j,i+1)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                  ( p(k,j-1,i) + p(k-1,j-1,i) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j-1,i) + p(k-2,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                  ( p(k,j+1,i) + p(k-1,j,i) + p(k-1,j+1,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                          + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                              )                                     &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          ) * p(k,j,i)                                                              &
                        )
             DO  k = nzb+3, nzt-1
                 p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                          ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          ) *                                                                         &
                          ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                          + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                    - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                            - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                    - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                            + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                    + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                       )                              &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                    - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                            - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                    - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                            + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                    + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                       )                              &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                              ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                              ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                              ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                              ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                    ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                    ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                    ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                    ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                               )                      &
                          + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                              ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                            + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                                )                                     &
                          - d(k,j,i)                                                                  &
                          - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                            + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            ) * p(k,j,i)                                                              &
                          )
             ENDDO
             k = nzt
             p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                      ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                        + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                      + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      ) *                                                                         &
                      ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                      + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                      + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                        - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                        + dzeta_dx_s(k,j,i)   * ( p(k,j,i+1)   - p(k,j,i-1)                       &
                                                + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                   )                              &
                      + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                        - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                        + dzeta_dy_s(k,j,i)   * ( p(k,j+1,i)   - p(k,j-1,i)                       &
                                                + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                   )                              &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                          ( p(k+1,j,i-1) + p(k+1,j,i) - p(k-1,j,i-1) - p(k-1,j,i) )               &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                          ( p(k+1,j,i) + p(k+1,j,i+1) - p(k-1,j,i) - p(k-1,j,i+1) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                          ( p(k+1,j-1,i) + p(k+1,j,i) - p(k-1,j-1,i) - p(k-1,j,i) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                          ( p(k+1,j,i) + p(k+1,j+1,i) - p(k-1,j,i) - p(k-1,j+1,i) )               &
                                                                           )                      &
                      - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                           )                      &
                      + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                          ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                            )                                     &
                      - d(k,j,i)                                                                  &
                      - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        ) * p(k,j,i)                                                              &
                      )
          ENDDO
       ENDDO

       DO  i = nxl2, nxr, 2
          DO  j = nys2, nyn, 2
             k = nzb+1 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                               &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                                &
                          ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                               &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k,j,i+1)   - p(k,j,i-1)   )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k,j+1,i)   - p(k,j-1,i)   )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k) * (                       &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                                                  ( p(k+1,j,i-1) + p(k+1,j,i)   + p(k,j,i-1)        &
                                                  - 2.0_wp * ( p(k-1,j,i-1) + p(k-1,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j,i+1) + p(k,j,i+1)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                                                  ( p(k+1,j-1,i) + p(k+1,j,i)   + p(k,j-1,i)        &
                                                  - 2.0_wp * ( p(k-1,j-1,i) + p(k-1,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                                                  ( p(k+1,j,i)   + p(k+1,j+1,i) + p(k,j+1,i)        &
                                                  - 2.0_wp * ( p(k-1,j,i)   + p(k-1,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) *                                       &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dx_u(k,j,i) - dzeta_dx_u(k,j,i+1) )                             &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) *                              &
                            ( dzeta_dy_v(k,j,i) - dzeta_dy_v(k,j+1,i) )                             &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k,j,i)   + dzeta_dx_u(k,j,i+1)   )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k,j,i)   + dzeta_dy_v(k,j+1,i)   ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) ) ) * p(k,j,i)                      &
                        )
             k = nzb+2 ! rho_air level needs to be checked
             p(k,j,i) = p(k,j,i) + omega_sor /                                                      &
                        ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                        ) *                                                                         &
                        ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                        + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                        + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                  - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                          - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                  - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                          + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                  + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                     )                              &
                        + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                            dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                  - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                          - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                  - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                          + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                  + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                     )                              &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                            ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                            ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                            ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                            ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                        - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                  ( p(k,j,i-1) + p(k-1,j,i-1) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i-1) + p(k-2,j,i)   ) )      &
                          + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                  ( p(k,j,i+1) + p(k-1,j,i) + p(k-1,j,i+1)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j,i+1) ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                  ( p(k,j-1,i) + p(k-1,j-1,i) + p(k-1,j,i)          &
                                                  - 2.0_wp * ( p(k-2,j-1,i) + p(k-2,j,i)   ) )      &
                          + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                  ( p(k,j+1,i) + p(k-1,j,i) + p(k-1,j+1,i)          &
                                                  - 2.0_wp * ( p(k-2,j,i)   + p(k-2,j+1,i) ) )      &
                                                                            )                       &
                        + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                            ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                          + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                              )                                     &
                        - d(k,j,i)                                                                  &
                        - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                          ) * p(k,j,i)                                                              &
                        )
             DO  k = nzb+3, nzt-1
                 p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                          ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) )   &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                              dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                            + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                          + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                            ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                          ) *                                                                         &
                          ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                          + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                          + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                    - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                            - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                    - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                            + dzeta_dx_s(k,j,i)   * ( p(k+1,j,i+1) - p(k+1,j,i-1)                     &
                                                    + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                       )                              &
                          + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                              dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                    - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                            - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                    - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                            + dzeta_dy_s(k,j,i)   * ( p(k+1,j+1,i) - p(k+1,j-1,i)                     &
                                                    + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                       )                              &
                          + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i) *                               &
                              ( p(k+2,j,i-1) + p(k+2,j,i) - p(k,j,i-1) )                              &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k+1,j,i+1) *                             &
                              ( p(k+2,j,i) + p(k+2,j,i+1) - p(k,j,i+1) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j,i) *                               &
                              ( p(k+2,j-1,i) + p(k+2,j,i) - p(k,j-1,i) )                              &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k+1,j+1,i) *                             &
                              ( p(k+2,j,i) + p(k+2,j+1,i) - p(k,j+1,i) )       )                      &
                          - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                    ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                            + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                    ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                    ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                            + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                    ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                               )                      &
                          + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                              ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                            + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                                )                                     &
                          - d(k,j,i)                                                                  &
                          - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k+1) * rho_air(k+1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k+1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k+1,j,i) + dzeta_dy_v(k+1,j+1,i) ) ) &
                            + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                                dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                              + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                            + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                              ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                            ) * p(k,j,i)                                                              &
                          )
             ENDDO
             k = nzt
             p(k,j,i) = p(k,j,i) + omega_sor /                                                    &
                      ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                     &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                          dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )     &
                        + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) )   &
                      + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                   &
                        ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                     &
                      ) *                                                                         &
                      ( ddx2 * rho_air(k) * ( p(k,j,i+1) + p(k,j,i-1) )                           &
                      + ddy2 * rho_air(k) * ( p(k,j+1,i) + p(k,j-1,i) )                           &
                      + 0.25_wp * ddx * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dx_u(k,j,i+1) * ( p(k+1,j,i)   + p(k+1,j,i+1)                     &
                                                - p(k-1,j,i)   - p(k-1,j,i+1) )                   &
                        - dzeta_dx_u(k,j,i)   * ( p(k+1,j,i-1) + p(k+1,j,i)                       &
                                                - p(k-1,j,i-1) - p(k-1,j,i)   )                   &
                        + dzeta_dx_s(k,j,i)   * ( p(k,j,i+1)   - p(k,j,i-1)                       &
                                                + p(k-1,j,i-1) - p(k-1,j,i+1) )                   &
                                                                   )                              &
                      + 0.25_wp * ddy * ddzzw(k) * rho_air_zw(k) * (                              &
                          dzeta_dy_v(k,j+1,i) * ( p(k+1,j,i)   + p(k+1,j+1,i)                     &
                                                - p(k-1,j,i)   - p(k-1,j+1,i) )                   &
                        - dzeta_dy_v(k,j,i)   * ( p(k+1,j-1,i) + p(k+1,j,i)                       &
                                                - p(k-1,j-1,i) - p(k-1,j,i)   )                   &
                        + dzeta_dy_s(k,j,i)   * ( p(k,j+1,i)   - p(k,j-1,i)                       &
                                                + p(k-1,j-1,i) - p(k-1,j+1,i) )                   &
                                                                   )                              &
                      + 0.0625_wp * ddzzw(k) * ddzzw(k) * rho_air_zw(k)  * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i) *                                 &
                          ( p(k+1,j,i-1) + p(k+1,j,i) - p(k-1,j,i-1) - p(k-1,j,i) )               &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k,j,i+1) *                               &
                          ( p(k+1,j,i) + p(k+1,j,i+1) - p(k-1,j,i) - p(k-1,j,i+1) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j,i) *                                 &
                          ( p(k+1,j-1,i) + p(k+1,j,i) - p(k-1,j-1,i) - p(k-1,j,i) )               &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k,j+1,i) *                               &
                          ( p(k+1,j,i) + p(k+1,j+1,i) - p(k-1,j,i) - p(k-1,j+1,i) )               &
                                                                           )                      &
                      - 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i) *                               &
                                                ( p(k,j,i-1) - p(k-2,j,i-1) - p(k-2,j,i)   )      &
                        + dzeta_dx_s(k,j,i) * dzeta_dx_u(k-1,j,i+1) *                             &
                                                ( p(k,j,i+1) - p(k-2,j,i) - p(k-2,j,i+1)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j,i) *                               &
                                                ( p(k,j-1,i) - p(k-2,j-1,i) - p(k-2,j,i)   )      &
                        + dzeta_dy_s(k,j,i) * dzeta_dy_v(k-1,j+1,i) *                             &
                                                ( p(k,j+1,i) - p(k-2,j,i) - p(k-2,j+1,i)   )      &
                                                                           )                      &
                      + dzeta_dz_s(j,i) * dzeta_dz_s(j,i) * (                                     &
                          ddzzw(k) * ddzzu(k+1) * rho_air_zw(k)   * p(k+1,j,i)                    &
                        + ddzzw(k) * ddzzu(k)   * rho_air_zw(k-1) * p(k-1,j,i)                    &
                                                            )                                     &
                      - d(k,j,i)                                                                  &
                      - ( 2.0_wp * ( ddx2 + ddy2 ) * rho_air(k)                                   &
                        + 0.0625_wp * ddzzw(k) * ddzzw(k-1) * rho_air(k-1) * (                    &
                            dzeta_dx_s(k,j,i) * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k-1,j,i+1) )   &
                          + dzeta_dy_s(k,j,i) * ( dzeta_dy_v(k-1,j,i) + dzeta_dy_v(k-1,j+1,i) ) ) &
                        + ddzzw(k) * ddzzu(k+1) * rho_air_zw(k) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        + ddzzw(k) * ddzzu(k) * rho_air_zw(k-1) *                                 &
                          ( dzeta_dz_s(j,i) * dzeta_dz_s(j,i) )                                   &
                        ) * p(k,j,i)                                                              &
                      )
          ENDDO
       ENDDO

!
!--    Exchange of boundary values for p.
       CALL exchange_horiz( p, nbgp )

!
!--    Boundary conditions top/bottom.
!--    Bottom boundary
       IF ( ibc_p_b == 1 )  THEN       ! Neumann
          p(nzb,:,:) = p(nzb+1,:,:)
       ELSE                            ! Dirichlet
          p(nzb,:,:) = 0.0_wp
       ENDIF

!
!--    Top boundary
       IF ( ibc_p_t == 1 )  THEN       ! Neumann
          p(nzt+1,:,:) = p(nzt,:,:)
       ELSE                            ! Dirichlet
          p(nzt+1,:,:) = 0.0_wp
       ENDIF

!
!--    Horizontal (Neumann) boundary conditions in case of non-cyclic boundaries
       IF ( .NOT. bc_lr_cyc )  THEN
          IF ( bc_dirichlet_l  .OR.  bc_radiation_l )  p(:,:,nxl-1) = p(:,:,nxl)
          IF ( bc_dirichlet_r  .OR.  bc_radiation_r )  p(:,:,nxr+1) = p(:,:,nxr)
       ENDIF
       IF ( .NOT. bc_ns_cyc )  THEN
          IF ( bc_dirichlet_n  .OR.  bc_radiation_n )  p(:,nyn+1,:) = p(:,nyn,:)
          IF ( bc_dirichlet_s  .OR.  bc_radiation_s )  p(:,nys-1,:) = p(:,nys,:)
       ENDIF


    ENDDO

 END SUBROUTINE sor
