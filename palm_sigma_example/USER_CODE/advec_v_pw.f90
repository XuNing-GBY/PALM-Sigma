!> @file advec_v_pw.f90
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
!> Advection term for v velocity-component using Piacsek and Williams.
!> Vertical advection at the first grid point above the surface is done with normal centred
!> differences, because otherwise no information from the surface would be communicated upwards due
!> to w=0 at K=nzb.
!--------------------------------------------------------------------------------------------------!
 MODULE advec_v_pw_mod
 

    PRIVATE
    PUBLIC advec_v_pw

    INTERFACE advec_v_pw
       MODULE PROCEDURE advec_v_pw
       MODULE PROCEDURE advec_v_pw_ij
    END INTERFACE advec_v_pw
 
 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_v_pw

    USE arrays_3d,                                                                                 &
        ONLY:  ddzw, tend, u, v, w

    USE control_parameters,                                                                        &
        ONLY:  u_gtrans, v_gtrans

    USE grid_variables,                                                                            &
        ONLY:  ddx, ddy

    USE indices,                                                                                   &
        ONLY:  nxl, nxr, nyn, nysv, nzb, nzt

    USE kinds


    IMPLICIT NONE

    INTEGER(iwp) ::  i !< grid index along x-direction
    INTEGER(iwp) ::  j !< grid index along y-direction
    INTEGER(iwp) ::  k !< grid index along z-direction
    
    REAL(wp)     ::  gu !< Galilei-transformation velocity along x
    REAL(wp)     ::  gv !< Galilei-transformation velocity along y


    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans
    DO  i = nxl, nxr
       DO  j = nysv, nyn
          DO  k = nzb+1, nzt
             tend(k,j,i) = tend(k,j,i) - 0.25_wp * (                                               &
                           ( v(k,j,i+1)   * ( u(k,j-1,i+1) + u(k,j,i+1) - gu )                     &
                           - v(k,j,i-1)   * ( u(k,j-1,i)   + u(k,j,i)   - gu ) ) * ddx             &
                           + ( v(k,j+1,i) * ( v(k,j+1,i)   + v(k,j,i)   - gv )                     &
                           - v(k,j-1,i)   * ( v(k,j,i)     + v(k,j-1,i) - gv ) ) * ddy             &
                           + ( v(k+1,j,i) * ( w(k,j-1,i)   + w(k,j,i) )                            &
                           - v(k-1,j,i)   * ( w(k-1,j-1,i) + w(k-1,j,i) ) )      * ddzw(k)         &
                                                   )
          ENDDO
       ENDDO
    ENDDO

 END SUBROUTINE advec_v_pw


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE advec_v_pw_ij( i, j )

    ! BOC - Xu
    USE arrays_3d,                                                                                 &
        ONLY:  ddzzu, ddzzw, dd2zzu, tend, u, v, w
    ! EOC - Xu

    USE control_parameters,                                                                        &
        ONLY:  u_gtrans, v_gtrans

    USE grid_variables,                                                                            &
        ONLY:  ddx, ddy

    USE indices,                                                                                   &
        ONLY:  nzb, nzt

    USE kinds

    ! BOC - Xu
    USE wind_wave_inta_mod,                                                                        &
        ONLY:  dzeta_dt_v, dzeta_dx_v, dzeta_dy_v, dzeta_dz_v
    ! EOC - Xu

    IMPLICIT NONE

    INTEGER(iwp) ::  i !< grid index along x-direction
    INTEGER(iwp) ::  j !< grid index along y-direction
    INTEGER(iwp) ::  k !< grid index along z-direction
    
    REAL(wp)     ::  gu !< Galilei-transformation velocity along x
    REAL(wp)     ::  gv !< Galilei-transformation velocity along y


    gu = 2.0_wp * u_gtrans
    gv = 2.0_wp * v_gtrans
    ! BOC - Xu
    DO  k = nzb+1, nzt
       IF ( k == nzb+1 ) THEN
          tend(k,j,i) = tend(k,j,i) - dzeta_dt_v(k,j,i) * ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zzu(k)  &
                      - 0.25_wp * (                                                                &
                        ( v(k,j+1,i) * ( v(k,j+1,i) + v(k,j,i) - gv )                              &
                        - v(k,j-1,i) * ( v(k,j-1,i) + v(k,j,i) - gv ) ) * ddy                      &
                      + ( v(k+1,j,i) * ( v(k+1,j,i) + v(k,j,i) - gv )                              &
                        - v(k-1,j,i) * ( v(k,j,i) + v(k-1,j,i) - gv ) )                            &
                        * ddzzw(k) * dzeta_dy_v(k,j,i) * 1.33333333333_wp                          &
                      + ( v(k,j,i+1) * ( u(k,j,i+1) + u(k,j-1,i+1) - gu )                          &
                        - v(k,j,i-1) * ( u(k,j-1,i) + u(k,j,i)     - gu ) ) * ddx                  &
                      + ( v(k+1,j,i) * ( u(k+1,j,i) + u(k+1,j-1,i) + u(k+1,j,i+1) + u(k+1,j-1,i+1) &
                                       + u(k,j,i)   + u(k,j-1,i)   + u(k,j,i+1)   + u(k,j-1,i+1)   &
                                       - 4.0_wp * gu ) * 0.25_wp                                   &
                        - v(k-1,j,i) * ( u(k-1,j,i) + u(k-1,j-1,i) + u(k-1,j,i+1) + u(k-1,j-1,i+1) &
                                       + u(k,j,i)   + u(k,j-1,i)   + u(k,j,i+1)   + u(k,j-1,i+1)   &
                                       - 4.0_wp * gu ) * 0.25_wp )                                 &
                        * ddzzw(k) * dzeta_dx_v(k,j,i) * 1.33333333333_wp                          &
                      + ( v(k+1,j,i) * ( w(k,j,i) + w(k,j-1,i) )                                   &
                        - v(k-1,j,i) * ( 0.75_wp * ( w(k-1,j,i) + w(k-1,j-1,i) )                   &
                                       + 0.25_wp * ( w(k,j,i)   + w(k,j-1,i)   ) ) )               &
                        * ddzzw(k) * dzeta_dz_v(j,i) * 1.33333333333_wp                            &
                                  )
       ELSE
          tend(k,j,i) = tend(k,j,i) - dzeta_dt_v(k,j,i) * ( v(k+1,j,i) - v(k-1,j,i) ) * dd2zzu(k)  &
                      - 0.25_wp * (                                                                &
                        ( v(k,j+1,i) * ( v(k,j+1,i) + v(k,j,i) - gv )                              &
                        - v(k,j-1,i) * ( v(k,j-1,i) + v(k,j,i) - gv ) ) * ddy                      &
                      + ( v(k+1,j,i) * ( v(k+1,j,i) + v(k,j,i) - gv )                              &
                        - v(k-1,j,i) * ( v(k,j,i) + v(k-1,j,i) - gv ) )                            &
                        * ddzzw(k) * dzeta_dy_v(k,j,i)                                             &
                      + ( v(k,j,i+1) * ( u(k,j,i+1) + u(k,j-1,i+1) - gu )                          &
                        - v(k,j,i-1) * ( u(k,j-1,i) + u(k,j,i)     - gu ) ) * ddx                  &
                      + ( v(k+1,j,i) * ( u(k+1,j,i) + u(k+1,j-1,i) + u(k+1,j,i+1) + u(k+1,j-1,i+1) &
                                       + u(k,j,i)   + u(k,j-1,i)   + u(k,j,i+1)   + u(k,j-1,i+1)   &
                                       - 4.0_wp * gu ) * 0.25_wp                                   &
                        - v(k-1,j,i) * ( u(k-1,j,i) + u(k-1,j-1,i) + u(k-1,j,i+1) + u(k-1,j-1,i+1) &
                                       + u(k,j,i)   + u(k,j-1,i)   + u(k,j,i+1)   + u(k,j-1,i+1)   &
                                       - 4.0_wp * gu ) * 0.25_wp )                                 &
                        * ddzzw(k) * dzeta_dx_v(k,j,i)                                             &
                      + ( v(k+1,j,i) * ( w(k,j,i)   + w(k,j-1,i)   )                               &
                        - v(k-1,j,i) * ( w(k-1,j,i) + w(k-1,j-1,i) ) )                             &
                        * ddzzw(k) * dzeta_dz_v(j,i)                                               &
                                  )
       ENDIF
    ENDDO
    ! EOC - Xu

 END SUBROUTINE advec_v_pw_ij

 END MODULE advec_v_pw_mod
 
