!> @file wind_wave_inta_mod.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of PALM-Sigma.
!
! PALM-Sigma consists of a series of user-defined modules and modified PALM modules, developed for
! wave-phase-resolved LES. PALM-Sigma should be used as USER_CODE based on the framework of PALM.
!
! PALM-Sigma is free software: you can redistribute it and/or modify it under the terms of the GNU
! General Public License as published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! PALM-Sigma is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! Copyright 2025 Xu Ning
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> A user-defined module for implementation of the sigma coordinate in PALM.
!--------------------------------------------------------------------------------------------------!


MODULE wind_wave_inta_mod
    USE basic_constants_and_equations_mod,                                                         &
        ONLY:  pi

    USE arrays_3d,                                                                                 &
        ONLY:  zu, zzu,                                                                            &
               zw, zzw,                                                                            &
               y,                                                                                  &
               xu,                                                                                 &
               yv,                                                                                 &
               x,                                                                                  &
               ddzzw,                                                                              &
               zzu_mg, zzw_mg,                                                                     &
               u_p, v_p, w_p

    USE control_parameters,                                                                        &
        ONLY:  simulated_time,                                                                     &
               intermediate_timestep_count,                                                        &
               dt_3d,                                                                              &
               initializing_actions,                                                               &
               psolver,                                                                            &
               maximum_grid_level,                                                                 &
               wind_wave_inta

    USE indices,                                                                                   &
        ONLY:  nxl, nxlg,                                                                          &
               nxr, nxrg,                                                                          &
               nys, nysg,                                                                          &
               nyn, nyng,                                                                          &
               nzb,                                                                                &
               nzt,                                                                                &
               nz,                                                                                 &
               nbgp,                                                                               &
               nxl_mg, nxr_mg, nys_mg, nyn_mg, nzt_mg

    USE kinds

    USE grid_variables,                                                                            &
        ONLY:  ddx, ddx2, ddx_mg, ddx2_mg,                                                         &
               ddy, ddy2, ddy_mg, ddy2_mg,                                                         &
               dx,                                                                                 &
               dy,                                                                                 &
               H_av

    USE netcdf_data_input_mod,                                                                     &
        ONLY:  input_pids_static,                                                                  &
               char_fill,                                                                          &
               check_existence,                                                                    &
               close_input_file,                                                                   &
               get_attribute,                                                                      &
               get_dimension_length,                                                               &
               get_variable,                                                                       &
               inquire_num_variables,                                                              &
               inquire_variable_names,                                                             &
               input_file_static,                                                                  &
               num_var_pids,                                                                       &
               open_read_file,                                                                     &
               pids_id,                                                                            &
               real_3d,                                                                            &
               vars_pids

    USE surface_mod,                                                                               &
        ONLY:  surf_def_h,                                                                         &
               surf_type

    USE interp,                                                                                    &
        ONLY:  interp3d 

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz,                                                                     &
               exchange_horiz_2d

    IMPLICIT NONE

    CHARACTER (LEN=16) :: wwinta_mode = 'rgh' ! or 'wbl', or 'sfc'
    CHARACTER (LEN=16) :: wave_type = 'regular'
    REAL(wp) :: wave_amplitude = 0.0_wp, wave_number = 0.0_wp, wave_direction = 0.0_wp

    REAL(wp) :: u_rws, v_rws, w_rws  

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: t_, y_, x_ 

    ! read eta and phi from netcdf file and derive u, v, w, cosx, cosy of wave surface from eta and phi
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: eta_field, phi_field
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: u_field, v_field, w_field
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: cosx_field, cosy_field, sinx_field, siny_field

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: eta_u, eta_v, eta_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: spv_u, spv_v, spv_w ! surface particle velocity, only defined in grid center.
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: sia_sin_x, sia_cos_x, sia_sin_y, sia_cos_y    ! surface inclination (angle) sin and cos along x and y direction
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: deta_dt_u, deta_dt_v, deta_dt_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: deta_dx_u, deta_dx_v, deta_dx_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: deta_dy_u, deta_dy_v, deta_dy_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: deta_dxx_u, deta_dxx_v, deta_dxx_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: deta_dyy_u, deta_dyy_v, deta_dyy_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dzeta_dt_s, dzeta_dt_u, dzeta_dt_v, dzeta_dt_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dzeta_dx_s, dzeta_dx_u, dzeta_dx_v, dzeta_dx_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dzeta_dy_s, dzeta_dy_u, dzeta_dy_v, dzeta_dy_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dzeta_dxx_s, dzeta_dxx_u, dzeta_dxx_v, dzeta_dxx_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: dzeta_dyy_s, dzeta_dyy_u, dzeta_dyy_v, dzeta_dyy_w
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: dzeta_dz_u, dzeta_dz_v, dzeta_dz_s
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: alpha, beta

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: eta_u_mg, eta_v_mg, eta_s_mg
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: deta_dx_mg, deta_dy_mg, deta_dxx_mg, deta_dyy_mg, dzeta_dz_mg
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: dzeta_dx_mg, dzeta_dy_mg, dzeta_dxx_mg, dzeta_dyy_mg, alpha_mg, beta_mg

    INTERFACE wwinta_exchange_horiz
       MODULE PROCEDURE wwinta_exchange_horiz
    END INTERFACE wwinta_exchange_horiz

    INTERFACE wwinta_prognostic_equations
       MODULE PROCEDURE wwinta_prognostic_equations_ij
    END INTERFACE wwinta_prognostic_equations

    INTERFACE wwinta_actions
       MODULE PROCEDURE wwinta_actions
       MODULE PROCEDURE wwinta_actions_ij
    END INTERFACE wwinta_actions

 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &wind_wave_inta_par for wind wave interaction model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wwinta_parin

    CHARACTER(LEN=100) ::  line  !< dummy string that contains the current line of the parameter file

    INTEGER(iwp) ::  io_status   !< status after reading the namelist file

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module
                                             !namelist appears in
                                             !< the namelist file

    NAMELIST /wind_wave_inta_parameters/                                    &
                               wwinta_mode,                                 &
                               wave_type,                                   &
                               wave_amplitude,                              &
                               wave_number,                                 &
                               wave_direction                       

!
!-- read (find) wind-wave interaction module namelist:
    REWIND ( 11 )
    READ( 11, wind_wave_inta_parameters, IOSTAT=io_status )

!
!-- Action depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    wind_wave_inta_parameters namelist was found and read correctly. Enable the
!--    wind wave interaction module.
       IF ( .NOT. switch_off_module )  wind_wave_inta = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    wind_wave_inta_parameters namelist was found but contained errors. Print an
!error message
!--    including the line that caused the problem.
       BACKSPACE( 11 )
       READ( 11 , '(A)' ) line
       CALL parin_fail_message( 'wind_wave_inta_parameters', line )

    ENDIF

 END SUBROUTINE wwinta_parin 

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> This routine writes the respective restart data.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wwinta_wrd_local


    !IMPLICIT NONE


    !IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

    !   IF ( ALLOCATED( deta_dt_u ) )  THEN
    !      CALL wrd_write_string( 'deta_dt_u' )
    !      WRITE ( 14 )  deta_dt_u
    !   ENDIF

    !ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

       !CALL wrd_mpi_io_global_array( 'deta_dt_u', deta_dt_u )

    !ENDIF

 END SUBROUTINE wwinta_wrd_local

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 !SUBROUTINE wwinta_rrd_global_ftn( found )


    !USE control_parameters,                                               &
    !    ONLY: length, restart_string


    !IMPLICIT NONE

    !LOGICAL, INTENT(OUT) ::  found


    !found = .TRUE.


    !SELECT CASE ( restart_string(1:length) )

    !   CASE ( 'deta_dt_u' )
    !      READ ( 13 )  deta_dt_u

    !   CASE DEFAULT
    !      found = .FALSE.

    !END SELECT

 !END SUBROUTINE wwinta_rrd_global_ftn

!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!------------------------------------------------------------------------------!
 SUBROUTINE wwinta_rrd_global_mpi

    !CALL rrd_mpi_io_global_array( 'deta_dt_u', deta_dt_u )

 END SUBROUTINE wwinta_rrd_global_mpi

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check &wind_wave_inta_par control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wwinta_check_parameters

 END SUBROUTINE wwinta_check_parameters

!
!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate wind wave interaction model arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wwinta_init_arrays

 END SUBROUTINE wwinta_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialization of the wind wave interaction model
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE wwinta_init

    INTEGER(iwp) ::  i, j

    IF ( TRIM( wwinta_mode ) == 'sfc' ) THEN

       IF (ALLOCATED(eta_u)) DEALLOCATE(eta_u); ALLOCATE( eta_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(eta_v)) DEALLOCATE(eta_v); ALLOCATE( eta_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(eta_s)) DEALLOCATE(eta_s); ALLOCATE( eta_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(spv_u)) DEALLOCATE(spv_u); ALLOCATE( spv_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(spv_v)) DEALLOCATE(spv_v); ALLOCATE( spv_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(spv_w)) DEALLOCATE(spv_w); ALLOCATE( spv_w(nysg:nyng,nxlg:nxrg) )

       IF ( TRIM( initializing_actions ) /= 'read_restart_data' )  THEN

          DO i = nxl, nxr
             !eta_u(:,i) = wa * SIN( - wn * xu(i) )
             !eta_v(:,i) = wa * SIN( - wn * x(i)  )
             !eta_s(:,i) = wa * SIN( - wn * x(i)  )
          ENDDO

       ELSEIF ( TRIM( wave_type ) == 'regular' ) THEN
          DO j = nys, nyn
             DO i = nxl, nxr
                eta_u(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * simulated_time                &
                                                 + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)     &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )              
                eta_v(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * simulated_time                &
                                                 + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + yv(j)*COS(wave_direction*pi/180.0_wp) ) )             
                eta_s(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * simulated_time                &
                                                 + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )

                spv_u(j,i) = - SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                               SIN( SQRT(9.81_wp*wave_number) * simulated_time                               &
                                  + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                                  +  y(j)*COS(wave_direction*pi/180.0_wp) ) ) *              &
                               SIN(wave_direction*pi/180.0_wp)
                spv_v(j,i) = - SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                               SIN( SQRT(9.81_wp*wave_number) * simulated_time                               &
                                  + wave_number * (  x(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                                  + yv(j)*COS(wave_direction*pi/180.0_wp) ) ) *              &
                               COS(wave_direction*pi/180.0_wp)
                spv_w(j,i) =   SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                               COS( SQRT(9.81_wp*wave_number) * simulated_time                               &
                                  + wave_number * (  x(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                                  +  y(j)*COS(wave_direction*pi/180.0_wp) ) )
             ENDDO
          ENDDO

       ELSEIF ( TRIM( wave_type ) == 'steady' ) THEN
          DO j = nys, nyn
             DO i = nxl, nxr
                eta_u(j,i) = wave_amplitude * SIN( wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)     &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )
                eta_v(j,i) = wave_amplitude * SIN( wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + yv(j)*COS(wave_direction*pi/180.0_wp) ) )
                eta_s(j,i) = wave_amplitude * SIN( wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )

                spv_u(j,i) = 0.0_wp
                spv_v(j,i) = 0.0_wp
                spv_w(j,i) = 0.0_wp
             ENDDO
          ENDDO

       ENDIF

       CALL exchange_horiz_2d( eta_u )
       CALL exchange_horiz_2d( eta_v )
       CALL exchange_horiz_2d( eta_s )

       CALL exchange_horiz_2d( spv_u )
       CALL exchange_horiz_2d( spv_v )
       CALL exchange_horiz_2d( spv_w )

       ! nys:nyn,nxl,nxr is enough ???
       IF (ALLOCATED(sia_sin_x)) DEALLOCATE(sia_sin_x); ALLOCATE( sia_sin_x(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(sia_cos_x)) DEALLOCATE(sia_cos_x); ALLOCATE( sia_cos_x(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(sia_sin_y)) DEALLOCATE(sia_sin_y); ALLOCATE( sia_sin_y(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(sia_cos_y)) DEALLOCATE(sia_cos_y); ALLOCATE( sia_cos_y(nysg:nyng,nxlg:nxrg) )

       ! nys:nyn,nxl,nxr is enough ???
       ! nzb+1:nzt,nys:nyn,nxl:nxr is enough ???
       IF (ALLOCATED(deta_dt_u)) DEALLOCATE(deta_dt_u); ALLOCATE( deta_dt_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dt_v)) DEALLOCATE(deta_dt_v); ALLOCATE( deta_dt_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dt_s)) DEALLOCATE(deta_dt_s); ALLOCATE( deta_dt_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dx_u)) DEALLOCATE(deta_dx_u); ALLOCATE( deta_dx_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dx_v)) DEALLOCATE(deta_dx_v); ALLOCATE( deta_dx_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dx_s)) DEALLOCATE(deta_dx_s); ALLOCATE( deta_dx_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dy_u)) DEALLOCATE(deta_dy_u); ALLOCATE( deta_dy_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dy_v)) DEALLOCATE(deta_dy_v); ALLOCATE( deta_dy_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dy_s)) DEALLOCATE(deta_dy_s); ALLOCATE( deta_dy_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dxx_u)) DEALLOCATE(deta_dxx_u); ALLOCATE( deta_dxx_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dxx_v)) DEALLOCATE(deta_dxx_v); ALLOCATE( deta_dxx_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dxx_s)) DEALLOCATE(deta_dxx_s); ALLOCATE( deta_dxx_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dyy_u)) DEALLOCATE(deta_dyy_u); ALLOCATE( deta_dyy_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dyy_v)) DEALLOCATE(deta_dyy_v); ALLOCATE( deta_dyy_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(deta_dyy_s)) DEALLOCATE(deta_dyy_s); ALLOCATE( deta_dyy_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dt_s)) DEALLOCATE(dzeta_dt_s); ALLOCATE( dzeta_dt_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dt_u)) DEALLOCATE(dzeta_dt_u); ALLOCATE( dzeta_dt_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dt_v)) DEALLOCATE(dzeta_dt_v); ALLOCATE( dzeta_dt_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dt_w)) DEALLOCATE(dzeta_dt_w); ALLOCATE( dzeta_dt_w(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dx_s)) DEALLOCATE(dzeta_dx_s); ALLOCATE( dzeta_dx_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dx_u)) DEALLOCATE(dzeta_dx_u); ALLOCATE( dzeta_dx_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dx_v)) DEALLOCATE(dzeta_dx_v); ALLOCATE( dzeta_dx_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dx_w)) DEALLOCATE(dzeta_dx_w); ALLOCATE( dzeta_dx_w(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dy_s)) DEALLOCATE(dzeta_dy_s); ALLOCATE( dzeta_dy_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dy_u)) DEALLOCATE(dzeta_dy_u); ALLOCATE( dzeta_dy_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dy_v)) DEALLOCATE(dzeta_dy_v); ALLOCATE( dzeta_dy_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dy_w)) DEALLOCATE(dzeta_dy_w); ALLOCATE( dzeta_dy_w(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dz_u)) DEALLOCATE(dzeta_dz_u); ALLOCATE( dzeta_dz_u(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dz_v)) DEALLOCATE(dzeta_dz_v); ALLOCATE( dzeta_dz_v(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dz_s)) DEALLOCATE(dzeta_dz_s); ALLOCATE( dzeta_dz_s(nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dxx_s)) DEALLOCATE(dzeta_dxx_s); ALLOCATE( dzeta_dxx_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dxx_u)) DEALLOCATE(dzeta_dxx_u); ALLOCATE( dzeta_dxx_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dxx_v)) DEALLOCATE(dzeta_dxx_v); ALLOCATE( dzeta_dxx_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dxx_w)) DEALLOCATE(dzeta_dxx_w); ALLOCATE( dzeta_dxx_w(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dyy_s)) DEALLOCATE(dzeta_dyy_s); ALLOCATE( dzeta_dyy_s(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dyy_u)) DEALLOCATE(dzeta_dyy_u); ALLOCATE( dzeta_dyy_u(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dyy_v)) DEALLOCATE(dzeta_dyy_v); ALLOCATE( dzeta_dyy_v(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(dzeta_dyy_w)) DEALLOCATE(dzeta_dyy_w); ALLOCATE( dzeta_dyy_w(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(alpha)) DEALLOCATE(alpha); ALLOCATE( alpha(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       IF (ALLOCATED(beta))  DEALLOCATE(beta);  ALLOCATE(  beta(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )

       CALL calc_metrics
    ENDIF

 END SUBROUTINE wwinta_init

 SUBROUTINE read_wave_field

     CHARACTER(LEN=128) :: filename

     INTEGER(iwp) :: t_dim, y_dim, x_dim
     INTEGER(iwp) :: t, j, i                 ! running index

     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: deta

     REAL(wp) :: dx, dy

     filename = "/cluster/work/users/xni001/JOBS/wwinta_gs20_regular/USER_CODE/waveData_regular.nc"

     ! open netcdf file that contains eta
     CALL open_read_file(filename, pids_id)

     ! get number of variables
     CALL inquire_num_variables( pids_id, num_var_pids )

     ! Allocate memory to store variable name, i.e. eta
     ALLOCATE( vars_pids(1:num_var_pids) )
     CALL inquire_variable_names( pids_id, vars_pids )

     ! get lengths of dimensions
     CALL get_dimension_length( pids_id, t_dim, 'time')
     CALL get_dimension_length( pids_id, y_dim, 'y')
     CALL get_dimension_length( pids_id, x_dim, 'x')

     ! get variable values of time, y, x
     IF (allocated(x_)) DEALLOCATE(x_); ALLOCATE(x_(x_dim))
     !CALL get_variable( pids_id, 'x', x_, 0, x_dim-1)
     CALL get_variable( pids_id, 'x', x_ )
     IF (allocated(y_)) DEALLOCATE(y_); ALLOCATE(y_(y_dim))
     !CALL get_variable( pids_id, 'y', y_, 0, y_dim-1)
     CALL get_variable( pids_id, 'y', y_)
     IF (allocated(t_)) DEALLOCATE(t_); ALLOCATE(t_(t_dim))
     !CALL get_variable( pids_id, 'time', t_, 0, t_dim-1)
     CALL get_variable( pids_id, 'time', t_)

     ! get variable values of eta
     IF (ALLOCATED(eta_field)) DEALLOCATE(eta_field); ALLOCATE( eta_field(1:t_dim,1:y_dim,1:x_dim) )
     CALL get_variable( pids_id, 'eta', eta_field, 0, x_dim-1, 0, y_dim-1, 0, t_dim-1)

     ! get variable values of phi
     IF (ALLOCATED(phi_field)) DEALLOCATE(phi_field); ALLOCATE( phi_field(1:t_dim,1:y_dim,1:x_dim) )
     CALL get_variable( pids_id, 'phi', phi_field, 0, x_dim-1, 0, y_dim-1, 0, t_dim-1)

     ! derive u,v,w,cosx,cosy from eta and phi
     IF (ALLOCATED(u_field)) DEALLOCATE(u_field); ALLOCATE(u_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(v_field)) DEALLOCATE(v_field); ALLOCATE(v_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(w_field)) DEALLOCATE(w_field); ALLOCATE(w_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(cosx_field)) DEALLOCATE(cosx_field); ALLOCATE(cosx_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(cosy_field)) DEALLOCATE(cosy_field); ALLOCATE(cosy_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(sinx_field)) DEALLOCATE(sinx_field); ALLOCATE(sinx_field(1:t_dim,1:y_dim,1:x_dim))
     IF (ALLOCATED(siny_field)) DEALLOCATE(siny_field); ALLOCATE(siny_field(1:t_dim,1:y_dim,1:x_dim))

     DO i = 1,x_dim
        IF (i == 1) THEN
           u_field(:,:,i) = (phi_field(:,:,i+1) - phi_field(:,:,i)) / (x_(i+1) - x_(i))
           deta = eta_field(:,:,i+1) - eta_field(:,:,i)
           dx = x_(i+1) - x_(i)
        ELSEIF (i == x_dim) THEN
           u_field(:,:,i) = (phi_field(:,:,i) - phi_field(:,:,i-1)) / (x_(i) - x_(i-1))
           deta = eta_field(:,:,i) - eta_field(:,:,i-1)
           dx = x_(i) - x_(i-1)
        ELSE
           u_field(:,:,i) = (phi_field(:,:,i+1) - phi_field(:,:,i-1)) / (x_(i+1) - x_(i-1))
           deta = eta_field(:,:,i+1) - eta_field(:,:,i-1)
           dx = x_(i+1) - x_(i-1)
        ENDIF
        cosx_field(:,:,i) = dx / SQRT(dx**2 + deta**2)
        sinx_field(:,:,i) = deta / SQRT(dx**2 + deta**2)  
     ENDDO

     DO j = 1,y_dim
        IF (j == 1) THEN
           v_field(:,j,:) = (phi_field(:,j+1,:) - phi_field(:,j,:)) / (y_(j+1) - y_(j))
           deta = eta_field(:,j+1,:) - eta_field(:,j,:)
           dy = y_(j+1) - y_(j)
        ELSEIF (j == y_dim) THEN
           v_field(:,j,:) = (phi_field(:,j,:) - phi_field(:,j-1,:)) / (y_(j) - y_(j-1))
           deta = eta_field(:,j,:) - eta_field(:,j-1,:)
           dy = y_(j) - y_(j-1)
        ELSE
           v_field(:,j,:) = (phi_field(:,j+1,:) - phi_field(:,j-1,:)) / (y_(j+1) - y_(j-1))
           deta = eta_field(:,j+1,:) - eta_field(:,j-1,:)
           dy = y_(j+1) - y_(j-1)
        ENDIF
        cosy_field(:,j,:) = dy / SQRT(dy**2 + deta**2)
        siny_field(:,j,:) = deta / SQRT(dy**2 + deta**2)
     ENDDO

     DO t = 1,t_dim
        IF (t == 1) THEN
           w_field(t,:,:) = (eta_field(t+1,:,:) - eta_field(t,:,:)) / (t_(t+1) - t_(t))
        ELSEIF (t == t_dim) THEN
           w_field(t,:,:) = (eta_field(t,:,:) - eta_field(t-1,:,:)) / (t_(t) - t_(t-1))
        ELSE
           w_field(t,:,:) = (eta_field(t+1,:,:) - eta_field(t-1,:,:)) / (t_(t+1) - t_(t-1))
        ENDIF
     ENDDO

 END SUBROUTINE read_wave_field

 SUBROUTINE change_roughness

    TYPE(surf_type), POINTER ::  surf

    INTEGER(iwp) :: m, i, j, k
    REAL(wp) :: eta

    surf => surf_def_h(0)

    DO m = 1, surf%ns
       i   = surf%i(m)
       j   = surf%j(m)
       k   = surf%k(m)

       !eta = interp3d(j,i)   
       surf%z0(m) = 0.1_wp
    ENDDO

 END SUBROUTINE change_roughness

 SUBROUTINE calc_u_rws(t,y,x,u,w)

    IMPLICIT NONE

    !INTEGER(iwp), INTENT(IN) ::  i,j     ! running index
    REAL(wp), INTENT(IN) :: t,y,x,u,w
    !REAL(wp), INTENT(OUT) :: u_rws       ! u-component relative velocity to wave surface
    REAL(wp) :: u_s, w_s, cosx, sinx    ! u-component and v-component velocity of wave surface, wave slope angle at x-direction

    u_s = interp3d(t_, y_, x_, u_field, t, y, x)
    w_s = interp3d(t_, y_, x_, w_field, t, y, x)
    cosx = interp3d(t_, y_, x_, cosx_field, t, y, x)
    sinx = interp3d(t_, y_, x_, sinx_field, t, y, x)

    u_rws = (u - u_s) * cosx + ( w - w_s ) * sinx

 END SUBROUTINE calc_u_rws

 SUBROUTINE calc_v_rws(t,y,x,v,w)

    IMPLICIT NONE

    !INTEGER(iwp), INTENT(IN) ::  i,j,k     ! running index
    REAL(wp), INTENT(IN) :: t,y,x,v,w
    !REAL(wp), INTENT(OUT) :: v_rws ! u-component relative velocity to wave surface
    REAL(wp) :: v_s, w_s, cosy, siny ! u-component and v-component velocity of wave surface, wave slope angle at x-direction
    
    !v_rws = ( v(k,j,i) - v_s ) * COS(cosy) + ( w(k,j,i) - w_s ) * SIN(cosy)
    v_s = interp3d(t_, y_, x_, v_field, t, y, x)
    w_s = interp3d(t_, y_, x_, w_field, t, y, x)
    cosy = interp3d(t_, y_, x_, cosy_field, t, y, x)
    siny = interp3d(t_, y_, x_, siny_field, t, y, x)

    v_rws = (v - v_s) * cosy + ( w - w_s ) * siny 

 END SUBROUTINE calc_v_rws

 
 SUBROUTINE wwinta_prognostic_equations_ij( i, j )

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j

    REAL(wp) :: dt_itm

    IF ( intermediate_timestep_count == 1 ) THEN
       dt_itm = 1.0_wp / 3.0_wp * dt_3d 
    ELSEIF ( intermediate_timestep_count == 2 ) THEN
       dt_itm = 3.0_wp / 4.0_wp * dt_3d
    ELSE
       dt_itm = dt_3d
    ENDIF

    IF ( TRIM( wave_type ) == 'regular' ) THEN
       eta_u(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                        + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)     &
                                                        + y(j)*COS(wave_direction*pi/180.0_wp) ) )        
       eta_v(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                        + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                        + yv(j)*COS(wave_direction*pi/180.0_wp) ) )        
       eta_s(j,i) = wave_amplitude * SIN( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                        + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                        + y(j)*COS(wave_direction*pi/180.0_wp) ) )
       spv_u(j,i) = - SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                      SIN( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                    &
                         + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                         +  y(j)*COS(wave_direction*pi/180.0_wp) ) ) *              &
                      SIN(wave_direction*pi/180.0_wp)
       spv_v(j,i) = - SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                      SIN( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                    &
                         + wave_number * (  x(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                         + yv(j)*COS(wave_direction*pi/180.0_wp) ) ) *              &
                               COS(wave_direction*pi/180.0_wp)
       spv_w(j,i) =   SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                      COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                    &
                         + wave_number * (  x(i)*SIN(wave_direction*pi/180.0_wp)                    &
                                         +  y(j)*COS(wave_direction*pi/180.0_wp) ) )

    ELSEIF ( TRIM( wave_type ) == 'steady' ) THEN
          eta_u(j,i) = wave_amplitude * SIN( wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)     &
                                                           + y(j)*COS(wave_direction*pi/180.0_wp) ) )
          eta_v(j,i) = wave_amplitude * SIN( wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                           + yv(j)*COS(wave_direction*pi/180.0_wp) ) )
          eta_s(j,i) = wave_amplitude * SIN( wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                           + y(j)*COS(wave_direction*pi/180.0_wp) ) )

          spv_u(j,i) = 0.0_wp
          spv_v(j,i) = 0.0_wp
          spv_w(j,i) = 0.0_wp

    ENDIF

 END SUBROUTINE wwinta_prognostic_equations_ij


 SUBROUTINE wwinta_actions( location )
 
    CHARACTER(LEN=*) ::  location

    SELECT CASE ( location )

       CASE ( 'after_prognostic_equations' )

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE wwinta_actions


 SUBROUTINE wwinta_actions_ij( i, j, location )

    CHARACTER(LEN=*) ::  location

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j         
         
    SELECT CASE ( location )
    
       CASE ( 'after_prognostic_equations' )
          CALL calc_metrics_ij( i, j )

       CASE DEFAULT
          CONTINUE

    END SELECT 

 END SUBROUTINE wwinta_actions_ij


 SUBROUTINE calc_metrics

    INTEGER(iwp) ::  i, j, k, l
    REAL(wp)     ::  deta_x, deta_y, dis_x, dis_y

    REAL(wp) :: dt_itm

    IF ( intermediate_timestep_count == 1 ) THEN
       dt_itm = 1.0_wp / 3.0_wp * dt_3d
    ELSEIF ( intermediate_timestep_count == 2 ) THEN
       dt_itm = 3.0_wp / 4.0_wp * dt_3d
    ELSE
       dt_itm = dt_3d
    ENDIF

    DO j = nys, nyn
       DO i = nxl, nxr
          deta_x    = eta_u(j,i+1) - eta_u(j,i)
          dis_x     = SQRT( deta_x**2 + dx**2 )
          sia_sin_x(j,i) = deta_x / dis_x
          sia_cos_x(j,i) = dx / dis_x

          deta_y    = eta_v(j+1,i) - eta_v(j,i)
          dis_y     = SQRT( deta_y**2 + dy**2 )
          sia_sin_y(j,i) = deta_y / dis_y
          sia_cos_y(j,i) = dy / dis_y
       ENDDO
    ENDDO

    CALL exchange_horiz_2d( sia_sin_x )
    CALL exchange_horiz_2d( sia_cos_x )
    CALL exchange_horiz_2d( sia_sin_y )
    CALL exchange_horiz_2d( sia_cos_y )


    IF ( TRIM( wave_type ) == 'regular' ) THEN
       DO j = nys, nyn
          DO i = nxl, nxr
             deta_dt_u(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                                              COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                                 + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)     &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )
             deta_dt_v(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                                              COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                                 + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + yv(j)*COS(wave_direction*pi/180.0_wp) ) ) 
             deta_dt_s(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                  &
                                              COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)     &
                                                 + wave_number * ( x(i)*SIN(wave_direction*pi/180.0_wp)      &
                                                                 + y(j)*COS(wave_direction*pi/180.0_wp) ) )
          ENDDO
       ENDDO

    ELSEIF ( TRIM( wave_type ) == 'steady' ) THEN
       DO j = nys, nyn
          DO i = nxl, nxr
             deta_dt_u(j,i)  = 0.0_wp
             deta_dt_v(j,i)  = 0.0_wp
             deta_dt_s(j,i)  = 0.0_wp
          ENDDO
       ENDDO

    ENDIF

    CALL exchange_horiz_2d( deta_dt_u )
    CALL exchange_horiz_2d( deta_dt_v )
    CALL exchange_horiz_2d( deta_dt_s )

    DO j = nys, nyn
       DO i = nxl, nxr
          deta_dx_u(j,i)  = ( eta_s(j,i)  - eta_s(j,i-1)  ) * ddx
          deta_dx_v(j,i)  = ( eta_u(j,i+1) + eta_u(j-1,i+1) -                                                &
                              eta_u(j,i)   - eta_u(j-1,i)   ) * ddx * 0.5_wp
          deta_dx_s(j,i)  = ( eta_u(j,i+1) - eta_u(j,i)     ) * ddx
          deta_dy_u(j,i)  = ( eta_v(j+1,i) + eta_v(j+1,i-1) -                                                &
                              eta_v(j,i)   - eta_v(j,i-1)   ) * ddy * 0.5_wp
          deta_dy_v(j,i)  = ( eta_s(j,i)  - eta_s(j-1,i)  ) * ddy
          deta_dy_s(j,i)  = ( eta_v(j+1,i) - eta_v(j,i)     ) * ddy

          dzeta_dz_u(j,i) = 1.0_wp / ( H_av - eta_u(j,i)  )
          dzeta_dz_v(j,i) = 1.0_wp / ( H_av - eta_v(j,i)  )
          dzeta_dz_s(j,i) = 1.0_wp / ( H_av - eta_s(j,i) )
       ENDDO
    ENDDO

    CALL exchange_horiz_2d( deta_dx_u )
    CALL exchange_horiz_2d( deta_dx_v )
    CALL exchange_horiz_2d( deta_dx_s )
    CALL exchange_horiz_2d( deta_dy_u )
    CALL exchange_horiz_2d( deta_dy_v )
    CALL exchange_horiz_2d( deta_dy_s )

    CALL exchange_horiz_2d( dzeta_dz_u )
    CALL exchange_horiz_2d( dzeta_dz_v )
    CALL exchange_horiz_2d( dzeta_dz_s )

    DO k = nzb, nzt+1  ! nzt+1 or nzt ???
       DO j = nys, nyn
          DO i = nxl, nxr
             dzeta_dt_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dt_s(j,i)
             dzeta_dt_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dt_u(j,i)
             dzeta_dt_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dt_v(j,i)
             dzeta_dt_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dt_s(j,i)
             dzeta_dx_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dx_s(j,i)
             dzeta_dx_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dx_u(j,i)
             dzeta_dx_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dx_v(j,i)
             dzeta_dx_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dx_s(j,i)
             dzeta_dy_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dy_s(j,i)
             dzeta_dy_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dy_u(j,i)
             dzeta_dy_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dy_v(j,i)
             dzeta_dy_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dy_s(j,i)
          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz( dzeta_dt_s, nbgp )
    CALL exchange_horiz( dzeta_dt_u, nbgp )
    CALL exchange_horiz( dzeta_dt_v, nbgp )
    CALL exchange_horiz( dzeta_dt_w, nbgp )
    CALL exchange_horiz( dzeta_dx_s, nbgp )
    CALL exchange_horiz( dzeta_dx_u, nbgp )
    CALL exchange_horiz( dzeta_dx_v, nbgp )
    CALL exchange_horiz( dzeta_dx_w, nbgp )
    CALL exchange_horiz( dzeta_dy_s, nbgp )
    CALL exchange_horiz( dzeta_dy_u, nbgp )
    CALL exchange_horiz( dzeta_dy_v, nbgp )
    CALL exchange_horiz( dzeta_dy_w, nbgp )

    DO j = nys, nyn
       DO i = nxl, nxr
         deta_dxx_u(j,i)  = ( deta_dx_s(j,i)  - deta_dx_s(j,i-1) ) * ddx
         deta_dxx_v(j,i)  = ( deta_dx_u(j,i+1) + deta_dx_u(j-1,i+1)                                       &
                            - deta_dx_u(j,i)   - deta_dx_u(j-1,i)  ) * ddx * 0.5_wp
         deta_dxx_s(j,i)  = ( deta_dx_u(j,i+1) - deta_dx_u(j,i)    ) * ddx
         deta_dyy_u(j,i)  = ( deta_dy_v(j+1,i) + deta_dy_v(j+1,i-1)                                       &
                            - deta_dy_v(j,i)   - deta_dy_v(j,i-1)  ) * ddy * 0.5_wp
         deta_dyy_v(j,i)  = ( deta_dy_s(j,i)  - deta_dy_s(j-1,i) ) * ddy
         deta_dyy_s(j,i)  = ( deta_dy_v(j+1,i) - deta_dy_v(j,i)    ) * ddy
       ENDDO
    ENDDO

    CALL exchange_horiz_2d( deta_dxx_u )
    CALL exchange_horiz_2d( deta_dxx_v )
    CALL exchange_horiz_2d( deta_dxx_s )
    CALL exchange_horiz_2d( deta_dyy_u )
    CALL exchange_horiz_2d( deta_dyy_v )
    CALL exchange_horiz_2d( deta_dyy_s )

    DO k = nzb, nzt+1  ! ???
       DO j = nys, nyn
          DO i = nxl, nxr
             dzeta_dxx_s(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp           &
                                          * deta_dx_s(j,i)**2.0_wp +                                      &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                   &
                                          * deta_dxx_s(j,i)
             dzeta_dxx_u(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )**2.0_wp           &
                                          * deta_dx_u(j,i)**2.0_wp  +                                     &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )                   &
                                          * deta_dxx_u(j,i)
             dzeta_dxx_v(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )**2.0_wp           &
                                          * deta_dx_v(j,i)**2.0_wp  +                                     &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )                   &
                                          * deta_dxx_v(j,i)
             dzeta_dxx_w(k,j,i)  = 2.0_wp * ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp           &
                                          * deta_dx_s(j,i)**2.0_wp +                                      &
                                            ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                   &
                                          * deta_dxx_s(j,i)
             dzeta_dyy_s(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp           &
                                          * deta_dy_s(j,i)**2.0_wp +                                      &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                   &
                                          * deta_dyy_s(j,i)
             dzeta_dyy_u(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )**2.0_wp           &
                                          * deta_dy_u(j,i)**2.0_wp  +                                     &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )                   &
                                          * deta_dyy_u(j,i)
             dzeta_dyy_v(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )**2.0_wp           &
                                          * deta_dy_v(j,i)**2.0_wp  +                                     &
                                            ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )                   &
                                          * deta_dyy_v(j,i)
             dzeta_dyy_w(k,j,i)  = 2.0_wp * ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp           &
                                          * deta_dy_s(j,i)**2.0_wp +                                      &
                                            ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                   &
                                          * deta_dyy_s(j,i)
          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz( dzeta_dxx_s, nbgp )
    CALL exchange_horiz( dzeta_dxx_u, nbgp )
    CALL exchange_horiz( dzeta_dxx_v, nbgp )
    CALL exchange_horiz( dzeta_dxx_w, nbgp )
    CALL exchange_horiz( dzeta_dyy_s, nbgp )
    CALL exchange_horiz( dzeta_dyy_u, nbgp )
    CALL exchange_horiz( dzeta_dyy_v, nbgp )
    CALL exchange_horiz( dzeta_dyy_w, nbgp )

    DO k = nzb, nzt+1
       DO j = nys, nyn
          DO i = nxl, nxr
             alpha(k,j,i) = dzeta_dx_s(k,j,i)**2.0_wp + dzeta_dy_s(k,j,i)**2.0_wp + dzeta_dz_s(j,i)**2.0_wp
             beta(k,j,i)  = dzeta_dxx_s(k,j,i) + dzeta_dyy_s(k,j,i)
          ENDDO
       ENDDO
    ENDDO

    CALL exchange_horiz( alpha, nbgp )
    CALL exchange_horiz( beta,  nbgp )

 END SUBROUTINE calc_metrics




 SUBROUTINE calc_metrics_ij( i, j )

    INTEGER(iwp) ::  i
    INTEGER(iwp) ::  j
    INTEGER(iwp) ::  k

    REAL(wp)     ::  dis_x, dis_y, deta_x, deta_y

    REAL(wp) :: dt_itm

    IF ( intermediate_timestep_count == 1 ) THEN
       dt_itm = 1.0_wp / 3.0_wp * dt_3d
    ELSEIF ( intermediate_timestep_count == 2 ) THEN
       dt_itm = 3.0_wp / 4.0_wp * dt_3d
    ELSE
       dt_itm = dt_3d
    ENDIF


    deta_x    = eta_u(j,i+1) - eta_u(j,i)
    dis_x     = SQRT( deta_x**2 + dx**2 )
    sia_sin_x(j,i) = deta_x / dis_x
    sia_cos_x(j,i) = dx / dis_x

    deta_y    = eta_v(j+1,i) - eta_v(j,i)
    dis_y     = SQRT( deta_y**2 + dy**2 )
    sia_sin_y(j,i) = deta_y / dis_y
    sia_cos_y(j,i) = dy / dis_y

    IF ( TRIM( wave_type ) == 'regular' ) THEN
       deta_dt_u(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                     &
                         COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                       &
                            + wave_number * ( xu(i)*SIN(wave_direction*pi/180.0_wp)                       &
                                            + y(j) *COS(wave_direction*pi/180.0_wp) ) )
       deta_dt_v(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                     &
                         COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                       &
                            + wave_number * ( x(i) *SIN(wave_direction*pi/180.0_wp)                       &
                                            + yv(j)*COS(wave_direction*pi/180.0_wp) ) )
       deta_dt_s(j,i)  = SQRT(9.81_wp*wave_number) * wave_amplitude *                                     &
                         COS( SQRT(9.81_wp*wave_number) * (simulated_time + dt_itm)                       &
                            + wave_number * ( x(i) *SIN(wave_direction*pi/180.0_wp)                       &
                                            + y(j) *COS(wave_direction*pi/180.0_wp) ) )

    ELSEIF ( TRIM( wave_type ) == 'steady' ) THEN
       deta_dt_u(j,i)  = 0.0_wp
       deta_dt_v(j,i)  = 0.0_wp
       deta_dt_s(j,i)  = 0.0_wp

    ENDIF

    deta_dx_u(j,i)  = ( eta_s(j,i)   - eta_s(j,i-1) ) * ddx
    deta_dx_v(j,i)  = ( eta_u(j,i+1) + eta_u(j-1,i+1) - eta_u(j,i) - eta_u(j-1,i) ) * ddx * 0.5_wp
    deta_dx_s(j,i)  = ( eta_u(j,i+1) - eta_u(j,i)   ) * ddx

    deta_dy_u(j,i)  = ( eta_v(j+1,i) + eta_v(j+1,i-1) - eta_v(j,i) - eta_v(j,i-1) ) * ddy * 0.5_wp
    deta_dy_v(j,i)  = ( eta_s(j,i)   - eta_s(j-1,i) ) * ddy
    deta_dy_s(j,i)  = ( eta_v(j+1,i) - eta_v(j,i)   ) * ddy

    dzeta_dz_u(j,i) = 1.0_wp / ( H_av - eta_u(j,i)  )
    dzeta_dz_v(j,i) = 1.0_wp / ( H_av - eta_v(j,i)  )
    dzeta_dz_s(j,i) = 1.0_wp / ( H_av - eta_s(j,i)  )

    deta_dxx_u(j,i)  = ( eta_s(j,i) - 2.0_wp * eta_u(j,i) + eta_s(j,i-1) ) * ddx2 * 4.0_wp
    deta_dxx_v(j,i)  = ( 0.5_wp * ( eta_u(j,i+1) + eta_u(j-1,i+1) + eta_u(j,i) + eta_u(j-1,i) )           &
                       - 2.0_wp * eta_v(j,i) ) * ddx2 * 4.0_wp
    deta_dxx_s(j,i)  = ( eta_u(j,i+1) - 2.0_wp * eta_s(j,i) + eta_u(j,i) ) * ddx2 * 4.0_wp

    deta_dyy_u(j,i)  = ( 0.5_wp * ( eta_v(j+1,i) + eta_v(j+1,i-1) + eta_v(j,i) + eta_v(j,i-1) )           &
                       - 2.0_wp * eta_u(j,i) ) * ddy2 * 4.0_wp
    deta_dyy_v(j,i)  = ( eta_s(j,i) - 2.0_wp * eta_v(j,i) + eta_s(j-1,i) ) * ddy2 * 4.0_wp
    deta_dyy_s(j,i)  = ( eta_v(j+1,i) - 2.0_wp * eta_s(j,i) + eta_v(j,i) ) * ddy2 * 4.0_wp

    DO k = nzb, nzt+1
       dzeta_dt_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dt_s(j,i)
       dzeta_dt_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dt_u(j,i)
       dzeta_dt_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dt_v(j,i)
       dzeta_dt_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dt_s(j,i)
       dzeta_dx_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dx_s(j,i)
       dzeta_dx_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dx_u(j,i)
       dzeta_dx_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dx_v(j,i)
       dzeta_dx_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dx_s(j,i)
       dzeta_dy_s(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dy_s(j,i)
       dzeta_dy_u(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) ) * deta_dy_u(j,i)
       dzeta_dy_v(k,j,i)  = ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) ) * deta_dy_v(j,i)
       dzeta_dy_w(k,j,i)  = ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) ) * deta_dy_s(j,i)

       dzeta_dxx_s(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp                 &
                                    * deta_dx_s(j,i)**2.0_wp +                                            &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                         &
                                    * deta_dxx_s(j,i)
       dzeta_dxx_u(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )**2.0_wp                 &
                                    * deta_dx_u(j,i)**2.0_wp  +                                           &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )                         &
                                    * deta_dxx_u(j,i)
       dzeta_dxx_v(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )**2.0_wp                 &
                                    * deta_dx_v(j,i)**2.0_wp  +                                           &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )                         &
                                    * deta_dxx_v(j,i)
       dzeta_dxx_w(k,j,i)  = 2.0_wp * ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp                 &
                                    * deta_dx_s(j,i)**2.0_wp +                                            &
                                      ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                         &
                                    * deta_dxx_s(j,i)
       dzeta_dyy_s(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp                 &
                                    * deta_dy_s(j,i)**2.0_wp +                                            &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                         &
                                    * deta_dyy_s(j,i)
       dzeta_dyy_u(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )**2.0_wp                 &
                                    * deta_dy_u(j,i)**2.0_wp  +                                           &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_u(j,i) )                         &
                                    * deta_dyy_u(j,i)
       dzeta_dyy_v(k,j,i)  = 2.0_wp * ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )**2.0_wp                 &
                                    * deta_dy_v(j,i)**2.0_wp  +                                           &
                                      ( zzu(k) - 1.0_wp ) / ( H_av - eta_v(j,i) )                         &
                                      * deta_dyy_v(j,i)
       dzeta_dyy_w(k,j,i)  = 2.0_wp * ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )**2.0_wp                 &
                                    * deta_dy_s(j,i)**2.0_wp +                                            &
                                      ( zzw(k) - 1.0_wp ) / ( H_av - eta_s(j,i) )                         &
                                    * deta_dyy_s(j,i)

       alpha(k,j,i) = dzeta_dx_s(k,j,i)**2.0_wp + dzeta_dy_s(k,j,i)**2.0_wp + dzeta_dz_s(j,i)**2.0_wp
       beta(k,j,i)  = dzeta_dxx_s(k,j,i) + dzeta_dyy_s(k,j,i)
    ENDDO
    
 END SUBROUTINE calc_metrics_ij


 SUBROUTINE wwinta_exchange_horiz( location )

    USE exchange_horiz_mod,                                                                        &
        ONLY:  exchange_horiz_2d, exchange_horiz

    CHARACTER (LEN=*), INTENT(IN) ::  location

    SELECT CASE ( location )

       CASE ( 'before_prognostic_equation' )

       CASE ( 'after_prognostic_equation' )

          CALL exchange_horiz_2d( eta_u )
          CALL exchange_horiz_2d( eta_v )
          CALL exchange_horiz_2d( eta_s )

          CALL exchange_horiz_2d( spv_u )
          CALL exchange_horiz_2d( spv_v )
          CALL exchange_horiz_2d( spv_w )

       CASE ( 'after_computing_metrics' )

          CALL exchange_horiz_2d( sia_sin_x )
          CALL exchange_horiz_2d( sia_cos_x )
          CALL exchange_horiz_2d( sia_sin_y )
          CALL exchange_horiz_2d( sia_cos_y )

          CALL exchange_horiz_2d( deta_dt_u )
          CALL exchange_horiz_2d( deta_dt_v )
          CALL exchange_horiz_2d( deta_dt_s )

          CALL exchange_horiz_2d( deta_dx_u )
          CALL exchange_horiz_2d( deta_dx_v )
          CALL exchange_horiz_2d( deta_dx_s )
          CALL exchange_horiz_2d( deta_dy_u )
          CALL exchange_horiz_2d( deta_dy_v )
          CALL exchange_horiz_2d( deta_dy_s )

          CALL exchange_horiz_2d( dzeta_dz_u )
          CALL exchange_horiz_2d( dzeta_dz_v )
          CALL exchange_horiz_2d( dzeta_dz_s )

          CALL exchange_horiz( dzeta_dt_s, nbgp )
          CALL exchange_horiz( dzeta_dt_u, nbgp )
          CALL exchange_horiz( dzeta_dt_v, nbgp )
          CALL exchange_horiz( dzeta_dt_w, nbgp )
          CALL exchange_horiz( dzeta_dx_s, nbgp )
          CALL exchange_horiz( dzeta_dx_u, nbgp )
          CALL exchange_horiz( dzeta_dx_v, nbgp )
          CALL exchange_horiz( dzeta_dx_w, nbgp )
          CALL exchange_horiz( dzeta_dy_s, nbgp )
          CALL exchange_horiz( dzeta_dy_u, nbgp )
          CALL exchange_horiz( dzeta_dy_v, nbgp )
          CALL exchange_horiz( dzeta_dy_w, nbgp )

          CALL exchange_horiz_2d( deta_dxx_u )
          CALL exchange_horiz_2d( deta_dxx_v )
          CALL exchange_horiz_2d( deta_dxx_s )
          CALL exchange_horiz_2d( deta_dyy_u )
          CALL exchange_horiz_2d( deta_dyy_v )
          CALL exchange_horiz_2d( deta_dyy_s )

          CALL exchange_horiz( dzeta_dxx_s, nbgp )
          CALL exchange_horiz( dzeta_dxx_u, nbgp )
          CALL exchange_horiz( dzeta_dxx_v, nbgp )
          CALL exchange_horiz( dzeta_dxx_w, nbgp )
          CALL exchange_horiz( dzeta_dyy_s, nbgp )
          CALL exchange_horiz( dzeta_dyy_u, nbgp )
          CALL exchange_horiz( dzeta_dyy_v, nbgp )
          CALL exchange_horiz( dzeta_dyy_w, nbgp )

          CALL exchange_horiz( alpha, nbgp )
          CALL exchange_horiz( beta,  nbgp )

       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE wwinta_exchange_horiz



 SUBROUTINE wwinta_boundary_conditions

    ! the velocity at bottom of the domain should equal to the particle velocity of wave surface.

    u_p(nzb,:,:) = spv_u(:,:)
    v_p(nzb,:,:) = spv_v(:,:)
    w_p(nzb,:,:) = spv_w(:,:)

 END SUBROUTINE wwinta_boundary_conditions


END MODULE wind_wave_inta_mod

