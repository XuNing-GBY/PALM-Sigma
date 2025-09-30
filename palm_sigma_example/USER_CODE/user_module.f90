!> @file user_module.f90
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
!> Declaration of user-defined variables. This module may only be used in the user-defined routines
!> (contained in user_interface.f90).
!--------------------------------------------------------------------------------------------------!
 MODULE user

    USE arrays_3d

    USE control_parameters

    USE cpulog

    USE indices

    USE kinds

    USE pegrid

    USE statistics

    USE surface_mod

    ! BOC - Xu
    USE wind_wave_inta_mod
    ! EOC - Xu

    IMPLICIT NONE

    INTEGER(iwp) ::  dots_num_palm      !<
    INTEGER(iwp) ::  dots_num_user = 0  !<
    INTEGER(iwp) ::  user_idummy        !<

    LOGICAL ::  user_module_enabled = .FALSE.  !<

    REAL(wp) ::  user_rdummy  !<

!
!-- Sample for user-defined output
!    REAL(wp) :: global_parameter  !< user defined global parameter
!
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u2      !< user defined array
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  u2_av   !< user defined array
!    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  ustvst  !< user defined array

    ! BOC - Xu
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wu_sgs       !< SGS-scale momentum flux w"u"
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wu_sgs_av    !< time-averaged SGS-scale momentum flux w"u"
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wv_sgs       !< SGS-scale momentum flux w"v"
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wv_sgs_av    !< time-averaged SGS-scale momentum flux w"v"
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wpt_sgs     !< SGS-scale heat flux w"theta"
    REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, PUBLIC ::  wpt_sgs_av  !< time-averaged SGS-scale heat flux w"theta"
    ! EOC - Xu

    SAVE

    PRIVATE

!
!- Public functions
    PUBLIC                                                                                         &
       user_actions,                                                                               &
       user_check_data_output,                                                                     &
       user_check_data_output_pr,                                                                  &
       user_check_data_output_ts,                                                                  &
       user_check_parameters,                                                                      &
       user_data_output_2d,                                                                        &
       user_data_output_3d,                                                                        &
       user_define_netcdf_grid,                                                                    &
       user_header,                                                                                &
       user_init,                                                                                  &
       user_init_arrays,                                                                           &
       user_last_actions,                                                                          &
       user_parin,                                                                                 &
       user_rrd_global,                                                                            &
       user_rrd_local,                                                                             &
       user_statistics,                                                                            &
       user_3d_data_averaging,                                                                     &
       user_wrd_global,                                                                            &
       user_wrd_local


!
!- Public parameters, constants and initial values
   PUBLIC                                                                                          &
      user_module_enabled

    INTERFACE user_parin
       MODULE PROCEDURE user_parin
    END INTERFACE user_parin

    INTERFACE user_check_parameters
       MODULE PROCEDURE user_check_parameters
    END INTERFACE user_check_parameters

    INTERFACE user_check_data_output_ts
       MODULE PROCEDURE user_check_data_output_ts
    END INTERFACE user_check_data_output_ts

    INTERFACE user_check_data_output_pr
       MODULE PROCEDURE user_check_data_output_pr
    END INTERFACE user_check_data_output_pr

    INTERFACE user_check_data_output
       MODULE PROCEDURE user_check_data_output
    END INTERFACE user_check_data_output

    INTERFACE user_define_netcdf_grid
       MODULE PROCEDURE user_define_netcdf_grid
    END INTERFACE user_define_netcdf_grid

    INTERFACE user_init
       MODULE PROCEDURE user_init
    END INTERFACE user_init

    INTERFACE user_init_arrays
       MODULE PROCEDURE user_init_arrays
    END INTERFACE user_init_arrays

    INTERFACE user_header
       MODULE PROCEDURE user_header
    END INTERFACE user_header

    INTERFACE user_actions
       MODULE PROCEDURE user_actions
       MODULE PROCEDURE user_actions_ij
    END INTERFACE user_actions

    INTERFACE user_3d_data_averaging
       MODULE PROCEDURE user_3d_data_averaging
    END INTERFACE user_3d_data_averaging

    INTERFACE user_data_output_2d
       MODULE PROCEDURE user_data_output_2d
    END INTERFACE user_data_output_2d

    INTERFACE user_data_output_3d
       MODULE PROCEDURE user_data_output_3d
    END INTERFACE user_data_output_3d

    INTERFACE user_statistics
       MODULE PROCEDURE user_statistics
    END INTERFACE user_statistics

    INTERFACE user_rrd_global
       MODULE PROCEDURE user_rrd_global_ftn
       MODULE PROCEDURE user_rrd_global_mpi
    END INTERFACE user_rrd_global

    INTERFACE user_rrd_local
       MODULE PROCEDURE user_rrd_local_ftn
       MODULE PROCEDURE user_rrd_local_mpi
    END INTERFACE user_rrd_local

    INTERFACE user_wrd_global
       MODULE PROCEDURE user_wrd_global
    END INTERFACE user_wrd_global

    INTERFACE user_wrd_local
       MODULE PROCEDURE user_wrd_local
    END INTERFACE user_wrd_local

    INTERFACE user_last_actions
       MODULE PROCEDURE user_last_actions
    END INTERFACE user_last_actions


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Parin for &user_parameters for user module
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_parin

    CHARACTER (LEN=80) ::  line  !< string containing the last line read from namelist file

    INTEGER(iwp) ::  i          !<
    INTEGER(iwp) ::  io_status  !< status after reading the namelist file
    INTEGER(iwp) ::  j          !<

    LOGICAL ::  switch_off_module = .FALSE.  !< local namelist parameter to switch off the module
                                             !< although the respective module namelist appears in
                                             !< the namelist file

    NAMELIST /user_parameters/  data_output_masks_user,                                            &
                                data_output_pr_user,                                               &
                                data_output_user,                                                  &
                                region,                                                            &
                                switch_off_module

!
!-- Next statement is to avoid compiler warnings about unused variables. Please remove in case
!-- that you are using them.
    IF ( dots_num_palm == 0  .OR.  dots_num_user == 0  .OR.  user_idummy == 0  .OR.                &
         user_rdummy == 0.0_wp )  CONTINUE

!
!-- Set revision number of this default interface version. It will be checked within the main
!-- program (palm). Please change the revision number in case that the current revision does not
!-- match with previous revisions (e.g. if routines have been added/deleted or if parameter lists
!-- in subroutines have been changed).
    user_interface_current_revision = 'r4495'

!
!-- Position the namelist-file at the beginning (it has already been opened in parin), and try to
!-- read (find) a namelist named "user_parameters".
    REWIND ( 11 )
    READ( 11, user_parameters, IOSTAT=io_status )

!
!-- Actions depending on the READ status
    IF ( io_status == 0 )  THEN
!
!--    User namelist found and correctly read. Set default module switch to true. This activates
!--    calls of the user-interface subroutines.
       IF ( .NOT. switch_off_module )  user_module_enabled = .TRUE.

    ELSEIF ( io_status > 0 )  THEN
!
!--    User namelist was found, but contained errors. Print an error message containing the line
!--    that caused the problem
       BACKSPACE( 11 )
       READ( 11 , '(A)') line
       CALL parin_fail_message( 'user_parameters', line )

    ENDIF

!
!-- Determine the number of user-defined profiles and append them to the standard data output
!-- (data_output_pr)
    IF ( user_module_enabled )  THEN
       IF ( data_output_pr_user(1) /= ' ' )  THEN
          i = 1
          DO WHILE ( data_output_pr(i) /= ' '  .AND.  i <= 200 )
             i = i + 1
          ENDDO
          j = 1
          DO WHILE ( data_output_pr_user(j) /= ' '  .AND.  j <= 300 )
             data_output_pr(i) = data_output_pr_user(j)
             max_pr_user_tmp   = max_pr_user_tmp + 1
             i = i + 1
             j = j + 1
          ENDDO
       ENDIF
    ENDIF


 END SUBROUTINE user_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check &userpar control parameters and deduce further quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_parameters

!
!-- Here the user may add code to check the validity of further &userpar control parameters or
!-- deduce further quantities.


 END SUBROUTINE user_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set module-specific timeseries units and labels
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp),      INTENT(IN)     ::  dots_max  !<
    INTEGER(iwp),      INTENT(INOUT)  ::  dots_num  !<

    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_label  !<
    CHARACTER(LEN=*), DIMENSION(dots_max), INTENT(INOUT)  ::  dots_unit   !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( dots_num == 0  .OR.  dots_label(1)(1:1) == ' '  .OR.  dots_unit(1)(1:1) == ' ' )  CONTINUE

!
!-- Sample for user-defined time series:
!-- For each time series quantity you have to give a label and a unit, which will be used for the
!-- NetCDF file. They must not contain more than seven characters. The value of dots_num has to be
!-- increased by the number of new time series quantities. Its old value has to be stored in
!-- dots_num_palm. See routine user_statistics on how to calculate and output these quantities.

!    dots_num_palm = dots_num

!    dots_num = dots_num + 1
!    dots_num_user = dots_num_user + 1
!    dots_label(dots_num) = 'abs_umx'
!    dots_unit(dots_num)  = 'm/s'

!    dots_num = dots_num + 1
!    dots_num_user = dots_num_user + 1
!    dots_label(dots_num) = 'abs_vmx'
!    dots_unit(dots_num)  = 'm/s'


 END SUBROUTINE user_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined profile output quantities. For those variables not recognized by the
!> user, the parameter unit is set to "illegal", which tells the calling routine that the
!> output variable is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output_pr( variable, var_count, unit, dopr_unit )


    USE profil_parameter


    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<
    CHARACTER (LEN=*) ::  dopr_unit !< local value of dopr_unit

    INTEGER(iwp) ::  user_pr_index  !<
    INTEGER(iwp) ::  var_count      !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( unit(1:1) == ' '  .OR.  dopr_unit(1:1) == ' '  .OR.  var_count == 0 )  CONTINUE

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    Add additional CASE statements depending on the number of quantities for which profiles are
!--    to be calculated. The respective calculations to be performed have to be added in routine
!--    user_statistics. The quantities are (internally) identified by a user-profile-number
!--    (see variable "user_pr_index" below). The first user-profile must be assigned the number
!--    "pr_palm+1", the second one "pr_palm+2", etc. The respective user-profile-numbers have also
!--    to be used in routine user_statistics!
!       CASE ( 'u*v*' )                      ! quantity string as given in data_output_pr_user
!          user_pr_index = pr_palm + 1
!          dopr_index(var_count)  = user_pr_index    ! quantities' user-profile-number
!          dopr_unit = 'm2/s2'  ! quantity unit
!          unit = dopr_unit
!          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 )
!                                            ! grid on which the quantity is defined (use zu or zw)
!

       CASE ( 'advcx' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 1
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'advcy' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 2
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'advcz' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 3
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid
 
       CASE ( 'dprex' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 4
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'dprez' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 5
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'diffx' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 6
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'diffy' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 7
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'diffz' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 8
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm/s2'                      ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'div' )                             ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 9
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = '1/s'                       ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE ( 'mflpz' )                           ! quantity string as given in data_output_pr_user
          user_pr_index = pr_palm + 10
          dopr_index(var_count)  = user_pr_index  ! quantities' user-profile-number
          dopr_unit = 'm2/s2'                     ! quantity unit
          unit = dopr_unit
          hom(:,2,user_pr_index,:) = SPREAD( zu, 2, statistic_regions+1 ) ! defined on u-grid

       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE user_check_data_output_pr


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the unit of user defined output quantities. For those variables not recognized by the user,
!> the parameter unit is set to "illegal", which tells the calling routine that the output variable
!> is not defined and leads to a program abort.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_check_data_output( variable, unit )


    CHARACTER (LEN=*) ::  unit      !<
    CHARACTER (LEN=*) ::  variable  !<


    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary
!       CASE ( 'u2' )
!          unit = 'm2/s2'
!
!       CASE ( 'u*v*' )
!          unit = 'm2/s2'
!
        ! BOC - Xu
        CASE ( 'wu_sgs' )
           unit = 'm2/s2'

        CASE ( 'wv_sgs' )
           unit = 'm2/s2'

        CASE ( 'wpt_sgs' )
           unit = 'K*m/s'
        ! EOC - Xu
       CASE DEFAULT
          unit = 'illegal'

    END SELECT


 END SUBROUTINE user_check_data_output


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize user-defined arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init_arrays

 ! BOC - Xu
    INTEGER(iwp) :: i       !< loop index
    INTEGER(iwp) :: j       !< loop index
    INTEGER(iwp) :: region  !< index for loop over statistic regions

!
!-- Allocate user-defined arrays and set flags for statistic regions.
!-- Sample for user-defined output
    ALLOCATE( wu_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wv_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
    ALLOCATE( wpt_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!
!-- Example for defining a statistic region:
!-- ATTENTION: rmask = 0 is required at the ghost boundaries to guarantee correct statistic
!--            evaluations (otherwise ghost points would be counted twice). This setting has
!--            already been cared for in routine init_3d_model. Please don't set the ghost points
!--            /= 0. i.e. run the following loop only over nxl,nxr and nys,nyn.
    !IF ( statistic_regions >= 1 )  THEN
    !   region = 1
    !
    !   rmask(:,:,region) = 0.0_wp
    !   DO  i = nxl, nxr
    !      IF ( i >= INT( 0.25 * nx ) .AND. i <= INT( 0.75 * nx ) )  THEN
    !         DO  j = nys, nyn
    !            IF ( i >= INT( 0.25 * ny ) .AND. i <= INT( 0.75 * ny ) )  THEN
    !               rmask(j,i,region) = 1.0_wp
    !            ENDIF
    !         ENDDO
    !      ENDIF
    !   ENDDO
    !ENDIF
 ! EOC - Xu

 END SUBROUTINE user_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined initializing actions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_init


!    CHARACTER(LEN=20) :: field_char  !<
!
!-- Here the user-defined initializing actions follow:
!-- Sample for user-defined output
!    ustvst = 0.0_wp

 ! BOC - Xu
 wu_sgs   = 0.0_wp
 wv_sgs   = 0.0_wp
 wpt_sgs  = 0.0_wp
 ! EOC - Xu

 END SUBROUTINE user_init


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Set the grids on which user-defined output quantities are defined. Allowed values for grid_x are
!> "x" and "xu", for grid_y "y" and "yv", and for grid_z "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_define_netcdf_grid( variable, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*) ::  grid_x     !<
    CHARACTER (LEN=*) ::  grid_y     !<
    CHARACTER (LEN=*) ::  grid_z     !<
    CHARACTER (LEN=*) ::  variable   !<

    LOGICAL ::  found   !<


    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary
!       CASE ( 'u2', 'u2_xy', 'u2_xz', 'u2_yz' )
!          found  = .TRUE.
!          grid_x = 'xu'
!          grid_y = 'y'
!          grid_z = 'zu'

!       CASE ( 'u*v*', 'u*v*_xy', 'u*v*_xz', 'u*v*_yz' )
!          found  = .TRUE.
!          grid_x = 'x'
!          grid_y = 'y'
!          grid_z = 'zu'

       ! BOC - Xu
       CASE ( 'wu_sgs', 'wu_sgs_xy', 'wu_sgs_xz', 'wu_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'wv_sgs', 'wv_sgs_xy', 'wv_sgs_xz', 'wv_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'

       CASE ( 'wpt_sgs', 'wpt_sgs_xy', 'wpt_sgs_xz', 'wpt_sgs_yz' )
          found  = .TRUE.
          grid_x = 'x'
          grid_y = 'y'
          grid_z = 'zw'
       ! EOC - Xu

       CASE DEFAULT
          found  = .FALSE.
          grid_x = 'none'
          grid_y = 'none'
          grid_z = 'none'

    END SELECT


 END SUBROUTINE user_define_netcdf_grid




!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Print a header with user-defined information.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_header( io )


    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  io  !<

!
!-- If no user-defined variables are read from the namelist-file, no information will be printed.
    IF ( .NOT. user_module_enabled )  THEN
       WRITE ( io, 100 )
       RETURN
    ENDIF

!
!-- Printing the information.
    WRITE ( io, 110 )

    IF ( statistic_regions /= 0 )  THEN
       WRITE ( io, 200 )
       DO  i = 0, statistic_regions
          WRITE ( io, 201 )  i, region(i)
       ENDDO
    ENDIF

!
!-- Format-descriptors
100 FORMAT (//' *** no user-defined variables found'/)
110 FORMAT (//1X,78('#') // ' User-defined variables and actions:' /                               &
            ' -----------------------------------'//)
200 FORMAT (' Output of profiles and time series for following regions:' /)
201 FORMAT (4X,'Region ',I1,':   ',A)


 END SUBROUTINE user_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for all grid points
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_actions( location )

    USE arrays_3d,                                                                              &
        ONLY:  ddzu

    USE grid_variables,                                                                         &
        ONLY:  ddx, ddy

    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    CALL cpu_log( log_point(24), 'user_actions', 'start' )

!
!-- Here the user-defined actions follow. No calls for single grid points are allowed at locations
!-- before and after the timestep, since these calls are not within an i,j-loop
    SELECT CASE ( location )

       CASE ( 'before_timestep' )
!
!--       Enter actions to be done before every timestep here

       CASE ( 'before_prognostic_equations' )
!
!--       Enter actions to be done before all prognostic equations here

       CASE ( 'after_integration' )
!
!--       Enter actions to be done after every time integration (before data output)
!--       Sample for user-defined output:
!          DO  i = nxlg, nxrg
!             DO  j = nysg, nyng
!                DO  k = nzb, nzt
!                   u2(k,j,i) = u(k,j,i)**2
!                ENDDO
!             ENDDO
!          ENDDO
!          DO  i = nxlg, nxr
!             DO  j = nysg, nyn
!                DO  k = nzb, nzt+1
!                   ustvst(k,j,i) =  &
!                      ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1,0) ) *                      &
!                      ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2,0) )
!                ENDDO
!             ENDDO
!          ENDDO
          ! BOC - Xu
          DO  i = nxlg, nxr
             DO  j = nysg, nyn
                DO  k = nzb, nzt+1
                   wu_sgs(k,j,i) =   &
                                                  - 0.25_wp * (                                    &
                                  km(k,j,i) + km(k+1,j,i) + km(k,j,i-1) + km(k+1,j,i-1)            &
                                                              ) * (                                &
                                      ( u(k+1,j,i) - u(k,j,i)   ) * ddzu(k+1)                      &
                                    + ( w(k,j,i)   - w(k,j,i-1) ) * ddx                            &
                                                                  )                                &
                                            * rho_air_zw(k)                                        &
                                            * momentumflux_output_conversion(k)                    
                   wv_sgs(k,j,i) =   &
                                                  - 0.25_wp * (                                    &
                                  km(k,j,i) + km(k+1,j,i) + km(k,j-1,i) + km(k+1,j-1,i)            &
                                                              ) * (                                &
                                      ( v(k+1,j,i) - v(k,j,i)   ) * ddzu(k+1)                      &
                                    + ( w(k,j,i)   - w(k,j-1,i) ) * ddy                            &
                                                                  )                                &
                                            * rho_air_zw(k)                                        &
                                            * momentumflux_output_conversion(k)

                   wpt_sgs(k,j,i) =  &
                                                  - 0.5_wp * ( kh(k,j,i) + kh(k+1,j,i) )           &
                                               * ( pt(k+1,j,i) - pt(k,j,i) )                       &
                                               * rho_air_zw(k)                                     &
                                               * heatflux_output_conversion(k)                     &
                                               * ddzu(k+1)
                ENDDO
             ENDDO
          ENDDO
          ! EOC - Xu

       CASE ( 'after_timestep' )
!
!--       Enter actions to be done after every timestep here


       CASE ( 'u-tendency' )
!
!--       Enter actions to be done in the u-tendency term here


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

    CALL cpu_log( log_point(24), 'user_actions', 'stop' )

 END SUBROUTINE user_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Call for grid point i,j
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_actions_ij( i, j, location )


    CHARACTER(LEN=*) ::  location  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<

!
!-- Here the user-defined actions follow
    SELECT CASE ( location )

       CASE ( 'u-tendency' )

!
!--       Next line is to avoid compiler warning about unused variables. Please remove.
          IF ( i == 0  .OR.  j == 0 )  CONTINUE

!
!--       Enter actions to be done in the u-tendency term here


       CASE ( 'v-tendency' )


       CASE ( 'w-tendency' )


       CASE ( 'pt-tendency' )


       CASE ( 'sa-tendency' )


       CASE ( 'e-tendency' )


       CASE ( 'q-tendency' )


       CASE ( 's-tendency' )


       CASE DEFAULT
          CONTINUE

    END SELECT

 END SUBROUTINE user_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Sum up and time-average user-defined output quantities as well as allocate the array necessary
!> for storing the average.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_3d_data_averaging( mode, variable )


    CHARACTER(LEN=*) ::  mode      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  i  !<
    INTEGER(iwp) ::  j  !<
    INTEGER(iwp) ::  k  !<

    IF ( mode == 'allocate' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have to be declared and
!--       defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( .NOT. ALLOCATED( u2_av ) )  THEN
!                ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!             ENDIF
!             u2_av = 0.0_wp

          ! BOC - Xu
          CASE ( 'wu_sgs' )
              IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                 ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
              ENDIF
              wu_sgs_av = 0.0_wp

          CASE ( 'wv_sgs' )
              IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                 ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
              ENDIF
              wv_sgs_av = 0.0_wp

          CASE ( 'wpt_sgs' )
              IF ( .NOT. ALLOCATED( wpt_sgs_av ) )  THEN
                 ALLOCATE( wpt_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
              ENDIF
              wpt_sgs_av = 0.0_wp    
          ! EOC - Xu

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'sum' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--       and defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( ALLOCATED( u2_av ) )  THEN
!                DO  i = nxlg, nxrg
!                   DO  j = nysg, nyng
!                      DO  k = nzb, nzt+1
!                         u2_av(k,j,i) = u2_av(k,j,i) + u2(k,j,i)
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF

          ! BOC - Xu
          CASE ( 'wu_sgs' )
             IF ( ALLOCATED( wu_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wu_sgs_av(k,j,i) = wu_sgs_av(k,j,i) + wu_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_sgs' )
             IF ( ALLOCATED( wv_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wv_sgs_av(k,j,i) = wv_sgs_av(k,j,i) + wv_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wpt_sgs' )
             IF ( ALLOCATED( wpt_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wpt_sgs_av(k,j,i) = wpt_sgs_av(k,j,i) + wpt_sgs(k,j,i)
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ! EOC - Xu

          CASE DEFAULT
             CONTINUE

       END SELECT

    ELSEIF ( mode == 'average' )  THEN

       SELECT CASE ( TRIM( variable ) )

!
!--       Uncomment and extend the following lines, if necessary.
!--       The arrays for storing the user defined quantities (here u2_av) have to be declared and
!--       defined by the user!
!--       Sample for user-defined output:
!          CASE ( 'u2' )
!             IF ( ALLOCATED( u2_av ) )  THEN
!                DO  i = nxlg, nxrg
!                   DO  j = nysg, nyng
!                      DO  k = nzb, nzt+1
!                         u2_av(k,j,i) = u2_av(k,j,i) / REAL( average_count_3d, KIND=wp )
!                      ENDDO
!                   ENDDO
!                ENDDO
!             ENDIF

          ! BOC - Xu
          CASE ( 'wu_sgs' )
             IF ( ALLOCATED( wu_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wu_sgs_av(k,j,i) = wu_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wv_sgs' )
             IF ( ALLOCATED( wv_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wv_sgs_av(k,j,i) = wv_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF

          CASE ( 'wpt_sgs' )
             IF ( ALLOCATED( wpt_sgs_av ) )  THEN
                DO  i = nxlg, nxrg
                   DO  j = nysg, nyng
                      DO  k = nzb, nzt+1
                         wpt_sgs_av(k,j,i) = wpt_sgs_av(k,j,i) / REAL( average_count_3d, KIND=wp )
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ! EOC - Xu

       END SELECT

    ENDIF


 END SUBROUTINE user_3d_data_averaging


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k) and sets the grid on which it is defined. Allowed values for grid are "zu" and "zw".
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )


    CHARACTER(LEN=*) ::  grid      !<
    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av      !< flag to control data output of instantaneous or time-averaged data
    INTEGER(iwp) ::  i       !< grid index along x-direction
    INTEGER(iwp) ::  j       !< grid index along y-direction
    INTEGER(iwp) ::  k       !< grid index along z-direction
!    INTEGER(iwp) ::  m       !< running index surface elements
    INTEGER(iwp) ::  nzb_do  !< lower limit of the domain (usually nzb)
    INTEGER(iwp) ::  nzt_do  !< upper limit of the domain (usually nzt+1)

    LOGICAL      ::  found  !<
    LOGICAL      ::  two_d  !< flag parameter that indicates 2D variables (horizontal cross sections)

    REAL(wp) ::  fill_value = -999.0_wp    !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp  .OR.  two_d )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--    and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2_xy', 'u2_xz', 'u2_yz' )
!          IF ( av == 0 )  THEN
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ELSE
!             IF ( .NOT. ALLOCATED( u2_av ) )  THEN
!                ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!                u2_av = REAL( fill_value, KIND = wp )
!             ENDIF
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2_av(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ENDIF
!
!          grid = 'zu'
!
!--    In case two-dimensional surface variables are output, the user has to access related
!--    surface-type. Uncomment and extend following lines appropriately (example output of vertical
!--    surface momentum flux of u-component). Please note, surface elements can be distributed over
!--    several data types, depending on their respective surface properties.
!       CASE ( 'usws_xy' )
!          IF ( av == 0 )  THEN
!
!--           Horizontal default-type surfaces
!             DO  m = 1, surf_def_h(0)%ns
!                i = surf_def_h(0)%i(m)
!                j = surf_def_h(0)%j(m)
!                local_pf(i,j,1) = surf_def_h(0)%usws(m)
!             ENDDO
!
!--           Horizontal natural-type surfaces
!             DO  m = 1, surf_lsm_h%ns
!                i = surf_lsm_h%i(m)
!                j = surf_lsm_h%j(m)
!                local_pf(i,j,1) = surf_lsm_h%usws(m)
!             ENDDO
!
!--           Horizontal urban-type surfaces
!             DO  m = 1, surf_usm_h%ns
!                i = surf_usm_h%i(m)
!                j = surf_usm_h%j(m)
!                local_pf(i,j,1) = surf_usm_h%usws(m)
!             ENDDO
!          ENDIF
!
!          grid = 'zu'
!--

       ! BOC - Xu
       CASE ( 'wu_sgs_xy', 'wu_sgs_xz', 'wu_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wv_sgs_xy', 'wv_sgs_xz', 'wv_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'

       CASE ( 'wpt_sgs_xy', 'wpt_sgs_xz', 'wpt_sgs_yz' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wpt_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wpt_sgs_av ) )  THEN
                ALLOCATE( wpt_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wpt_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wpt_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

          grid = 'zw'
       ! EOC - Xu

       CASE DEFAULT
          found = .FALSE.
          grid  = 'none'

    END SELECT


 END SUBROUTINE user_data_output_2d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Resorts the user-defined output quantity with indices (k,j,i) to a temporary array with indices
!> (i,j,k).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )


    CHARACTER(LEN=*) ::  variable  !<

    INTEGER(iwp) ::  av     !<
    INTEGER(iwp) ::  i      !<
    INTEGER(iwp) ::  j      !<
    INTEGER(iwp) ::  k      !<
    INTEGER(iwp) ::  nzb_do !< lower limit of the data output (usually 0)
    INTEGER(iwp) ::  nzt_do !< vertical upper limit of the data output (usually nz_do3d)

    LOGICAL      ::  found  !<

    REAL(wp) ::  fill_value = -999.0_wp  !< value for the _FillValue attribute

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do) ::  local_pf  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( av == 0  .OR.  local_pf(nxl,nys,nzb_do) == 0.0_wp )  CONTINUE


    found = .TRUE.

    SELECT CASE ( TRIM( variable ) )

!
!--    Uncomment and extend the following lines, if necessary.
!--    The arrays for storing the user defined quantities (here u2 and u2_av) have to be declared
!--    and defined by the user!
!--    Sample for user-defined output:
!       CASE ( 'u2' )
!          IF ( av == 0 )  THEN
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ELSE
!             IF ( .NOT. ALLOCATED( u2_av ) )  THEN
!                ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!                u2_av = REAL( fill_value, KIND = wp )
!             ENDIF
!             DO  i = nxl, nxr
!                DO  j = nys, nyn
!                   DO  k = nzb_do, nzt_do
!                      local_pf(i,j,k) = u2_av(k,j,i)
!                   ENDDO
!                ENDDO
!             ENDDO
!          ENDIF
!

       ! BOC - Xu
       CASE ( 'wu_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
                ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wu_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wu_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wv_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
                ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wv_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wv_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF

       CASE ( 'wpt_sgs' )
          IF ( av == 0 )  THEN
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wpt_sgs(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF ( .NOT. ALLOCATED( wpt_sgs_av ) )  THEN
                ALLOCATE( wpt_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
                wpt_sgs_av = REAL( fill_value, KIND = wp )
             ENDIF
             DO  i = nxl, nxr
                DO  j = nys, nyn
                   DO  k = nzb_do, nzt_do
                      local_pf(i,j,k) = wpt_sgs_av(k,j,i)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF          
       ! EOC - Xu

       CASE DEFAULT
          found = .FALSE.

    END SELECT


 END SUBROUTINE user_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Calculation of user-defined statistics, i.e. horizontally averaged profiles and time series.
!> This routine is called for every statistic region sr defined by the user, but at least for the
!> region "total domain" (sr=0). See section 3.5.4 on how to define, calculate, and output user
!> defined quantities.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_statistics( mode, sr, tn )


    CHARACTER(LEN=*) ::  mode  !<
    INTEGER(iwp) ::  i   !<
    INTEGER(iwp) ::  j   !<
    INTEGER(iwp) ::  k   !<
    INTEGER(iwp) ::  sr  !<
    INTEGER(iwp) ::  tn  !<

    ! BOC - Xu
    INTEGER(iwp) :: m
    INTEGER(iwp) :: surf_s
    INTEGER(iwp) :: surf_e

    REAL(wp)     ::  kmym
    REAL(wp)     ::  kmyp
    REAL(wp)     ::  kmzm
    REAL(wp)     ::  kmzp
    ! EOC - Xu

!    REAL(wp), DIMENSION(:), ALLOCATABLE ::  ts_value_l  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    IF ( sr == 0  .OR.  tn == 0 )  CONTINUE

    IF ( mode == 'profiles' )  THEN

!
!--    Sample on how to calculate horizontally averaged profiles of user-defined quantities. Each
!--    quantity is identified by the index "pr_palm+#" where "#" is an integer starting from 1.
!--    These user-profile-numbers must also be assigned to the respective strings given by
!--    data_output_pr_user in routine user_check_data_output_pr.
!       !$OMP DO
!       DO  i = nxl, nxr
!          DO  j = nys, nyn
!             DO  k = nzb+1, nzt
!!
!!--             Sample on how to calculate the profile of the resolved-scale horizontal momentum
!!--             flux u*v*
!                sums_l(k,pr_palm+1,tn) = sums_l(k,pr_palm+1,tn) +                                  &
!                                         ( 0.5_wp * ( u(k,j,i) + u(k,j,i+1) ) - hom(k,1,1,sr) ) *  &
!                                         ( 0.5_wp * ( v(k,j,i) + v(k,j+1,i) ) - hom(k,1,2,sr) ) *  &
!                                         rmask(j,i,sr) * MERGE( 1.0_wp, 0.0_wp,                    &
!                                         BTEST( topo_flags(k,j,i), 0 ) )
!!
!!--             Further profiles can be defined and calculated by increasing the second index of
!!--             array sums_l (replace ... appropriately)
!                sums_l(k,pr_palm+2,tn) = sums_l(k,pr_palm+2,tn) + ...   * rmask(j,i,sr)
!             ENDDO
!          ENDDO
!       ENDDO

       !$OMP DO
       DO  i = nxl, nxr
          DO  j = nys, nyn
             DO  k = nzb+1, nzt
                ! BOC - Xu
                ! advcx
                sums_l(k,pr_palm+1,tn) = sums_l(k,pr_palm+1,tn)                                    &
                                       + ( ( ( u(k,j,i) + u(k,j,i+1) ) *                           &
                                             ( u(k,j,i) + u(k,j,i+1) ) * 0.25_wp                   &
                                           - ( u(k,j,i) + u(k,j,i-1) ) *                           &
                                             ( u(k,j,i) + u(k,j,i-1) ) * 0.25_wp                   &
                                           ) * ddx                                                 &
                                         ) * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! advcy
                sums_l(k,pr_palm+2,tn) = sums_l(k,pr_palm+2,tn)                                    &
                                       + ( ( ( v(k,j+1,i-1)   + v(k,j+1,i) ) *                     &
                                             ( u(k,j+1,i)     + u(k,j,i)   ) * 0.25_wp             &
                                           - ( v(k,j,i-1)     + v(k,j,i)   ) *                     &
                                             ( u(k,j-1,i)     + u(k,j,i)   ) * 0.25_wp             &
                                           ) * ddy                                                 &
                                         ) * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! advcz
                sums_l(k,pr_palm+3,tn) = sums_l(k,pr_palm+3,tn)                                    &
                                       + ( ( ( u(k,j,i) + u(k+1,j,i) ) *                           &
                                             ( u(k,j,i) + u(k+1,j,i) ) * 0.25_wp                   &
                                           - ( u(k,j,i) + u(k-1,j,i) ) *                           &
                                             ( u(k,j,i) + u(k-1,j,i) ) * 0.25_wp                   &
                                             * MERGE( 1.0_wp, 0.0_wp,                              &
                                                      BTEST( topo_flags(k-1,j,i), 8 ) )            &
                                           ) * ddzzw(k) * dzeta_dx_u(k,j,i)                        &
                                         + ( ( v(k,j,i)       + v(k,j+1,i)                         &
                                             + v(k,j+1,i-1)   + v(k,j,i-1)                         &
                                             + v(k+1,j,i)     + v(k+1,j+1,i)                       &
                                             + v(k+1,j+1,i-1) + v(k+1,j,i-1) ) *                   &
                                             ( u(k+1,j,i) + u(k,j,i) ) * 0.0625_wp                 &
                                           - ( v(k,j,i)       + v(k,j+1,i)                         &
                                             + v(k,j+1,i-1)   + v(k,j,i-1)                         &
                                             + v(k-1,j,i)     + v(k-1,j+1,i)                       &
                                             + v(k-1,j+1,i-1) + v(k-1,j,i-1) ) *                   &
                                             ( u(k-1,j,i) + u(k,j,i) ) * 0.0625_wp                 &
                                             * MERGE( 1.0_wp, 0.0_wp,                              &
                                                      BTEST( topo_flags(k-1,j,i), 8 ) )            &
                                           ) * ddzzw(k) * dzeta_dy_u(k,j,i)                        &
                                         + ( ( w(k,j,i-1)     + w(k,j,i)   ) *                     &
                                             ( u(k+1,j,i)     + u(k,j,i)   ) * 0.25_wp             &
                                           - ( w(k-1,j,i-1)   + w(k-1,j,i) ) *                     &
                                             ( u(k-1,j,i)     + u(k,j,i)   ) * 0.25_wp             &
                                             * MERGE( 1.0_wp, 0.0_wp,                              &
                                                      BTEST( topo_flags(k-1,j,i), 8 ) )            &
                                           ) * ddzzw(k) * dzeta_dz_u(j,i)                          &
                                         ) * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! dprex
                sums_l(k,pr_palm+4,tn) = sums_l(k,pr_palm+4,tn)                                    &
                                       + ( ( p(k,j,i) - p(k,j,i-1)                                 &
                                           ) * drho_air(k) * ddx                                   &
                                         ) * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! dprez
                sums_l(k,pr_palm+5,tn) = sums_l(k,pr_palm+5,tn)                                    &
                                       + ( ( ( p(k+1,j,i) + p(k+1,j,i-1)                           &
                                             + p(k,j,i)   + p(k,j,i-1)   ) * 0.25_wp               &
                                             * drho_air_zw(k)                                      &
                                           - ( p(k-1,j,i) + p(k-1,j,i-1)                           &
                                             + p(k,j,i)   + p(k,j,i-1)   ) * 0.25_wp               &
                                             * drho_air_zw(k-1)                                    &
                                             * MERGE( 1.0_wp, 0.0_wp,                              &
                                                      BTEST( topo_flags(k-1,j,i), 8 ) )            &
                                           ) * ddzzw(k) * dzeta_dx_u(k,j,i)                        &
                                         ) * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                kmyp = 0.25_wp * ( km(k,j,i)+km(k,j+1,i)+km(k,j,i-1)+km(k,j+1,i-1) )
                kmym = 0.25_wp * ( km(k,j,i)+km(k,j-1,i)+km(k,j,i-1)+km(k,j-1,i-1) )
                kmzp = 0.25_wp * ( km(k,j,i)+km(k+1,j,i)+km(k,j,i-1)+km(k+1,j,i-1) )
                kmzm = 0.25_wp * ( km(k,j,i)+km(k-1,j,i)+km(k,j,i-1)+km(k-1,j,i-1) )

                ! diffx
                sums_l(k,pr_palm+6,tn) = sums_l(k,pr_palm+6,tn)                                    &
                      + ( km(k,j,i)   * ( ( u(k,j,i+1) - u(k,j,i) ) * ddx                          &
                                        + ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1) + u(k,j,i) )  &
                                          - ( u(k-1,j,i+1) + u(k-1,j,i) + u(k,j,i+1) + u(k,j,i) )  &
                                            * MERGE( 1.0_wp, 0.0_wp,                               &
                                                     BTEST( topo_flags(k-1,j,i), 8 ) )             &
                                          ) * 0.25_wp * ddzzw(k) * dzeta_dx_s(k,j,i)               &
                                        )                                                          &
                        - km(k,j,i-1) * ( ( u(k,j,i) - u(k,j,i-1) ) * ddx                          &
                                        + ( ( u(k+1,j,i-1) + u(k+1,j,i) + u(k,j,i-1) + u(k,j,i) )  &
                                          - ( u(k-1,j,i-1) + u(k-1,j,i) + u(k,j,i-1) + u(k,j,i) )  &
                                            * MERGE( 1.0_wp, 0.0_wp,                               &
                                                     BTEST( topo_flags(k-1,j,i), 8 ) )             &
                                          ) * 0.25_wp * ddzzw(k) * dzeta_dx_s(k,j,i-1)             &
                                        )                                                          &
                        ) * 2.0_wp * ddx                                                           &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                 ! diffy
                 sums_l(k,pr_palm+7,tn) = sums_l(k,pr_palm+7,tn)                                   &
                      + ( kmyp * ( ( u(k,j+1,i) - u(k,j,i) ) * ddy                                 &
                                 + ( ( u(k,j,i) + u(k,j+1,i) + u(k+1,j,i) + u(k+1,j+1,i) )         &
                                   - ( u(k,j,i) + u(k,j+1,i) + u(k-1,j,i) + u(k-1,j+1,i) )         &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )    &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dy_u(k,j,i) + dzeta_dy_u(k,j+1,i) ) * 0.5_wp          &
                                 + ( v(k,j+1,i) - v(k,j+1,i-1) ) * ddx                             &
                                 + ( ( v(k,j+1,i) + v(k,j+1,i-1) + v(k+1,j+1,i) + v(k+1,j+1,i-1) ) &
                                   - ( v(k,j+1,i) + v(k,j+1,i-1) + v(k-1,j+1,i) + v(k-1,j+1,i-1) ) &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )    &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k,j+1,i) ) * 0.5_wp          &
                                 )                                                                 &
                        - kmym * ( ( u(k,j,i) - u(k,j-1,i) ) * ddy                                 &
                                 + ( ( u(k,j,i) + u(k,j-1,i) + u(k+1,j,i) + u(k+1,j-1,i) )         &
                                   - ( u(k,j,i) + u(k,j-1,i) + u(k-1,j,i) + u(k-1,j-1,i) )         &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )    &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dy_u(k,j,i) + dzeta_dy_u(k,j-1,i) ) * 0.5_wp          &
                                 + ( v(k,j,i) - v(k,j,i-1) ) * ddx                                 &
                                 + ( ( v(k,j,i) + v(k,j,i-1) + v(k+1,j,i) + v(k+1,j,i-1) )         &
                                   - ( v(k,j,i) + v(k,j,i-1) + v(k-1,j,i) + v(k-1,j,i-1) )         &
                                     * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )    &
                                   ) * 0.25_wp * ddzzw(k)                                          &
                                   * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k,j-1,i) ) * 0.5_wp          &
                                 )                                                                 &
                        ) * ddy * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! diffz
                sums_l(k,pr_palm+8,tn) = sums_l(k,pr_palm+8,tn)                                    &
                      + ( kmzp * ( ( ( u(k+1,j,i+1) + u(k+1,j,i) + u(k,j,i+1) + u(k,j,i) )         &
                                   - ( u(k+1,j,i-1) + u(k+1,j,i) + u(k,j,i-1) + u(k,j,i) )         &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( u(k+1,j,i) - u(k,j,i) ) * ddzzu(k+1)                          &
                                   * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k+1,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k)                                                 &
                        - kmzm * ( ( ( u(k-1,j,i+1) + u(k-1,j,i) + u(k,j,i+1) + u(k,j,i) )         &
                                   - ( u(k-1,j,i-1) + u(k-1,j,i) + u(k,j,i-1) + u(k,j,i) )         &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( u(k,j,i) - u(k-1,j,i) ) * ddzzu(k)                            &
                                   * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k-1,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1)                                               &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )      &
                        ) * 2.0_wp * ddzzw(k) * drho_air(k) * dzeta_dx_u(k,j,i)                    &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                 &
                      + ( kmzp * ( ( ( u(k,j,i) + u(k,j+1,i) + u(k+1,j,i) + u(k+1,j+1,i) )         &
                                   - ( u(k,j,i) + u(k,j-1,i) + u(k+1,j,i) + u(k+1,j-1,i) )         &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( u(k+1,j,i) - u(k,j,i) ) * ddzzu(k+1)                          &
                                   * ( dzeta_dy_u(k+1,j,i) + dzeta_dy_u(k,j,i) ) * 0.5_wp          &
                                 + ( ( v(k,j,i) + v(k,j+1,i) + v(k+1,j,i) + v(k+1,j+1,i) )         &
                                   - ( v(k,j,i-1) + v(k,j+1,i-1) + v(k+1,j,i-1) + v(k+1,j+1,i-1) ) &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( ( v(k+1,j,i) + v(k+1,j,i-1) + v(k+1,j+1,i) + v(k+1,j+1,i-1) ) &
                                   - ( v(k,j,i) + v(k,j,i-1) + v(k,j+1,i) + v(k,j+1,i-1) )         &
                                   ) * 0.25_wp * ddzzu(k+1)                                        &
                                   * ( dzeta_dy_u(k+1,j,i) + dzeta_dy_u(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k)                                                 &
                        - kmzm * ( ( ( u(k,j,i) + u(k,j+1,i) + u(k-1,j,i) + u(k-1,j+1,i) )         &
                                   - ( u(k,j,i) + u(k,j-1,i) + u(k-1,j,i) + u(k-1,j-1,i) )         &
                                   ) * 0.25_wp * ddy                                               &
                                 + ( u(k,j,i) - u(k-1,j,i) ) * ddzzu(k)                            &
                                   * ( dzeta_dy_u(k-1,j,i) + dzeta_dy_u(k,j,i) ) * 0.5_wp          &
                                 + ( ( v(k,j,i) + v(k,j+1,i) + v(k-1,j,i) + v(k-1,j+1,i) )         &
                                   - ( v(k,j,i-1) + v(k,j+1,i-1) + v(k-1,j,i-1) + v(k-1,j+1,i-1) ) &
                                   ) * 0.25_wp * ddx                                               &
                                 + ( ( v(k,j,i) + v(k,j,i-1) + v(k,j+1,i) + v(k,j+1,i-1) )         &
                                   - ( v(k-1,j,i) + v(k-1,j,i-1) + v(k-1,j+1,i) + v(k-1,j+1,i-1) ) &
                                   ) * 0.25_wp * ddzzu(k)                                          &
                                   * ( dzeta_dy_u(k-1,j,i) + dzeta_dy_u(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1)                                               &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )      &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dy_u(k,j,i)                             &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                 &
                      + ( kmzp * ( ( u(k+1,j,i) - u(k,j,i) ) * ddzzu(k+1) * dzeta_dz_u(j,i)        &
                                 + ( w(k,j,i) - w(k,j,i-1) ) * ddx                                 &
                                 + ( ( w(k+1,j,i) + w(k+1,j,i-1) + w(k,j,i) + w(k,j,i-1) )         &
                                   - ( w(k-1,j,i) + w(k-1,j,i-1) + w(k,j,i) + w(k,j,i-1) )         &
                                   ) * 0.25_wp * ddzzu(k+1)                                        &
                                   * ( dzeta_dx_u(k+1,j,i) + dzeta_dx_u(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k)                                                 &
                        - kmzm * ( ( u(k,j,i) - u(k-1,j,i) ) * ddzzu(k) * dzeta_dz_u(j,i)          &
                                 + ( w(k-1,j,i) - w(k-1,j,i-1) ) * ddx                             &
                                 + ( ( w(k-1,j,i) + w(k-1,j,i-1) + w(k,j,i) + w(k,j,i-1) )         &
                                   - ( w(k-1,j,i) + w(k-1,j,i-1) + w(k-2,j,i) + w(k-2,j,i-1) )     &
                                   ) * 0.25_wp * ddzzu(k)                                          &
                                   * ( dzeta_dx_u(k-1,j,i) + dzeta_dx_u(k,j,i) ) * 0.5_wp          &
                                 ) * rho_air_zw(k-1)                                               &
                                   * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )      &
                        ) * ddzzw(k) * drho_air(k) * dzeta_dz_u(j,i)                               &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! div
                sums_l(k,pr_palm+9,tn) = sums_l(k,pr_palm+9,tn)                                    &
                      + ( ( u(k,j,i+1)   + u(k,j,i)   ) * 0.5_wp                                   &
                        - ( u(k,j,i-1)   + u(k,j,i)   ) * 0.5_wp                                   &
                        ) * ddx                                                                    &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                   &
                      + ( ( u(k+1,j,i)   + u(k,j,i)   ) * 0.5_wp                                   &
                        - ( u(k-1,j,i)   + u(k,j,i)   ) * 0.5_wp                                   &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )               &
                        ) * ddzzw(k) * dzeta_dx_u(k,j,i)                                           &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                   &
                      + ( ( v(k,j+1,i-1) + v(k,j+1,i) ) * 0.5_wp                                   &
                        - ( v(k,j,i-1)   + v(k,j,i)   ) * 0.5_wp                                   &
                        ) * ddy                                                                    &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                   &
                      + ( ( v(k+1,j+1,i-1) + v(k+1,j+1,i) + v(k+1,j,i-1) + v(k+1,j,i)              &
                          + v(k,j+1,i-1)   + v(k,j+1,i)   + v(k,j,i-1)   + v(k,j,i)   ) * 0.125_wp &
                        - ( v(k-1,j+1,i-1) + v(k-1,j+1,i) + v(k-1,j,i-1) + v(k-1,j,i)              &
                          + v(k,j+1,i-1)   + v(k,j+1,i)   + v(k,j,i-1)   + v(k,j,i)   ) * 0.125_wp &
                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k-1,j,i), 8 ) )               &
                        ) * ddzzw(k) * dzeta_dy_u(k,j,i)                                           &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )                   &
                      + ( ( w(k,j,i-1)   + w(k,j,i)   ) * 0.5_wp                                   &
                        - ( w(k-1,j,i-1) + w(k-1,j,i) ) * 0.5_wp                                   &
                        ) * ddzzw(k) * dzeta_dz_u(j,i)                                             &
                        * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                ! mflpz
                sums_l(k,pr_palm+10,tn) = sums_l(k,pr_palm+10,tn)                                  &
                                        + p(k,j,i) * drho_air(k) * dzeta_dx_s(k,j,i)

                IF ( k == nzb+1 ) THEN
                   sums_l(k,pr_palm+3,tn) = sums_l(k,pr_palm+3,tn)                                 &
                                          - u(k-1,j,i) * u(k-1,j,i)                                &
                                            * ddzzw(k) * dzeta_dx_u(k,j,i)                         &
                                           * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )  &
                                          - ( v(k-1,j,i)     + v(k-1,j+1,i)                        &
                                            + v(k-1,j+1,i-1) + v(k-1,j,i-1) ) *                    &
                                            u(k-1,j,i) * 0.25_wp                                   &
                                            * ddzzw(k) * dzeta_dy_u(k,j,i)                         &
                                           * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )  &
                                          - ( w(k-1,j,i-1) + w(k-1,j,i) ) *                        &
                                            u(k-1,j,i) * 0.5_wp                                    &
                                            * ddzzw(k) * dzeta_dz_u(j,i)                           &
                                           * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   sums_l(k,pr_palm+5,tn) = sums_l(k,pr_palm+5,tn)                                 &
                                          - ( p(k-1,j,i) * p(k-1,j,i-1) ) * 0.5_wp                 &
                                            * drho_air_zw(k-1)                                     &
                                            * ddzzw(k) * dzeta_dx_u(k,j,i)                         &
                                           * rmask(j,i,sr)                                         &
                                         * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   sums_l(k,pr_palm+6,tn) = sums_l(k,pr_palm+6,tn)                                 &
                         + ( - km(k,j,i)   * ( u(k-1,j,i) + u(k-1,j,i+1) )                         &
                                           * 0.5_wp * ddzzw(k) * dzeta_dx_s(k,j,i)                 &
                             + km(k,j,i-1) * ( u(k-1,j,i) + u(k-1,j,i-1) )                         &
                                           * 0.5_wp * ddzzw(k) * dzeta_dx_s(k,j,i-1)               &
                           ) * 2.0_wp * ddx                                                        &
                             * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   sums_l(k,pr_palm+7,tn) = sums_l(k,pr_palm+7,tn)                                 &
                      + ( - kmyp * ( u(k-1,j,i)   + u(k-1,j+1,i) )                                 &
                                 * 0.5_wp * ddzzw(k)                                               &
                                 * ( dzeta_dy_u(k,j,i) + dzeta_dy_u(k,j+1,i) ) * 0.5_wp            &
                          - kmyp * ( v(k-1,j+1,i) + v(k-1,j+1,i-1) )                               &
                                 * 0.5_wp * ddzzw(k)                                               &
                                 * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k,j+1,i) ) * 0.5_wp            &
                          + kmym * ( u(k-1,j,i)   + u(k-1,j-1,i) )                                 &
                                 * 0.5_wp * ddzzw(k)                                               &
                                 * ( dzeta_dy_u(k,j,i) + dzeta_dy_u(k,j-1,i) ) * 0.5_wp            &
                          + kmym * ( v(k-1,j,i)   + v(k-1,j,i-1) )                                 &
                                 * 0.5_wp * ddzzw(k)                                               &
                                 * ( dzeta_dx_u(k,j,i) + dzeta_dx_u(k,j-1,i) ) * 0.5_wp            &
                        ) * ddy * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )

                   sums_l(k,pr_palm+9,tn) = sums_l(k,pr_palm+9,tn)                                 &
                                          + ( - u(k-1,j,i) * ddzzw(k) * dzeta_dx_u(k,j,i)          &
                                              - ( v(k-1,j+1,i-1) + v(k-1,j+1,i)                    &
                                                + v(k-1,j,i-1)   + v(k-1,j,i)   ) * 0.25_wp        &
                                                * ddzzw(k) * dzeta_dy_u(k,j,i)                     &
                                            )                                                      &
                                          * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) ) 
                ENDIF
                ! EOC - Xu
             ENDDO

             IF ( use_surface_fluxes )  THEN
                surf_s = surf_def_h(0)%start_index(j,i)
                surf_e = surf_def_h(0)%end_index(j,i)
                DO  m = surf_s, surf_e
                   k = surf_def_h(0)%k(m)
                   sums_l(k,pr_palm+8,tn) = sums_l(k,pr_palm+8,tn)                                 &
                        + ( - ( - surf_def_h(0)%usws(m) ) * rho_air_zw(k-1)                        &
                              * ddzzw(k) * drho_air(k) * dzeta_dz_u(j,i)                           &
                          ) * MERGE( 1.0_wp, 0.0_wp, BTEST( topo_flags(k,j,i), 1 ) )
                ENDDO
             ENDIF
          ENDDO
       ENDDO


    ELSEIF ( mode == 'time_series' )  THEN


!       ALLOCATE ( ts_value_l(dots_num_user) )
!
!--    Sample on how to add values for the user-defined time series quantities.
!--    These have to be defined before in routine user_init. This sample creates two time series for
!--    the absolut values of the horizontal velocities u and v.
!       ts_value_l = 0.0_wp
!       ts_value_l(1) = ABS( u_max )
!       ts_value_l(2) = ABS( v_max )
!
!--     Collect / send values to PE0, because only PE0 outputs the time series.
!--     CAUTION: Collection is done by taking the sum over all processors. You may have to normalize
!--              this sum, depending on the quantity that you like to calculate. For serial runs,
!--              nothing has to be done.
!--     HINT: If the time series value that you are calculating has the same value on all PEs, you
!--           can omit the MPI_ALLREDUCE call and assign ts_value(dots_num_palm+1:,sr) = ts_value_l directly.
!#if defined( __parallel )
!       IF ( collective_wait )  CALL MPI_BARRIER( comm2d, ierr )
!       CALL MPI_ALLREDUCE( ts_value_l(1), ts_value(dots_num_palm+1,sr), dots_num_user, MPI_REAL,   &
!                           MPI_MAX, comm2d, ierr )
!#else
!       ts_value(dots_num_palm+1:dots_num_palm+dots_num_user,sr) = ts_value_l
!#endif

    ENDIF

 END SUBROUTINE user_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_global_ftn( found )


    LOGICAL, INTENT(OUT)  ::  found  !<


    found = .TRUE.


    SELECT CASE ( restart_string(1:length) )

       CASE ( 'global_paramter' )
!          READ ( 13 )  global_parameter

       CASE DEFAULT

          found = .FALSE.

    END SELECT


 END SUBROUTINE user_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific global restart data (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_global_mpi

!    USE restart_data_mpi_io_mod,                                                                   &
!        ONLY:  rrd_mpi_io

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )
    CONTINUE

 END SUBROUTINE user_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!> Subdomain
!> index limits on file are given by nxl_on_file, etc. Indices nxlc, etc. indicate the range of
!> gridpoints to be mapped from the subdomain on file (f) to the subdomain of the current PE (c).
!> They have been calculated in routine rrd_local.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_local_ftn( k, nxlf, nxlc, nxl_on_file, nxrf, nxrc, nxr_on_file, nynf, nync,   &
                                nyn_on_file, nysf, nysc, nys_on_file, tmp_3d, found )


    INTEGER(iwp) ::  idum            !<
    INTEGER(iwp) ::  k               !<
    INTEGER(iwp) ::  nxlc            !<
    INTEGER(iwp) ::  nxlf            !<
    INTEGER(iwp) ::  nxl_on_file     !<
    INTEGER(iwp) ::  nxrc            !<
    INTEGER(iwp) ::  nxrf            !<
    INTEGER(iwp) ::  nxr_on_file     !<
    INTEGER(iwp) ::  nync            !<
    INTEGER(iwp) ::  nynf            !<
    INTEGER(iwp) ::  nyn_on_file     !<
    INTEGER(iwp) ::  nysc            !<
    INTEGER(iwp) ::  nysf            !<
    INTEGER(iwp) ::  nys_on_file     !<

    LOGICAL, INTENT(OUT)  ::  found  !<

    REAL(wp), DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp) :: tmp_3d  !<

!
!-- Next line is to avoid compiler warning about unused variables. Please remove.
    idum = k + nxlc + nxlf + nxrc + nxrf + nync + nynf + nysc + nysf +                             &
           INT( tmp_3d(nzb,nys_on_file,nxl_on_file) )

!
!-- Here the reading of user-defined restart data follows:
!-- Sample for user-defined output

    found = .TRUE.

    SELECT CASE ( restart_string(1:length) )

!       CASE ( 'u2_av' )
!          IF ( .NOT. ALLOCATED( u2_av ) )  THEN
!               ALLOCATE( u2_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
!          ENDIF
!          IF ( k == 1 )  READ ( 13 )  tmp_3d
!             u2_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                    &
!             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
!
       ! BOC - Xu
       CASE ( 'wu_sgs' )
          IF ( .NOT. ALLOCATED( wu_sgs ) )  THEN
               ALLOCATE( wu_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wu_sgs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wu_sgs_av' )
          IF ( .NOT. ALLOCATED( wu_sgs_av ) )  THEN
               ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wu_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wv_sgs' )
          IF ( .NOT. ALLOCATED( wv_sgs ) )  THEN
               ALLOCATE( wv_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wv_sgs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                   &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wv_sgs_av' )
          IF ( .NOT. ALLOCATED( wv_sgs_av ) )  THEN
               ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wv_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wpt_sgs' )
          IF ( .NOT. ALLOCATED( wpt_sgs ) )  THEN
               ALLOCATE( wpt_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wpt_sgs(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                                 &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)

       CASE ( 'wpt_sgs_av' )
          IF ( .NOT. ALLOCATED( wpt_sgs_av ) )  THEN
               ALLOCATE( wpt_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
          ENDIF
          IF ( k == 1 )  READ ( 13 )  tmp_3d
             wpt_sgs_av(:,nysc-nbgp:nync+nbgp,nxlc-nbgp:nxrc+nbgp) =                              &
             tmp_3d(:,nysf-nbgp:nynf+nbgp,nxlf-nbgp:nxrf+nbgp)
       ! EOC - Xu

       CASE DEFAULT

          found = .FALSE.

    END SELECT

 END SUBROUTINE user_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_rrd_local_mpi

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  rd_mpi_io_check_array, rrd_mpi_io

    IMPLICIT NONE

    LOGICAL ::  array_found

    ! BOC - Xu
    CALL rd_mpi_io_check_array( 'wu_sgs' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wu_sgs ) )  ALLOCATE( wu_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wu_sgs', wu_sgs )
    ENDIF

    CALL rd_mpi_io_check_array( 'wu_sgs_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wu_sgs_av ) )  ALLOCATE( wu_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wu_sgs_av', wu_sgs_av )
    ENDIF

    CALL rd_mpi_io_check_array( 'wv_sgs' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wv_sgs ) )  ALLOCATE( wv_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wv_sgs', wv_sgs )
    ENDIF

    CALL rd_mpi_io_check_array( 'wv_sgs_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wv_sgs_av ) )  ALLOCATE( wv_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wv_sgs_av', wv_sgs_av )
    ENDIF

    CALL rd_mpi_io_check_array( 'wpt_sgs' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wpt_sgs ) )  ALLOCATE( wpt_sgs(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wpt_sgs', wpt_sgs )
    ENDIF

    CALL rd_mpi_io_check_array( 'wpt_sgs_av' , found = array_found )
    IF ( array_found )  THEN
       IF ( .NOT. ALLOCATED( wpt_sgs_av ) )  ALLOCATE( wpt_sgs_av(nzb:nzt+1,nysg:nyng,nxlg:nxrg) )
       CALL rrd_mpi_io( 'wpt_sgs_av', wpt_sgs_av )
    ENDIF
    ! EOC - Xu

    CONTINUE

 END SUBROUTINE user_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes global and user-defined restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_wrd_global

!    USE restart_data_mpi_io_mod,                                                                   &
!        ONLY:  wrd_mpi_io

    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

!       CALL wrd_write_string( 'global_parameter' )
!       WRITE ( 14 )  global_parameter

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!    CALL rrd_mpi_io( 'global_parameter', global_parameter )

    ENDIF

 END SUBROUTINE user_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Writes processor specific and user-defined restart data into binary file(s) for restart runs.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_wrd_local

    USE restart_data_mpi_io_mod,                                                                   &
        ONLY:  wrd_mpi_io

!
!-- Here the user-defined actions at the end of a job follow.
!-- Sample for user-defined output:

    ! BOC - Xu
    IF ( TRIM( restart_data_format_output ) == 'fortran_binary' )  THEN

!       IF ( ALLOCATED( u2_av ) )  THEN
!          CALL wrd_write_string( 'u2_av' )
!          WRITE ( 14 )  u2_av
!       ENDIF

       IF ( ALLOCATED( wu_sgs ) )  THEN
          CALL wrd_write_string( 'wu_sgs' )
          WRITE ( 14 )  wu_sgs
       ENDIF

       IF ( ALLOCATED( wu_sgs_av ) )  THEN
          CALL wrd_write_string( 'wu_sgs_av' )
          WRITE ( 14 )  wu_sgs_av
       ENDIF

       IF ( ALLOCATED( wv_sgs ) )  THEN
          CALL wrd_write_string( 'wv_sgs' )
          WRITE ( 14 )  wv_sgs
       ENDIF

       IF ( ALLOCATED( wv_sgs_av ) )  THEN
          CALL wrd_write_string( 'wv_sgs_av' )
          WRITE ( 14 )  wv_sgs_av
       ENDIF       

       IF ( ALLOCATED( wpt_sgs ) )  THEN
          CALL wrd_write_string( 'wpt_sgs' )
          WRITE ( 14 )  wpt_sgs
       ENDIF

       IF ( ALLOCATED( wpt_sgs_av ) )  THEN
          CALL wrd_write_string( 'wpt_sgs_av' )
          WRITE ( 14 )  wpt_sgs_av
       ENDIF

    ELSEIF ( restart_data_format_output(1:3) == 'mpi' )  THEN

!       IF ( ALLOCATED( u2_av ) )  CALL wrd_mpi_io( 'u2_av', u2_av )

       IF ( ALLOCATED( wu_sgs ) )       CALL wrd_mpi_io( 'wu_sgs', wu_sgs )
       IF ( ALLOCATED( wu_sgs_av ) )    CALL wrd_mpi_io( 'wu_sgs_av', wu_sgs_av )
       IF ( ALLOCATED( wv_sgs ) )       CALL wrd_mpi_io( 'wv_sgs', wv_sgs )
       IF ( ALLOCATED( wv_sgs_av ) )    CALL wrd_mpi_io( 'wv_sgs_av', wv_sgs_av )
       IF ( ALLOCATED( wpt_sgs ) )      CALL wrd_mpi_io( 'wpt_sgs', wpt_sgs )
       IF ( ALLOCATED( wpt_sgs_av ) )   CALL wrd_mpi_io( 'wpt_sgs_av', wpt_sgs_av )
    ENDIF
    ! EOC - Xu

 END SUBROUTINE user_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Execution of user-defined actions at the end of a job.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE user_last_actions

!
!-- Here the user-defined actions at the end of a job might follow.

    ! BOC - Xu
    !IF ( ALLOCATED( wu_sgs_av ) ) THEN
    !WRITE ( 14 ) 'wu_sgs_av'
    !WRITE ( 14 ) wu_sgs_av
    !ENDIF
    ! EOC - Xu

 END SUBROUTINE user_last_actions

 END MODULE user
