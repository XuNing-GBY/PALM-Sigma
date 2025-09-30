!> @file module_interface.f90
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
!> This is the interface between the PALM model core and all its modules.
!>
!> @todo Re-format module to be consistent with coding standard
!--------------------------------------------------------------------------------------------------!
 MODULE module_interface

    ! BOC - DEBUG
    USE arrays_3d,                                                                                 &
        ONLY:  w, w_1, w_2, w_3, u, u_1, u_2, u_3
    ! EOC - DEBUG

    USE indices,                                                                                   &
        ONLY:  nbgp, nxl, nxlg, nxr, nxrg, nys, nysg, nyn, nyng, nzb, nzt

    USE kinds

    USE pegrid,                                                                                    &
        ONLY:  comm2d

!
!-- load module-specific control parameters.
!-- ToDo: move all of them to respective module or a dedicated central module
    USE data_output_module,                                                                        &
        ONLY:  dom_def_end,                                                                        &
               dom_finalize_output,                                                                &
               dom_init

    USE data_output_topo_and_surface_setup_mod,                                                    &
        ONLY:  do_topo_and_surface_actions,                                                        &
               do_topo_and_surface_init_output

    USE dynamics_mod,                                                                              &
        ONLY:  dynamics_3d_data_averaging,                                                         &
               dynamics_actions,                                                                   &
               dynamics_boundary_conditions,                                                       &
               dynamics_check_data_output,                                                         &
               dynamics_check_data_output_pr,                                                      &
               dynamics_check_data_output_ts,                                                      &
               dynamics_check_parameters,                                                          &
               dynamics_data_output_2d,                                                            &
               dynamics_data_output_3d,                                                            &
               dynamics_define_netcdf_grid,                                                        &
               dynamics_exchange_horiz,                                                            &
               dynamics_header,                                                                    &
               dynamics_init,                                                                      &
               dynamics_init_arrays,                                                               &
               dynamics_init_checks,                                                               &
               dynamics_init_masks,                                                                &
               dynamics_last_actions,                                                              &
               dynamics_non_advective_processes,                                                   &
               dynamics_parin,                                                                     &
               dynamics_prognostic_equations,                                                      &
               dynamics_rrd_global,                                                                &
               dynamics_rrd_local,                                                                 &
               dynamics_statistics,                                                                &
               dynamics_swap_timelevel,                                                            &
               dynamics_wrd_global,                                                                &
               dynamics_wrd_local

    USE turbulence_closure_mod,                                                                    &
        ONLY:  tcm_3d_data_averaging,                                                              &
               tcm_actions,                                                                        &
               tcm_boundary_conds,                                                                 &
               tcm_check_data_output,                                                              &
               tcm_check_parameters,                                                               &
               tcm_data_output_2d,                                                                 &
               tcm_data_output_3d,                                                                 &
               tcm_init,                                                                           &
               tcm_init_arrays,                                                                    &
               tcm_prognostic_equations,                                                           &
               tcm_swap_timelevel

    ! BOC - Xu
    USE control_parameters,                                                                        &
        ONLY:  air_chemistry,                                                                      &
               biometeorology,                                                                     &
               coupling_char,                                                                      &
               debug_output,                                                                       &
               debug_output_timestep,                                                              &
               indoor_model,                                                                       &
               land_surface,                                                                       &
               large_scale_forcing,                                                                &
               nesting_offline,                                                                    &
               nudging,                                                                            &
               ocean_mode,                                                                         &
               plant_canopy,                                                                       &
               salsa,                                                                              &
               surface_output,                                                                     &
               syn_turb_gen,                                                                       &
               urban_surface,                                                                      &
               vdi_checks,                                                                         &
               virtual_flight,                                                                     &
               virtual_measurement,                                                                &
               wind_turbine,                                                                       &
               wind_wave_inta
    ! EOC - Xu

!
!-- Load interface routines of all PALM modules
    USE biometeorology_mod,                                                                        &
        ONLY:  bio_3d_data_averaging,                                                              &
               bio_check_data_output,                                                              &
               bio_data_output_2d,                                                                 &
               bio_data_output_3d,                                                                 &
               bio_init,                                                                           &
               bio_init_checks,                                                                    &
               bio_header,                                                                         &
               bio_parin,                                                                          &
               bio_rrd_global,                                                                     &
               bio_rrd_local,                                                                      &
               bio_wrd_global,                                                                     &
               bio_wrd_local

    USE bulk_cloud_model_mod,                                                                      &
        ONLY:  bcm_3d_data_averaging,                                                              &
               bcm_actions,                                                                        &
               bcm_boundary_conditions,                                                            &
               bcm_check_data_output_pr,                                                           &
               bcm_check_data_output,                                                              &
               bcm_check_parameters,                                                               &
               bcm_data_output_2d,                                                                 &
               bcm_data_output_3d,                                                                 &
               bcm_exchange_horiz,                                                                 &
               bcm_header,                                                                         &
               bcm_init_arrays,                                                                    &
               bcm_init,                                                                           &
               bcm_non_advective_processes,                                                        &
               bcm_parin,                                                                          &
               bcm_prognostic_equations,                                                           &
               bcm_rrd_global,                                                                     &
               bcm_rrd_local,                                                                      &
               bcm_swap_timelevel,                                                                 &
               bcm_wrd_global,                                                                     &
               bcm_wrd_local,                                                                      &
               bulk_cloud_model

   USE chemistry_model_mod,                                                                        &
       ONLY:  chem_3d_data_averaging,                                                              &
              chem_actions,                                                                        &
              chem_boundary_conditions,                                                            &
              chem_check_data_output,                                                              &
              chem_check_data_output_pr,                                                           &
              chem_check_parameters,                                                               &
              chem_data_output_2d,                                                                 &
              chem_data_output_3d,                                                                 &
              chem_exchange_horiz_bounds,                                                          &
              chem_header,                                                                         &
              chem_init,                                                                           &
              chem_init_arrays,                                                                    &
              chem_non_advective_processes,                                                        &
              chem_parin,                                                                          &
              chem_prognostic_equations,                                                           &
              chem_rrd_local,                                                                      &
              chem_statistics,                                                                     &
              chem_swap_timelevel,                                                                 &
              chem_wrd_local

    USE diagnostic_output_quantities_mod,                                                          &
        ONLY:  doq_3d_data_averaging,                                                              &
               doq_actions,                                                                        &
               doq_check_data_output,                                                              &
               doq_check_data_output_pr,                                                           &
               doq_define_netcdf_grid,                                                             &
               doq_init,                                                                           &
               doq_output_2d,                                                                      &
               doq_output_3d,                                                                      &
               doq_statistics,                                                                     &
               doq_rrd_local,                                                                      &
               doq_wrd_local

    USE flight_mod,                                                                                &
        ONLY:  flight_header,                                                                      &
               flight_init,                                                                        &
               flight_parin,                                                                       &
               flight_rrd_global,                                                                  &
               flight_wrd_global

    USE gust_mod,                                                                                  &
        ONLY:  gust_3d_data_averaging,                                                             &
               gust_actions,                                                                       &
               gust_check_parameters,                                                              &
               gust_check_data_output,                                                             &
               gust_check_data_output_pr,                                                          &
               gust_data_output_2d,                                                                &
               gust_data_output_3d,                                                                &
               gust_header,                                                                        &
               gust_init,                                                                          &
               gust_init_arrays,                                                                   &
               gust_module_enabled,                                                                &
               gust_parin,                                                                         &
               gust_prognostic_equations,                                                          &
               gust_rrd_global,                                                                    &
               gust_rrd_local,                                                                     &
               gust_statistics,                                                                    &
               gust_swap_timelevel,                                                                &
               gust_wrd_global,                                                                    &
               gust_wrd_local

    USE indoor_model_mod,                                                                          &
        ONLY:  im_check_data_output,                                                               &
               im_check_parameters,                                                                &
               im_data_output_3d,                                                                  &
               im_init,                                                                            &
               im_parin

    USE lagrangian_particle_model_mod,                                                             &
        ONLY:  lpm_actions,                                                                        &
               lpm_check_parameters,                                                               &
               lpm_exchange_horiz_bounds,                                                          &
               lpm_header,                                                                         &
               lpm_init,                                                                           &
               lpm_init_arrays,                                                                    &
               lpm_last_actions,                                                                   &
               lpm_parin,                                                                          &
               lpm_rrd_global,                                                                     &
               lpm_rrd_local,                                                                      &
               lpm_wrd_local,                                                                      &
               lpm_wrd_global

    USE land_surface_model_mod,                                                                    &
        ONLY:  lsm_3d_data_averaging,                                                              &
               lsm_check_parameters,                                                               &
               lsm_check_data_output,                                                              &
               lsm_check_data_output_pr,                                                           &
               lsm_data_output_2d,                                                                 &
               lsm_header,                                                                         &
               lsm_init,                                                                           &
               lsm_init_arrays,                                                                    &
               lsm_parin,                                                                          &
               lsm_rrd_local,                                                                      &
               lsm_swap_timelevel,                                                                 &
               lsm_wrd_local

    USE lsf_nudging_mod,                                                                           &
        ONLY:  lsf_nudging_check_parameters,                                                       &
               lsf_nudging_check_data_output_pr,                                                   &
               lsf_nudging_header,                                                                 &
               lsf_init,                                                                           &
               nudge_init

    USE multi_agent_system_mod,                                                                    &
        ONLY:  mas_parin

    USE nesting_offl_mod,                                                                          &
        ONLY:  nesting_offl_check_parameters,                                                      &
               nesting_offl_header,                                                                &
               nesting_offl_parin

    USE ocean_mod,                                                                                 &
        ONLY:  ocean_3d_data_averaging,                                                            &
               ocean_actions,                                                                      &
               ocean_boundary_conditions,                                                          &
               ocean_check_parameters,                                                             &
               ocean_check_data_output,                                                            &
               ocean_check_data_output_pr,                                                         &
               ocean_data_output_2d,                                                               &
               ocean_data_output_3d,                                                               &
               ocean_exchange_horiz,                                                               &
               ocean_init_arrays,                                                                  &
               ocean_init,                                                                         &
               ocean_header,                                                                       &
               ocean_parin,                                                                        &
               ocean_prognostic_equations,                                                         &
               ocean_rrd_global,                                                                   &
               ocean_rrd_local,                                                                    &
               ocean_swap_timelevel,                                                               &
               ocean_wrd_global,                                                                   &
               ocean_wrd_local

    USE particle_attributes,                                                                       &
        ONLY:  particle_advection

    USE poismg_noopt_mod,                                                                          &
        ONLY:  poismg_noopt_init

    USE plant_canopy_model_mod,                                                                    &
         ONLY: pcm_3d_data_averaging,                                                              &
               pcm_check_parameters,                                                               &
               pcm_check_data_output,                                                              &
               pcm_data_output_3d,                                                                 &
               pcm_header,                                                                         &
               pcm_init,                                                                           &
               pcm_parin,                                                                          &
               pcm_rrd_global,                                                                     &
               pcm_rrd_local,                                                                      &
               pcm_wrd_global,                                                                     &
               pcm_wrd_local

    USE radiation_model_mod,                                                                       &
        ONLY:  radiation,                                                                          &
               radiation_3d_data_averaging,                                                        &
               radiation_check_parameters,                                                         &
               radiation_check_data_output,                                                        &
               radiation_check_data_output_pr,                                                     &
               radiation_check_data_output_ts,                                                     &
               radiation_data_output_2d,                                                           &
               radiation_data_output_3d,                                                           &
               radiation_header,                                                                   &
               radiation_init,                                                                     &
               radiation_parin,                                                                    &
               radiation_rrd_local,                                                                &
               radiation_wrd_local

    USE salsa_mod,                                                                                 &
        ONLY:  salsa_3d_data_averaging,                                                            &
               salsa_actions,                                                                      &
               salsa_boundary_conditions,                                                          &
               salsa_check_data_output,                                                            &
               salsa_check_data_output_pr,                                                         &
               salsa_check_parameters,                                                             &
               salsa_data_output_2d,                                                               &
               salsa_data_output_3d,                                                               &
               salsa_exchange_horiz_bounds,                                                        &
               salsa_header,                                                                       &
               salsa_init,                                                                         &
               salsa_init_arrays,                                                                  &
               salsa_non_advective_processes,                                                      &
               salsa_parin,                                                                        &
               salsa_prognostic_equations,                                                         &
               salsa_rrd_global,                                                                   &
               salsa_rrd_local,                                                                    &
               salsa_statistics,                                                                   &
               salsa_swap_timelevel,                                                               &
               salsa_wrd_global,                                                                   &
               salsa_wrd_local

    USE spectra_mod,                                                                               &
        ONLY:  calculate_spectra,                                                                  &
               spectra_check_parameters,                                                           &
               spectra_header,                                                                     &
               spectra_parin

    USE surface_data_output_mod,                                                                   &
        ONLY:  surface_data_output_check_parameters,                                               &
               surface_data_output_init_arrays,                                                    &
               surface_data_output_parin,                                                          &
               surface_data_output_rrd_global,                                                     &
               surface_data_output_rrd_local,                                                      &
               surface_data_output_wrd_global,                                                     &
               surface_data_output_wrd_local

    USE surface_mod,                                                                               &
        ONLY:  init_bc

    USE synthetic_turbulence_generator_mod,                                                        &
        ONLY:  stg_check_parameters,                                                               &
               stg_header,                                                                         &
               stg_parin,                                                                          &
               stg_rrd_global,                                                                     &
               stg_wrd_global

    USE urban_surface_mod,                                                                         &
        ONLY:  usm_3d_data_averaging,                                                              &
               usm_check_data_output,                                                              &
               usm_check_parameters,                                                               &
               usm_init,                                                                           &
               usm_init_arrays,                                                                    &
               usm_parin,                                                                          &
               usm_swap_timelevel,                                                                 &
               usm_rrd_local,                                                                      &
               usm_wrd_local

    USE vdi_internal_controls,                                                                     &
        ONLY:  vdi_actions

    USE virtual_measurement_mod,                                                                   &
        ONLY:  vm_check_parameters,                                                                &
               vm_init,                                                                            &
               vm_init_output,                                                                     &
               vm_parin

    USE wind_turbine_model_mod,                                                                    &
        ONLY:  wtm_actions,                                                                        &
               wtm_check_parameters,                                                               &
               wtm_init,                                                                           &
               wtm_init_arrays,                                                                    &
               wtm_init_output,                                                                    &
               wtm_parin,                                                                          &
               wtm_rrd_global,                                                                     &
               wtm_wrd_global

    ! BOC - Xu
    USE wind_wave_inta_mod,                                                                        &
        ONLY: wwinta_parin,                                                                        &
              wwinta_init,                                                                         &
              wwinta_exchange_horiz,                                                               &
              wwinta_prognostic_equations,                                                         &
              wwinta_actions,                                                                      &
              wwinta_boundary_conditions
    ! EOC - Xu

    USE user,                                                                                      &
        ONLY:  user_3d_data_averaging,                                                             &
               user_actions,                                                                       &
               user_check_data_output,                                                             &
               user_check_data_output_pr,                                                          &
               user_check_data_output_ts,                                                          &
               user_check_parameters,                                                              &
               user_data_output_2d,                                                                &
               user_data_output_3d,                                                                &
               user_header,                                                                        &
               user_init,                                                                          &
               user_init_arrays,                                                                   &
               user_last_actions,                                                                  &
               user_module_enabled,                                                                &
               user_parin,                                                                         &
               user_statistics,                                                                    &
               user_rrd_global,                                                                    &
               user_rrd_local,                                                                     &
               user_wrd_global,                                                                    &
               user_wrd_local

    IMPLICIT NONE

    PRIVATE

!
!-- Public functions
    PUBLIC                                                                                         &
       module_interface_parin,                                                                     &
       module_interface_check_parameters,                                                          &
       module_interface_check_data_output_ts,                                                      &
       module_interface_check_data_output_pr,                                                      &
       module_interface_check_data_output,                                                         &
       module_interface_init_masks,                                                                &
       module_interface_define_netcdf_grid,                                                        &
       module_interface_init_arrays,                                                               &
       module_interface_init,                                                                      &
       module_interface_init_checks,                                                               &
       module_interface_init_numerics,                                                             &
       module_interface_init_output,                                                               &
       module_interface_header,                                                                    &
       module_interface_actions,                                                                   &
       module_interface_non_advective_processes,                                                   &
       module_interface_exchange_horiz,                                                            &
       module_interface_prognostic_equations,                                                      &
       module_interface_boundary_conditions,                                                       &
       module_interface_swap_timelevel,                                                            &
       module_interface_3d_data_averaging,                                                         &
       module_interface_data_output_2d,                                                            &
       module_interface_data_output_3d,                                                            &
       module_interface_statistics,                                                                &
       module_interface_rrd_global,                                                                &
       module_interface_wrd_global,                                                                &
       module_interface_rrd_local,                                                                 &
       module_interface_wrd_local,                                                                 &
       module_interface_last_actions


    INTERFACE module_interface_parin
       MODULE PROCEDURE module_interface_parin
    END INTERFACE module_interface_parin

    INTERFACE module_interface_check_parameters
       MODULE PROCEDURE module_interface_check_parameters
    END INTERFACE module_interface_check_parameters

    INTERFACE module_interface_check_data_output_ts
       MODULE PROCEDURE module_interface_check_data_output_ts
    END INTERFACE module_interface_check_data_output_ts

    INTERFACE module_interface_check_data_output_pr
       MODULE PROCEDURE module_interface_check_data_output_pr
    END INTERFACE module_interface_check_data_output_pr

    INTERFACE module_interface_check_data_output
       MODULE PROCEDURE module_interface_check_data_output
    END INTERFACE module_interface_check_data_output

    INTERFACE module_interface_init_masks
       MODULE PROCEDURE module_interface_init_masks
    END INTERFACE module_interface_init_masks

    INTERFACE module_interface_define_netcdf_grid
       MODULE PROCEDURE module_interface_define_netcdf_grid
    END INTERFACE module_interface_define_netcdf_grid

    INTERFACE module_interface_init_arrays
       MODULE PROCEDURE module_interface_init_arrays
    END INTERFACE module_interface_init_arrays

    INTERFACE module_interface_init
       MODULE PROCEDURE module_interface_init
    END INTERFACE module_interface_init

    INTERFACE module_interface_init_checks
       MODULE PROCEDURE module_interface_init_checks
    END INTERFACE module_interface_init_checks

    INTERFACE module_interface_init_numerics
       MODULE PROCEDURE module_interface_init_numerics
    END INTERFACE module_interface_init_numerics

    INTERFACE module_interface_init_output
       MODULE PROCEDURE module_interface_init_output
    END INTERFACE module_interface_init_output

    INTERFACE module_interface_header
       MODULE PROCEDURE module_interface_header
    END INTERFACE module_interface_header

    INTERFACE module_interface_actions
       MODULE PROCEDURE module_interface_actions
       MODULE PROCEDURE module_interface_actions_ij
    END INTERFACE module_interface_actions

    INTERFACE module_interface_non_advective_processes
       MODULE PROCEDURE module_interface_non_advective_processes
       MODULE PROCEDURE module_interface_non_advective_processes_ij
    END INTERFACE module_interface_non_advective_processes

    INTERFACE module_interface_exchange_horiz
       MODULE PROCEDURE module_interface_exchange_horiz
    END INTERFACE module_interface_exchange_horiz

    INTERFACE module_interface_prognostic_equations
       MODULE PROCEDURE module_interface_prognostic_equations
       MODULE PROCEDURE module_interface_prognostic_equations_ij
    END INTERFACE module_interface_prognostic_equations

    INTERFACE module_interface_swap_timelevel
       MODULE PROCEDURE module_interface_swap_timelevel
    END INTERFACE module_interface_swap_timelevel

    INTERFACE module_interface_boundary_conditions
       MODULE PROCEDURE module_interface_boundary_conditions
    END INTERFACE module_interface_boundary_conditions

    INTERFACE module_interface_3d_data_averaging
       MODULE PROCEDURE module_interface_3d_data_averaging
    END INTERFACE module_interface_3d_data_averaging

    INTERFACE module_interface_data_output_2d
       MODULE PROCEDURE module_interface_data_output_2d
    END INTERFACE module_interface_data_output_2d

    INTERFACE module_interface_data_output_3d
       MODULE PROCEDURE module_interface_data_output_3d
    END INTERFACE module_interface_data_output_3d

    INTERFACE module_interface_statistics
       MODULE PROCEDURE module_interface_statistics
    END INTERFACE module_interface_statistics

    INTERFACE module_interface_rrd_global
       MODULE PROCEDURE module_interface_rrd_global_ftn
       MODULE PROCEDURE module_interface_rrd_global_mpi
    END INTERFACE module_interface_rrd_global

    INTERFACE module_interface_wrd_global
       MODULE PROCEDURE module_interface_wrd_global
    END INTERFACE module_interface_wrd_global

    INTERFACE module_interface_rrd_local
       MODULE PROCEDURE module_interface_rrd_local_ftn
       MODULE PROCEDURE module_interface_rrd_local_mpi
    END INTERFACE module_interface_rrd_local

    INTERFACE module_interface_wrd_local
       MODULE PROCEDURE module_interface_wrd_local
    END INTERFACE module_interface_wrd_local

    INTERFACE module_interface_last_actions
       MODULE PROCEDURE module_interface_last_actions
    END INTERFACE module_interface_last_actions


 CONTAINS


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific parameter namelists
!------------------------------------------------------------------------------                    !
 SUBROUTINE module_interface_parin


    IF ( debug_output )  CALL debug_message( 'reading module-specific parameters', 'start' )

    CALL dynamics_parin

    CALL bio_parin
    CALL bcm_parin
    CALL chem_parin
    CALL flight_parin ! ToDo: rename module to match filename
    CALL gust_parin
    CALL im_parin
    CALL lpm_parin
    CALL lsm_parin
    ! ToDo: create parin routine for large_scale_forcing and nudging (should be seperate modules or new module switch)
    CALL mas_parin
    CALL nesting_offl_parin
    CALL ocean_parin
    CALL pcm_parin
    CALL radiation_parin
    CALL salsa_parin
    CALL spectra_parin
    CALL surface_data_output_parin
    CALL stg_parin
    CALL usm_parin
    CALL vm_parin
    CALL wtm_parin

    ! BOC - Xu
    CALL wwinta_parin
    ! EOC - Xu

    CALL user_parin

    IF ( debug_output )  CALL debug_message( 'reading module-specific parameters', 'end' )


 END SUBROUTINE module_interface_parin


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific initialization checks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_check_parameters


    IF ( debug_output )  CALL debug_message( 'checking module-specific parameters', 'start' )

    CALL dynamics_check_parameters
    CALL tcm_check_parameters

    IF ( bulk_cloud_model )     CALL bcm_check_parameters
    IF ( air_chemistry )        CALL chem_check_parameters
    IF ( gust_module_enabled )  CALL gust_check_parameters
    IF ( indoor_model )         CALL im_check_parameters
    IF ( particle_advection )   CALL lpm_check_parameters
    IF ( land_surface )         CALL lsm_check_parameters
    IF ( large_scale_forcing  .OR.  nudging )  CALL lsf_nudging_check_parameters ! ToDo: create single module switch
    IF ( nesting_offline )      CALL nesting_offl_check_parameters
    IF ( ocean_mode )           CALL ocean_check_parameters
    IF ( plant_canopy )         CALL pcm_check_parameters
    IF ( radiation )            CALL radiation_check_parameters
    IF ( salsa )                CALL salsa_check_parameters
    IF ( calculate_spectra )    CALL spectra_check_parameters
    IF ( surface_output )       CALL surface_data_output_check_parameters
    IF ( syn_turb_gen )         CALL stg_check_parameters
    IF ( urban_surface )        CALL usm_check_parameters
    IF ( virtual_measurement )  CALL vm_check_parameters
    IF ( wind_turbine )         CALL wtm_check_parameters

    IF ( user_module_enabled )  CALL user_check_parameters

    IF ( debug_output )  CALL debug_message( 'checking module-specific parameters', 'end' )


 END SUBROUTINE module_interface_check_parameters


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check module-specific data output of timeseries
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    INTEGER(iwp),      INTENT(IN)       ::  dots_max !< variable output array index
    INTEGER(iwp),      INTENT(INOUT)    ::  dots_num !< variable output array index

    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT) :: dots_label
    CHARACTER (LEN=*), DIMENSION(dots_max), INTENT(INOUT) :: dots_unit


    IF ( debug_output )  CALL debug_message( 'checking module-specific data output ts', 'start' )

    CALL dynamics_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )

    IF ( radiation )  THEN
       CALL radiation_check_data_output_ts( dots_max, dots_num )
    ENDIF

    IF ( user_module_enabled )  THEN
       CALL user_check_data_output_ts( dots_max, dots_num, dots_label, dots_unit )
    ENDIF

    IF ( debug_output )  CALL debug_message( 'checking module-specific data output ts', 'end' )


 END SUBROUTINE module_interface_check_data_output_ts


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check module-specific data output of profiles
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_check_data_output_pr( variable, var_count, unit, dopr_unit )


    CHARACTER (LEN=*), INTENT(IN)    ::  variable  !< variable name
    CHARACTER (LEN=*), INTENT(INOUT) ::  unit      !< physical unit of variable
    CHARACTER (LEN=*), INTENT(OUT)   ::  dopr_unit !< local value of dopr_unit

    INTEGER(iwp),      INTENT(IN)    ::  var_count !< variable output array index


    IF ( debug_output )  CALL debug_message( 'checking module-specific data output pr', 'start' )

    CALL dynamics_check_data_output_pr( variable, var_count, unit, dopr_unit )
    CALL doq_check_data_output_pr( variable, var_count, unit, dopr_unit )

    IF ( unit == 'illegal'  .AND.  bulk_cloud_model )  THEN
       CALL bcm_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  air_chemistry )  THEN
       CALL chem_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  gust_module_enabled  )  THEN
       CALL gust_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal' )  THEN ! ToDo: add module switch if possible
       CALL lsm_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal' )  THEN ! ToDo: add module switch if possible
       CALL lsf_nudging_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  ocean_mode )  THEN
       CALL ocean_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  radiation )  THEN
       CALL radiation_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  salsa )  THEN
       CALL salsa_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  user_module_enabled )  THEN
       unit = '' ! ToDo: Seems like a hack. Find a general soultion!
       CALL user_check_data_output_pr( variable, var_count, unit, dopr_unit )
    ENDIF

    IF ( debug_output )  CALL debug_message( 'checking module-specific data output pr', 'end' )


 END SUBROUTINE module_interface_check_data_output_pr

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Check module-specific 2D and 3D data output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_check_data_output( variable, unit, i, j, ilen, k )


    CHARACTER (LEN=*), INTENT(IN)    ::  variable !< variable name
    CHARACTER (LEN=*), INTENT(INOUT) ::  unit     !< physical unit of variable

    INTEGER(iwp),      INTENT(IN)    :: i         !< ToDo: remove dummy argument, instead pass string from data_output
    INTEGER(iwp),      INTENT(IN)    :: j         !< average quantity? 0 = no, 1 = yes
    INTEGER(iwp),      INTENT(IN)    :: ilen      !< ToDo: remove dummy argument, instead pass string from data_output
    INTEGER(iwp),      INTENT(IN)    :: k         !< ToDo: remove dummy argument, instead pass string from data_output


    IF ( debug_output )  CALL debug_message( 'checking module-specific data output 2d/3d', 'start' )

    CALL dynamics_check_data_output( variable, unit )

    CALL tcm_check_data_output( variable, unit )

    IF ( unit == 'illegal'  .AND.  biometeorology )  THEN
       CALL bio_check_data_output( variable, unit, i, j, ilen, k )
    ENDIF

    IF ( unit == 'illegal'  .AND.  bulk_cloud_model  )  THEN
       CALL bcm_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  air_chemistry  .AND.  (variable(1:3) == 'kc_'  .OR.             &
         variable(1:3) == 'em_') )  THEN  ! ToDo: remove aditional conditions
       CALL chem_check_data_output( variable, unit, i, ilen, k )
    ENDIF

    IF ( unit == 'illegal' )  THEN
       CALL doq_check_data_output( variable, unit, i, ilen, k )
    ENDIF

    IF ( unit == 'illegal'  .AND.  gust_module_enabled  )  THEN
       CALL gust_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal' )  THEN  ! ToDo: add module switch if possible
       CALL lsm_check_data_output( variable, unit, i, ilen, k )
    ENDIF

    IF ( unit == 'illegal'  .AND.  ocean_mode )  THEN
       CALL ocean_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  plant_canopy  .AND.  variable(1:4) == 'pcm_' )  THEN  ! ToDo: remove aditional conditions
       CALL pcm_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  radiation )  THEN
       CALL radiation_check_data_output( variable, unit, i, ilen, k )
    ENDIF

    IF ( unit == 'illegal' .AND. salsa ) THEN
       CALL salsa_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal' .AND. indoor_model ) THEN
       CALL im_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  urban_surface  .AND.  variable(1:4) == 'usm_' )  THEN  ! ToDo: remove aditional conditions
       CALL usm_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  user_module_enabled )  THEN
       unit = ''
       CALL user_check_data_output( variable, unit )
    ENDIF

    IF ( debug_output )  CALL debug_message( 'checking module-specific data output 2d/3d', 'end' )


 END SUBROUTINE module_interface_check_data_output


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Interface for init_masks. ToDo: get rid of these redundant calls!
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init_masks( variable, unit )


    CHARACTER (LEN=*), INTENT(IN)    ::  variable !< variable name
    CHARACTER (LEN=*), INTENT(INOUT) ::  unit     !< physical unit of variable


    IF ( debug_output )  CALL debug_message( 'initializing module-specific masks', 'start' )

    CALL dynamics_init_masks( variable, unit )

    IF ( unit == 'illegal'  .AND.  air_chemistry  .AND.  (variable(1:3) == 'kc_'  .OR.             &
         variable(1:3) == 'em_') )  THEN  ! ToDo: remove aditional conditions
       CALL chem_check_data_output( variable, unit, 0, 0, 0 )
    ENDIF

    IF ( unit == 'illegal' )  THEN
       CALL doq_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  radiation )  THEN
       CALL radiation_check_data_output( variable, unit, 0, 0, 0 )
    ENDIF

    IF ( unit == 'illegal'  .AND.  salsa )  THEN
       CALL salsa_check_data_output( variable, unit )
    ENDIF

    IF ( unit == 'illegal'  .AND.  user_module_enabled )  THEN
       unit = ''
       CALL user_check_data_output( variable, unit )
    ENDIF

    IF ( debug_output )  CALL debug_message( 'initializing module-specific masks', 'end' )


 END SUBROUTINE module_interface_init_masks


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Define appropriate grid for module-specific netcdf output variables.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_define_netcdf_grid( var, found, grid_x, grid_y, grid_z )


    CHARACTER (LEN=*), INTENT(IN)  ::  var    !< variable name
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_x !< netcdf dimension in x-direction
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_y !< netcdf dimension in y-direction
    CHARACTER (LEN=*), INTENT(OUT) ::  grid_z !< netcdf dimension in z-direction

    LOGICAL,           INTENT(OUT) ::  found  !< indicates if variable was found


    IF ( debug_output )  CALL debug_message( 'defining module-specific netcdf grids', 'start' )
!
!-- As long as no action is done in this subroutine, initialize strings with intent(out) attribute,
!-- in order to avoid compiler warnings.
    found  = .FALSE.
    grid_x = 'none'
    grid_y = 'none'
    grid_z = 'none'
!
!-- Use var to avoid compiler warning about unused variable
    IF ( var == ' ' )  RETURN

    IF ( debug_output )  CALL debug_message( 'defining module-specific netcdf grids', 'end' )


 END SUBROUTINE module_interface_define_netcdf_grid


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Allocate module-specific arrays and pointers
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init_arrays


    IF ( debug_output )  CALL debug_message( 'initializing module-specific arrays', 'start' )

    CALL dynamics_init_arrays
    CALL tcm_init_arrays

    IF ( bulk_cloud_model    )  CALL bcm_init_arrays
    IF ( air_chemistry       )  CALL chem_init_arrays
    IF ( gust_module_enabled )  CALL gust_init_arrays
    IF ( particle_advection  )  CALL lpm_init_arrays
    IF ( land_surface        )  CALL lsm_init_arrays
    IF ( ocean_mode          )  CALL ocean_init_arrays
    IF ( salsa               )  CALL salsa_init_arrays
    IF ( urban_surface       )  CALL usm_init_arrays
    IF ( surface_output      )  CALL surface_data_output_init_arrays
    IF ( wind_turbine        )  CALL wtm_init_arrays

    IF ( user_module_enabled )  CALL user_init_arrays

    IF ( debug_output )  CALL debug_message( 'initializing module-specific arrays', 'end' )


 END SUBROUTINE module_interface_init_arrays


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific initialization
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init


    IF ( debug_output )  CALL debug_message( 'module-specific initialization', 'start' )

    ! BOC - Xu
    IF ( wind_wave_inta ) CALL wwinta_init
    ! EOC - Xu

    CALL dynamics_init
    CALL tcm_init

    IF ( biometeorology      )  CALL bio_init
    IF ( bulk_cloud_model    )  CALL bcm_init
    IF ( air_chemistry       )  CALL chem_init
    IF ( virtual_flight      )  CALL flight_init
    IF ( gust_module_enabled )  CALL gust_init
    IF ( indoor_model        )  CALL im_init
    IF ( particle_advection  )  CALL lpm_init
    IF ( large_scale_forcing )  CALL lsf_init
    IF ( land_surface        )  CALL lsm_init
    IF ( nudging             )  CALL nudge_init
    IF ( ocean_mode          )  CALL ocean_init
    IF ( plant_canopy        )  CALL pcm_init
    IF ( salsa               )  CALL salsa_init
    IF ( urban_surface       )  CALL usm_init
    IF ( virtual_measurement )  CALL vm_init
    IF ( wind_turbine        )  CALL wtm_init
    IF ( radiation           )  CALL radiation_init

    CALL doq_init

    IF ( user_module_enabled )  CALL user_init

    IF ( debug_output )  CALL debug_message( 'module-specific initialization', 'end' )

 END SUBROUTINE module_interface_init

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize boundary conditions and numerical schemes.
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init_numerics

!
!-- Initialize boundary conditions via surface type
    CALL init_bc
!
!-- Calculate wall flag arrays for the multigrid method.
!-- Please note, wall flags are only applied in the non-optimized version.
    CALL poismg_noopt_init

 END SUBROUTINE module_interface_init_numerics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Initialize data output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init_output

    INTEGER(iwp) ::  return_value  !< returned status value of called function

!
!-- Initialize data-output module
    CALL dom_init( file_suffix_of_output_group=coupling_char,                                      &
                   mpi_comm_of_output_group=comm2d,                                                &
                   program_debug_output_unit=9,                                                    &
                   debug_output=debug_output )
!
!-- Initialze output of topography and surface setup
    CALL do_topo_and_surface_init_output
!
!-- Define module-specific output quantities
    IF ( virtual_measurement )  CALL vm_init_output
    IF ( wind_turbine )         CALL wtm_init_output
!
!-- Leave output-definition state
    return_value = dom_def_end()

 END SUBROUTINE module_interface_init_output

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific post-initialization checks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_init_checks


    IF ( debug_output )  CALL debug_message( 'module-specific post-initialization checks', 'start' )

    CALL dynamics_init_checks

    IF ( biometeorology )  CALL bio_init_checks

    IF ( debug_output )  CALL debug_message( 'module-specific post-initialization checks', 'end' )


 END SUBROUTINE module_interface_init_checks


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Gather module-specific header output
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_header( io )


    INTEGER(iwp), INTENT(IN) ::  io  !< unit of the output file


    IF ( debug_output )  CALL debug_message( 'module-specific header output', 'start' )

    CALL dynamics_header( io )

    IF ( biometeorology      )  CALL bio_header ( io )
    IF ( bulk_cloud_model    )  CALL bcm_header( io )
    IF ( air_chemistry       )  CALL chem_header ( io )
    IF ( virtual_flight      )  CALL flight_header( io )
    IF ( gust_module_enabled )  CALL gust_header( io )
    IF ( particle_advection  )  CALL lpm_header( io )
    IF ( land_surface        )  CALL lsm_header( io )
    IF ( large_scale_forcing )  CALL lsf_nudging_header( io )
    IF ( nesting_offline     )  CALL nesting_offl_header( io )
    IF ( ocean_mode          )  CALL ocean_header( io )
    IF ( plant_canopy        )  CALL pcm_header( io )
    IF ( radiation           )  CALL radiation_header( io )
    IF ( salsa               )  CALL salsa_header( io )
    IF ( calculate_spectra   )  CALL spectra_header( io )
    IF ( syn_turb_gen        )  CALL stg_header( io )

    IF ( user_module_enabled )  CALL user_header( io )

    IF ( debug_output )  CALL debug_message( 'module-specific header output', 'end' )


 END SUBROUTINE module_interface_header


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific actions while in time-integration (vector-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_actions( location )


    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string

    CALL do_topo_and_surface_actions( location )
    CALL dynamics_actions( location )
    CALL tcm_actions( location )
    CALL doq_actions( location )

    ! BOC - Xu
    CALL wwinta_actions( location )
    ! EOC - Xu

    IF ( bulk_cloud_model    )  CALL bcm_actions( location )
    IF ( air_chemistry       )  CALL chem_actions( location )
    IF ( gust_module_enabled )  CALL gust_actions( location )
    IF ( particle_advection  )  CALL lpm_actions( location )
    IF ( ocean_mode          )  CALL ocean_actions( location )
    IF ( salsa               )  CALL salsa_actions( location )
    IF ( wind_turbine        )  CALL wtm_actions( location )

    IF ( user_module_enabled )  CALL user_actions( location )
    IF ( vdi_checks          )  CALL vdi_actions( location )


 END SUBROUTINE module_interface_actions


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific actions while in time-integration (cache-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_actions_ij( i, j, location )

    CHARACTER (LEN=*), INTENT(IN) ::  location  !< call location string

    INTEGER(iwp),      INTENT(IN) ::  i         !< grid index in x-direction
    INTEGER(iwp),      INTENT(IN) ::  j         !< grid index in y-direction

    CALL dynamics_actions( i, j, location )
    CALL tcm_actions( i, j, location )
    CALL doq_actions( i, j, location )

    ! BOC - Xu
    IF ( wind_wave_inta      )  CALL wwinta_actions( i, j, location )
    ! EOC - Xu

    IF ( bulk_cloud_model    )  CALL bcm_actions( i, j, location )
    IF ( air_chemistry       )  CALL chem_actions( i, j, location )
    IF ( gust_module_enabled )  CALL gust_actions( i, j, location )
    IF ( ocean_mode          )  CALL ocean_actions( i, j, location )
    IF ( salsa               )  CALL salsa_actions( i, j, location )
    IF ( wind_turbine        )  CALL wtm_actions( i, j, location )

    IF ( user_module_enabled )  CALL user_actions( i, j, location )


 END SUBROUTINE module_interface_actions_ij


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non_advective_processes (vector-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_non_advective_processes


    CALL dynamics_non_advective_processes

    IF ( bulk_cloud_model    )  CALL bcm_non_advective_processes
    IF ( air_chemistry       )  CALL chem_non_advective_processes
    IF ( salsa               )  CALL salsa_non_advective_processes


 END SUBROUTINE module_interface_non_advective_processes


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific non_advective_processes (cache-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_non_advective_processes_ij( i, j )


    INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction

    CALL dynamics_non_advective_processes( i, j )

    IF ( bulk_cloud_model    )  CALL bcm_non_advective_processes( i, j )
    IF ( air_chemistry       )  CALL chem_non_advective_processes( i, j )
    IF ( salsa               )  CALL salsa_non_advective_processes( i, j )


 END SUBROUTINE module_interface_non_advective_processes_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Exchange horiz for module-specific quantities
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_exchange_horiz( location )

    CHARACTER (LEN=*), INTENT(IN) ::  location !< call location string


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific exchange_horiz', 'start' )

    ! BOC - Xu
    IF ( wind_wave_inta      )  CALL wwinta_exchange_horiz( location )
    ! EOC - Xu

    CALL dynamics_exchange_horiz( location )

    IF ( bulk_cloud_model    )  CALL bcm_exchange_horiz( location )
    IF ( air_chemistry       )  CALL chem_exchange_horiz_bounds( location )
    IF ( ocean_mode          )  CALL ocean_exchange_horiz( location )
    IF ( particle_advection  )  CALL lpm_exchange_horiz_bounds( location )
    IF ( salsa               )  CALL salsa_exchange_horiz_bounds( location )

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific exchange_horiz', 'end' )


 END SUBROUTINE module_interface_exchange_horiz


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic_equations (vector-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_prognostic_equations


    CALL dynamics_prognostic_equations
    CALL tcm_prognostic_equations

    ! BOC - Xu
    ! to do
    ! EOC - Xu

    IF ( bulk_cloud_model    )  CALL bcm_prognostic_equations
    IF ( air_chemistry       )  CALL chem_prognostic_equations
    IF ( gust_module_enabled )  CALL gust_prognostic_equations
    IF ( ocean_mode          )  CALL ocean_prognostic_equations
    IF ( salsa               )  CALL salsa_prognostic_equations


 END SUBROUTINE module_interface_prognostic_equations


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific prognostic_equations (cache-optimized)
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_prognostic_equations_ij( i, j, i_omp_start, tn )


    INTEGER(iwp), INTENT(IN) ::  i            !< grid index in x-direction
    INTEGER(iwp), INTENT(IN) ::  j            !< grid index in y-direction
    INTEGER(iwp), INTENT(IN) ::  i_omp_start  !< first loop index of i-loop in prognostic_equations
    INTEGER(iwp), INTENT(IN) ::  tn           !< task number of openmp task

    CALL dynamics_prognostic_equations( i, j, i_omp_start, tn )
    CALL tcm_prognostic_equations( i, j, i_omp_start, tn )

    ! BOC - Xu
    CALL wwinta_prognostic_equations( i, j )
    ! EOC - Xu

    IF ( bulk_cloud_model    )  CALL bcm_prognostic_equations( i, j, i_omp_start, tn )
    IF ( air_chemistry       )  CALL chem_prognostic_equations( i, j, i_omp_start, tn )
    IF ( gust_module_enabled )  CALL gust_prognostic_equations( i, j, i_omp_start, tn )
    IF ( ocean_mode          )  CALL ocean_prognostic_equations( i, j, i_omp_start, tn )
    IF ( salsa               )  CALL salsa_prognostic_equations( i, j, i_omp_start, tn )


 END SUBROUTINE module_interface_prognostic_equations_ij

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific boundary conditions
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_boundary_conditions


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific boundary_conditions', 'start' )

    CALL dynamics_boundary_conditions
    CALL tcm_boundary_conds

    ! BOC - Xu
    IF ( wind_wave_inta      )  CALL wwinta_boundary_conditions 
    ! EOC - Xu

    IF ( bulk_cloud_model    )  CALL bcm_boundary_conditions
    IF ( air_chemistry       )  CALL chem_boundary_conditions
    IF ( ocean_mode          )  CALL ocean_boundary_conditions
    IF ( salsa               )  CALL salsa_boundary_conditions

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific boundary_conditions', 'end' )


 END SUBROUTINE module_interface_boundary_conditions

!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Swap the timelevel pointers for module-specific arrays
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_swap_timelevel ( swap_mode )


    INTEGER(iwp), INTENT(IN) :: swap_mode !< determines procedure of pointer swap


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific swap timelevel', 'start' )

    ! BOC - DEBUG
    !WRITE(*,*) "w_bf_dyswap: ", w(nzb,5,5), w_1(nzb,5,5), w_2(nzb,5,5), w_3(nzb,5,5)
    !WRITE(*,*) "u_bf_dyswap: ", u(nzb,5,5), u_1(nzb,5,5), u_2(nzb,5,5), u_3(nzb,5,5)
    CALL dynamics_swap_timelevel( swap_mode )
    !WRITE(*,*) "w_aft_dyswap: ", w(nzb,5,5), w_1(nzb,5,5), w_2(nzb,5,5), w_3(nzb,5,5)
    !WRITE(*,*) "u_aft_dyswap: ", u(nzb,5,5), u_1(nzb,5,5), u_2(nzb,5,5), u_3(nzb,5,5)
    ! EOC - DEBUG
    CALL tcm_swap_timelevel( swap_mode )

    IF ( bulk_cloud_model    )  CALL bcm_swap_timelevel( swap_mode )
    IF ( air_chemistry       )  CALL chem_swap_timelevel( swap_mode )
    IF ( gust_module_enabled )  CALL gust_swap_timelevel( swap_mode )
    IF ( land_surface        )  CALL lsm_swap_timelevel( swap_mode )
    IF ( ocean_mode          )  CALL ocean_swap_timelevel( swap_mode )
    IF ( salsa               )  CALL salsa_swap_timelevel( swap_mode )
    IF ( urban_surface       )  CALL usm_swap_timelevel( swap_mode )

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific swap timelevel', 'end' )


 END SUBROUTINE module_interface_swap_timelevel


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Perform module-specific averaging of 3D data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_3d_data_averaging( mode, variable )


    CHARACTER (LEN=*), INTENT(IN) ::  mode     !< averaging interface mode
    CHARACTER (LEN=*), INTENT(IN) ::  variable !< variable name


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 3d data averaging', 'start' )

    CALL dynamics_3d_data_averaging( mode, variable )
    CALL tcm_3d_data_averaging( mode, variable )

    IF ( biometeorology      )  CALL bio_3d_data_averaging( mode, variable )
    IF ( bulk_cloud_model    )  CALL bcm_3d_data_averaging( mode, variable )
    IF ( air_chemistry       )  CALL chem_3d_data_averaging( mode, variable )
    CALL doq_3d_data_averaging( mode, variable )  ! ToDo: this seems to be not according to the design
    IF ( gust_module_enabled )  CALL gust_3d_data_averaging( mode, variable )
    IF ( land_surface        )  CALL lsm_3d_data_averaging( mode, variable )
    IF ( ocean_mode          )  CALL ocean_3d_data_averaging( mode, variable )
    IF ( plant_canopy        )  CALL pcm_3d_data_averaging( mode, variable )
    IF ( radiation           )  CALL radiation_3d_data_averaging( mode, variable )
    IF ( salsa               )  CALL salsa_3d_data_averaging( mode, variable )
    IF ( urban_surface       )  CALL usm_3d_data_averaging( mode, variable )

    IF ( user_module_enabled )  CALL user_3d_data_averaging( mode, variable )

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 3d data averaging', 'end' )


 END SUBROUTINE module_interface_3d_data_averaging

!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Define module-specific 2D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_data_output_2d( av, variable, found, grid, mode, local_pf, two_d,     &
                                             nzb_do, nzt_do, fill_value )

    CHARACTER (LEN=*), INTENT(IN)    ::  mode       !< either 'xy', 'xz' or 'yz'
    CHARACTER (LEN=*), INTENT(IN)    ::  variable   !< variable name
    CHARACTER (LEN=*), INTENT(INOUT) ::  grid       !< name of vertical grid

    INTEGER(iwp),      INTENT(IN)    ::  av         !< flag for (non-)average output
    INTEGER(iwp),      INTENT(IN)    ::  nzb_do     !< vertical output index (bottom) (usually 0)
    INTEGER(iwp),      INTENT(IN)    ::  nzt_do     !< vertical output index (top) (usually nz_do3d)

    LOGICAL,           INTENT(INOUT) ::  found      !< flag if output variable is found
    LOGICAL,           INTENT(OUT)   ::  two_d      !< flag for 2D variables

    REAL(wp),          INTENT(IN)    ::  fill_value !< to be removed

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf !< ToDo: can also be kind=sp


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 2d data output', 'start' )

    CALL dynamics_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do,&
                                  fill_value )

    IF ( .NOT. found )  THEN
       CALL tcm_data_output_2d( av, variable, found, grid, mode, local_pf, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  biometeorology )  THEN
       CALL bio_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  bulk_cloud_model )  THEN
       CALL bcm_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  air_chemistry )  THEN
       CALL chem_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do, &
                                 fill_value )
    ENDIF

    IF ( .NOT. found )  THEN
       CALL doq_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do,       &
                           fill_value )
    ENDIF

    IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
       CALL gust_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do, &
                                 fill_value )
    ENDIF

    IF ( .NOT. found  .AND.  land_surface )  THEN
       CALL lsm_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  ocean_mode )  THEN
       CALL ocean_data_output_2d( av, variable, found, grid, mode, local_pf, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  radiation )  THEN
       CALL radiation_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do,    &
                                      nzt_do, fill_value )
    ENDIF

    IF ( .NOT. found  .AND.  salsa )  THEN
       CALL salsa_data_output_2d( av, variable, found, grid, mode, local_pf, two_d, nzb_do, nzt_do )
    ENDIF

    IF ( .NOT. found  .AND.  user_module_enabled )  THEN
       CALL user_data_output_2d( av, variable, found, grid, local_pf, two_d, nzb_do, nzt_do )
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 2d data output', 'end' )


 END SUBROUTINE module_interface_data_output_2d


!--------------------------------------------------------------------------------------------------!
!
! Description:
! ------------
!> Define module-specific 3D output variables
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_data_output_3d( av, variable, found, local_pf, fill_value, resorted,  &
                                             nzb_do, nzt_do )

    CHARACTER (LEN=*), INTENT(IN)    ::  variable   !< variable name

    INTEGER(iwp),      INTENT(IN)    ::  av         !< flag for (non-)average output
    INTEGER(iwp),      INTENT(IN)    ::  nzb_do     !< vertical output index (bottom) (usually 0)
    INTEGER(iwp),      INTENT(IN)    ::  nzt_do     !< vertical output index (top) (usually nz_do3d)

    LOGICAL,           INTENT(INOUT) ::  found      !< flag if output variable is found
    LOGICAL,           INTENT(OUT)   ::  resorted   !< flag if output has been resorted

    REAL(wp),          INTENT(IN)    ::  fill_value !< ToDo: refactor

    REAL(wp), DIMENSION(nxl:nxr,nys:nyn,nzb_do:nzt_do), INTENT(INOUT) ::  local_pf


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 3d data output', 'start' )

    CALL dynamics_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
    resorted = .FALSE.

    IF ( .NOT. found )  THEN
       CALL tcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  biometeorology )  THEN
       CALL bio_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .FALSE.
    ENDIF

    IF ( .NOT. found  .AND.  bulk_cloud_model )  THEN
       CALL bcm_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  air_chemistry )  THEN
       CALL chem_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found )  THEN
       CALL doq_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  gust_module_enabled )  THEN
       CALL gust_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  indoor_model )  THEN
       CALL im_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  ocean_mode )  THEN
       CALL ocean_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  plant_canopy )  THEN
       CALL pcm_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  radiation )  THEN
       CALL radiation_data_output_3d( av, variable, found, local_pf, fill_value, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  salsa )  THEN
       CALL salsa_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( .NOT. found  .AND.  user_module_enabled )  THEN
       CALL user_data_output_3d( av, variable, found, local_pf, nzb_do, nzt_do )
       resorted = .TRUE.
    ENDIF

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific 3d data output', 'end' )


 END SUBROUTINE module_interface_data_output_3d


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Compute module-specific profile and timeseries data
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_statistics( mode, sr, tn, dots_max )


    CHARACTER (LEN=*), INTENT(IN) ::  mode     !< statistical analysis mode

    INTEGER(iwp),      INTENT(IN) ::  dots_max !< maximum number of timeseries
    INTEGER(iwp),      INTENT(IN) ::  sr       !<
    INTEGER(iwp),      INTENT(IN) ::  tn       !<


    IF ( debug_output_timestep )  CALL debug_message( 'module-specific statistics', 'start' )

    CALL dynamics_statistics( mode, sr, tn )
    CALL doq_statistics( mode, sr, tn )

    IF ( gust_module_enabled )  CALL gust_statistics( mode, sr, tn, dots_max )
    IF ( air_chemistry       )  CALL chem_statistics( mode, sr, tn )
    IF ( salsa               )  CALL salsa_statistics( mode, sr, tn )

    IF ( user_module_enabled )  CALL user_statistics( mode, sr, tn )

    IF ( debug_output_timestep )  CALL debug_message( 'module-specific statistics', 'end' )


 END SUBROUTINE module_interface_statistics


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific restart data in Fortran binary format globaly shared by all MPI ranks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_rrd_global_ftn( found )


    LOGICAL, INTENT(INOUT) ::  found    !< flag if variable was found


    IF ( debug_output )  CALL debug_message( 'module-specific read global restart data', 'start' )

    CALL dynamics_rrd_global( found ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL bio_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL bcm_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL flight_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL gust_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL lpm_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL pcm_rrd_global( found )
    IF ( .NOT. found )  CALL ocean_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL salsa_rrd_global( found )
    IF ( .NOT. found )  CALL stg_rrd_global ( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL wtm_rrd_global( found ) ! ToDo: change interface to pass variable
    IF ( .NOT. found )  CALL surface_data_output_rrd_global( found )

    IF ( .NOT. found )  CALL user_rrd_global( found ) ! ToDo: change interface to pass variable

    IF ( debug_output )  CALL debug_message( 'module-specific read global restart data', 'end' )


 END SUBROUTINE module_interface_rrd_global_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific restart data in MPI format globaly shared by all MPI ranks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_rrd_global_mpi


    IF ( debug_output )  CALL debug_message( 'module-specific read global restart data', 'start' )

    CALL dynamics_rrd_global

    IF ( biometeorology )       CALL bio_rrd_global
    IF ( bulk_cloud_model )     CALL bcm_rrd_global
    IF ( virtual_flight )       CALL flight_rrd_global
    IF ( gust_module_enabled )  CALL gust_rrd_global
    IF ( particle_advection )   CALL lpm_rrd_global
    IF ( plant_canopy )         CALL pcm_rrd_global
    IF ( ocean_mode )           CALL ocean_rrd_global
    IF ( salsa )                CALL salsa_rrd_global
    IF ( syn_turb_gen )         CALL stg_rrd_global
    IF ( wind_turbine )         CALL wtm_rrd_global
    IF ( surface_output )       CALL surface_data_output_rrd_global

    IF ( user_module_enabled )  CALL user_rrd_global

    IF ( debug_output )  CALL debug_message( 'module-specific read global restart data', 'end' )


 END SUBROUTINE module_interface_rrd_global_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write module-specific restart data globaly shared by all MPI ranks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_wrd_global


    IF ( debug_output )  CALL debug_message( 'module-specific write global restart data', 'start' )

    CALL dynamics_wrd_global

    IF ( biometeorology )       CALL bio_wrd_global
    IF ( bulk_cloud_model )     CALL bcm_wrd_global
    IF ( virtual_flight )       CALL flight_wrd_global
    IF ( gust_module_enabled )  CALL gust_wrd_global
    IF ( particle_advection )   CALL lpm_wrd_global
    IF ( plant_canopy )         CALL pcm_wrd_global
    IF ( ocean_mode )           CALL ocean_wrd_global
    IF ( salsa )                CALL salsa_wrd_global
    IF ( syn_turb_gen )         CALL stg_wrd_global
    IF ( wind_turbine )         CALL wtm_wrd_global
    IF ( surface_output )       CALL surface_data_output_wrd_global

    IF ( user_module_enabled )  CALL user_wrd_global

    IF ( debug_output )  CALL debug_message( 'module-specific write global restart data', 'end' )


 END SUBROUTINE module_interface_wrd_global


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (Fortran binary format).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_rrd_local_ftn( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxrc, nxr_on_file,                               &
                                            nynf, nync, nyn_on_file,                               &
                                            nysf, nysc, nys_on_file,                               &
                                            tmp_2d, tmp_3d, found )


    INTEGER(iwp), INTENT(IN)  ::  map_index    !<
    INTEGER(iwp), INTENT(IN)  ::  nxlc         !<
    INTEGER(iwp), INTENT(IN)  ::  nxlf         !<
    INTEGER(iwp), INTENT(IN)  ::  nxl_on_file  !<
    INTEGER(iwp), INTENT(IN)  ::  nxrc         !<
    INTEGER(iwp), INTENT(IN)  ::  nxrf         !<
    INTEGER(iwp), INTENT(IN)  ::  nxr_on_file  !<
    INTEGER(iwp), INTENT(IN)  ::  nync         !<
    INTEGER(iwp), INTENT(IN)  ::  nynf         !<
    INTEGER(iwp), INTENT(IN)  ::  nyn_on_file  !<
    INTEGER(iwp), INTENT(IN)  ::  nysc         !<
    INTEGER(iwp), INTENT(IN)  ::  nysf         !<
    INTEGER(iwp), INTENT(IN)  ::  nys_on_file  !<

    LOGICAL,      INTENT(INOUT) ::  found        !< flag if variable was found

    REAL(wp),                                                                                      &
       DIMENSION(nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp),             &
       INTENT(OUT) :: tmp_2d   !<
    REAL(wp),                                                                                      &
       DIMENSION(nzb:nzt+1,nys_on_file-nbgp:nyn_on_file+nbgp,nxl_on_file-nbgp:nxr_on_file+nbgp),   &
       INTENT(OUT) :: tmp_3d   !<


    IF ( debug_output )  THEN
       CALL debug_message( 'module-specific read local restart data (Fortran binary)', 'start' )
    ENDIF

    CALL dynamics_rrd_local( map_index,                                                            &
                             nxlf, nxlc, nxl_on_file,                                              &
                             nxrf, nxrc, nxr_on_file,                                              &
                             nynf, nync, nyn_on_file,                                              &
                             nysf, nysc, nys_on_file,                                              &
                             tmp_2d, tmp_3d, found                                                 &
                           ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL bio_rrd_local( found )

    IF ( .NOT. found )  CALL bcm_rrd_local( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxrc, nxr_on_file,                               &
                                            nynf, nync, nyn_on_file,                               &
                                            nysf, nysc, nys_on_file,                               &
                                            tmp_2d, tmp_3d, found                                  &
                                          ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL chem_rrd_local( map_index,                                            &
                                             nxlf, nxlc, nxl_on_file,                              &
                                             nxrf, nxrc, nxr_on_file,                              &
                                             nynf, nync, nyn_on_file,                              &
                                             nysf, nysc, nys_on_file,                              &
                                             tmp_3d, found                                         &
                                           ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL doq_rrd_local( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxrc, nxr_on_file,                               &
                                            nynf, nync, nyn_on_file,                               &
                                            nysf, nysc, nys_on_file,                               &
                                            tmp_2d, tmp_3d, found                                  &
                                          ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL gust_rrd_local( map_index,                                            &
                                             nxlf, nxlc, nxl_on_file,                              &
                                             nxrf, nxrc, nxr_on_file,                              &
                                             nynf, nync, nyn_on_file,                              &
                                             nysf, nysc, nys_on_file,                              &
                                             tmp_2d, tmp_3d, found                                 &
                                           ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL lpm_rrd_local( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxrc, nxr_on_file,                               &
                                            nynf, nync, nyn_on_file,                               &
                                            nysf, nysc, nys_on_file,                               &
                                            tmp_3d, found                                          &
                                          ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL lsm_rrd_local( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxrc, nxr_on_file,                               &
                                            nynf, nync, nyn_on_file,                               &
                                            nysf, nysc, nys_on_file,                               &
                                            tmp_2d, found                                          &
                                          ) ! ToDo: change interface to pass variable

     IF ( .NOT. found )  CALL pcm_rrd_local( map_index,                                            &
                                             nxlf, nxlc, nxl_on_file,                              &
                                             nxrf, nxrc, nxr_on_file,                              &
                                             nynf, nync, nyn_on_file,                              &
                                             nysf, nysc, nys_on_file,                              &
                                             found                                                 &
                                           )

    IF ( .NOT. found )  CALL ocean_rrd_local( map_index,                                           &
                                              nxlf, nxlc, nxl_on_file,                             &
                                              nxrf, nxrc, nxr_on_file,                             &
                                              nynf, nync, nyn_on_file,                             &
                                              nysf, nysc, nys_on_file,                             &
                                              tmp_3d, found                                        &
                                            ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL radiation_rrd_local( map_index,                                       &
                                                  nxlf, nxlc, nxl_on_file,                         &
                                                  nxrf, nxrc, nxr_on_file,                         &
                                                  nynf, nync, nyn_on_file,                         &
                                                  nysf, nysc, nys_on_file,                         &
                                                  tmp_2d, tmp_3d, found                            &
                                                ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL salsa_rrd_local( map_index,                                           &
                                              nxlf, nxlc, nxl_on_file,                             &
                                              nxrf, nxrc, nxr_on_file,                             &
                                              nynf, nync, nyn_on_file,                             &
                                              nysf, nysc, nys_on_file,                             &
                                              tmp_3d, found                                        &
                                            ) ! ToDo: change interface to pass variable

    IF ( .NOT. found )  CALL usm_rrd_local( map_index,                                             &
                                            nxlf, nxlc, nxl_on_file,                               &
                                            nxrf, nxr_on_file,                                     &
                                            nynf, nyn_on_file,                                     &
                                            nysf, nysc, nys_on_file,                               &
                                            found                                                  &
                                          ) ! ToDo: change interface to pass variable
!
!-- Surface data do not need overlap data, so do not pass these information.
    IF ( .NOT. found )  CALL surface_data_output_rrd_local( found )

    IF ( .NOT. found )  CALL user_rrd_local( map_index,                                            &
                                             nxlf, nxlc, nxl_on_file,                              &
                                             nxrf, nxrc, nxr_on_file,                              &
                                             nynf, nync, nyn_on_file,                              &
                                             nysf, nysc, nys_on_file,                              &
                                             tmp_3d, found                                         &
                                           ) ! ToDo: change interface to pass variable

    IF ( debug_output )  THEN
       CALL debug_message( 'module-specific read local restart data (Fortran binary)', 'end' )
    ENDIF

 END SUBROUTINE module_interface_rrd_local_ftn


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Read module-specific local restart data arrays (MPI-IO).
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_rrd_local_mpi

    IF ( debug_output )  THEN
       CALL debug_message( 'module-specific read local restart data (MPI-IO)', 'start' )
    ENDIF

    CALL dynamics_rrd_local

    IF ( biometeorology )       CALL bio_rrd_local
    IF ( bulk_cloud_model )     CALL bcm_rrd_local
    IF ( air_chemistry )        CALL chem_rrd_local
    CALL doq_rrd_local
    IF ( gust_module_enabled )  CALL gust_rrd_local
    IF ( particle_advection )   CALL lpm_rrd_local
    IF ( land_surface )         CALL lsm_rrd_local
    IF ( plant_canopy )         CALL pcm_rrd_local
    IF ( ocean_mode )           CALL ocean_rrd_local
    IF ( radiation )            CALL radiation_rrd_local
    IF ( salsa )                CALL salsa_rrd_local
    IF ( urban_surface )        CALL usm_rrd_local
    IF ( surface_output )       CALL surface_data_output_rrd_local

    IF ( user_module_enabled )  CALL user_rrd_local

    IF ( debug_output )  THEN
       CALL debug_message( 'module-specific read local restart data (Fortran binary)', 'end' )
    ENDIF

 END SUBROUTINE module_interface_rrd_local_mpi


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Write module-specific restart data specific to local MPI ranks
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_wrd_local


    IF ( debug_output )  CALL debug_message( 'module-specific write local restart data', 'start' )

    CALL dynamics_wrd_local

    IF ( biometeorology )       CALL bio_wrd_local
    IF ( bulk_cloud_model )     CALL bcm_wrd_local
    IF ( air_chemistry )        CALL chem_wrd_local
    CALL doq_wrd_local
    IF ( gust_module_enabled )  CALL gust_wrd_local
    IF ( particle_advection )   CALL lpm_wrd_local
    IF ( land_surface )         CALL lsm_wrd_local
    IF ( plant_canopy )         CALL pcm_wrd_local
    IF ( ocean_mode )           CALL ocean_wrd_local
    IF ( radiation )            CALL radiation_wrd_local
    IF ( salsa )                CALL salsa_wrd_local
    IF ( urban_surface )        CALL usm_wrd_local
    IF ( surface_output )       CALL surface_data_output_wrd_local

    IF ( user_module_enabled )  CALL user_wrd_local

    IF ( debug_output )  CALL debug_message( 'module-specific write local restart data', 'end' )


 END SUBROUTINE module_interface_wrd_local


!--------------------------------------------------------------------------------------------------!
! Description:
! ------------
!> Perform module-specific last actions before the program terminates
!--------------------------------------------------------------------------------------------------!
 SUBROUTINE module_interface_last_actions

    INTEGER ::  return_value  !< returned status value of a called function


    IF ( debug_output )  CALL debug_message( 'module-specific last actions', 'start' )

    return_value = dom_finalize_output()

    CALL dynamics_last_actions
    IF ( particle_advection )   CALL lpm_last_actions
    IF ( user_module_enabled )  CALL user_last_actions

    IF ( debug_output )  CALL debug_message( 'module-specific last actions', 'end' )


 END SUBROUTINE module_interface_last_actions


 END MODULE module_interface
