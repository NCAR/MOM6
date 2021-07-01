!> A tracer package for tracers computed in the MARBL library
!!
!! Currently configured for use with marbl0.36.0
!! https://github.com/marbl-ecosys/MARBL/releases/tag/marbl0.36.0
!! (clone entire repo into pkg/MARBL)
module MARBL_tracers

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,            only : root_PE, broadcast
use MOM_diag_mediator,   only : diag_ctrl
use MOM_diag_vkernels,   only : reintegrate_column
use MOM_error_handler,   only : is_root_PE, MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : file_exists, MOM_read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, MOM_restart_CS, register_restart_field
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,   only : tracer_Z_init, read_Z_edges
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface, thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_diag_mediator,   only : register_diag_field, post_data!, safe_alloc_ptr

use MARBL_interface,              only : MARBL_interface_class
use MARBL_interface_public_types, only : marbl_diagnostics_type, marbl_saved_state_type

use coupler_types_mod,      only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_MARBL_tracers, initialize_MARBL_tracers
public MARBL_tracers_column_physics, MARBL_tracers_surface_state
public MARBL_tracer_stock, MARBL_tracers_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> Temporary type for diagnostic variables coming from MARBL
!! Allocate exactly one of field_[23]d
type :: temp_MARBL_diag
  integer :: id !< index into MOM diagnostic structure
  real, allocatable :: field_2d(:,:) !< memory for 2D field
  real, allocatable :: field_3d(:,:,:) !< memory for 3D field
end type temp_MARBL_diag

!> MOM6 needs to know the index of some MARBL tracers to properly apply river fluxes
type :: tracer_ind_type
  integer :: no3_ind  !< NO3 index
  integer :: po4_ind  !< PO4 index
  integer :: don_ind  !< DON index
  integer :: donr_ind  !< DONr index
  integer :: dop_ind  !< DOP index
  integer :: dopr_ind  !< DOPr index
  integer :: sio3_ind  !< SiO3 index
  integer :: fe_ind  !< Fe index
  integer :: doc_ind  !< DOC index
  integer :: docr_ind  !< DOCr index
  integer :: alk_ind  !< ALK index
  integer :: alk_alt_co2_ind  !< ALK_ALT_CO2 index
  integer :: dic_ind  !< DIC index
  integer :: dic_alt_co2_ind  !< DIC_ALT_CO2 index
end type tracer_ind_type

!> MOM needs to store some information about saved_state; besides providing these
!! fields to MARBL, they are also written to restart files
type :: saved_state_for_MARBL_type
  character(len=200) :: short_name !< name of variable being saved
  character(len=200) :: file_varname !< name of variable in restart file
  character(len=200) :: units !< variable units
  real, pointer :: field_2d(:,:) !< memory for 2D field
  real, pointer :: field_3d(:,:,:) !< memory for 3D field
end type saved_state_for_MARBL_type

!> All calls to MARBL are done via the interface class
type(MARBL_interface_class) :: MARBL_instances

!> The control structure for the MARBL tracer package
type, public :: MARBL_tracers_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  character(len=200) :: IC_file !< The file in which the age-tracer initial values cam be found.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine, in g m-3?

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by aof_set_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .false. !< If true the tracers may be initialized if not found in a restart file

  character(len=200) :: fesedflux_file  !< name of [netCDF] file containing iron sediment flux
  character(len=200) :: feventflux_file  !< name of [netCDF] file containing iron vent flux
  character(len=35) :: marbl_settings_file  !< name of [text] file containing MARBL settings
  real :: atm_co2_const  !< atmospheric CO2 (if specifying a constant value)
  real :: atm_alt_co2_const  !< alternate atmospheric CO2 for _ALT_CO2 tracers (if specifying a constant value)
  real :: ndep_scale_factor  !< scale factor to apply to nitrogen deposition

  type(temp_MARBL_diag), allocatable :: surface_flux_diags(:)  !< collect surface flux diagnostics from all columns
                                                               !! before posting
  type(temp_MARBL_diag), allocatable :: interior_tendency_diags(:)  !< collect tendency diagnostics from all columns
                                                                    !! before posting
  type(saved_state_for_MARBL_type), allocatable :: surface_flux_saved_state(:)  !< surface_flux saved state
  type(saved_state_for_MARBL_type), allocatable :: interior_tendency_saved_state(:)  !< interior_tendency saved state

  ! TODO: If we can post data column by column, all we need are integer arrays for ids
  ! integer, allocatable :: id_surface_flux_diags(:)  !< array of indices for surface_flux diagnostics
  ! integer, allocatable :: id_interior_tendency_diags(:)  !< array of indices for interior_tendency diagnostics

  type(tracer_ind_type) :: tracer_inds  !< Indices to tracers that will have river fluxes added to STF

  !> Need to store global output from both marbl_instance%surface_flux_compute() and
  !! marbl_instance%interior_tendency_compute(). For the former, just need id to register
  !! because we already copy data into CS%STF; latter requires copying data and indices
  !! so currently using temp_MARBL_diag for that.
  integer, allocatable :: id_surface_flux_out(:)  !< register_diag indices for surface_flux output
  type(temp_MARBL_diag), allocatable :: interior_tendency_out(:)  !< collect interior tendencies for diagnostic output
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint(:)  !< vertical integral of interior tendencies
                                                                       !! (full column)
  type(temp_MARBL_diag), allocatable :: interior_tendency_out_zint_100m(:)  !< vertical integral of interior tendencies
                                                                            !! (top 100m)

  ! NOTE: MARBL will return in cgs so we need to convert to mks
  real, allocatable :: STF(:,:,:) !< surface fluxes returned from MARBL to use in tracer_vertdiff [i, j, tracer]

  integer :: u10_sqr_ind  !< index of MARBL forcing field array to copy 10-m wind (squared) into
  integer :: sss_ind  !< index of MARBL forcing field array to copy sea surface salinity into
  integer :: sst_ind  !< index of MARBL forcing field array to copy sea surface temperature into
  integer :: ifrac_ind  !< index of MARBL forcing field array to copy ice fraction into
  integer :: dust_dep_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: fe_dep_ind  !< index of MARBL forcing field array to copy iron flux into
  integer :: nox_flux_ind  !< index of MARBL forcing field array to copy NOx flux into
  integer :: nhy_flux_ind  !< index of MARBL forcing field array to copy NHy flux into
  integer :: atmpress_ind  !< index of MARBL forcing field array to copy atmospheric pressure into
  integer :: xco2_ind  !< index of MARBL forcing field array to copy CO2 flux into
  integer :: xco2_alt_ind  !< index of MARBL forcing field array to copy CO2 flux (alternate CO2) into

  !> Indices for forcing fields required to compute interior tendencies
  integer :: dustflux_ind  !< index of MARBL forcing field array to copy dust flux into
  integer :: PAR_col_frac_ind  !< index of MARBL forcing field array to copy PAR column fraction into
  integer :: surf_shortwave_ind  !< index of MARBL forcing field array to copy surface shortwave into
  integer :: potemp_ind  !< index of MARBL forcing field array to copy potential temperature into
  integer :: salinity_ind  !< index of MARBL forcing field array to copy salinity into
  integer :: pressure_ind  !< index of MARBL forcing field array to copy pressure into
  integer :: fesedflux_ind  !< index of MARBL forcing field array to copy iron sediment flux into
  integer :: o2_scalef_ind  !< index of MARBL forcing field array to copy O2 scale length into
  integer :: remin_scalef_ind  !< index of MARBL forcing field array to copy remin scale length into

  ! TODO: create generic 3D forcing input type to read z coordinate + values
  real    :: fesedflux_scale_factor !< scale factor for iron sediment flux
  integer :: fesedflux_nz  !< number of levels in iron sediment flux file
  real, allocatable, dimension(:,:,:) :: fesedflux_in  !< Field to read iron sediment flux into
  real, allocatable, dimension(:,:,:) :: feventflux_in  !< Field to read iron vent flux into
  real, allocatable, dimension(:) :: &
    fesedflux_z  !< The depths of the cell interfaces in the input data [Z ~> m]
  ! TODO: this thickness does not need to be 3D, but that's a problem for future Mike
  real, allocatable, dimension(:,:,:) :: &
    fesedflux_dz  !< The thickness of the cell layers in the input data [Z ~> m]
end type MARBL_tracers_CS

  ! Module parameters
  real, parameter :: cm_per_m = 100.  !< convert from m -> cm (MARBL is cgs)
  real, parameter :: g_per_kg = 1000. !< convert from kg -> g (MARBL is cgs)
  real, parameter :: m_per_cm = 0.01  !< convert from cm -> m
  real, parameter :: atm_per_Pa = 1./101325.  !< convert from Pa -> atm

contains

!> This subroutine is used to read marbl_in, configure MARBL accordingly, and then
!! call MARBL's initialization routine
subroutine configure_MARBL_tracers(GV, param_file, CS)
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module

#include "version_variable.h"
  character(len=40)  :: mdl = "MARBL_tracers" ! This module's name.
  character(len=256) :: log_message
  character(len=256) :: marbl_in_line(1)
  integer :: m, nz, marbl_settings_in, read_error
  nz = GV%ke
  marbl_settings_in = 615

  ! (1) Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "MARBL_SETTINGS_FILE", CS%marbl_settings_file, &
                 "The name of a file from which to read the run-time "//&
                 "settings for MARBL.", default="marbl_in")
  call get_param(param_file, mdl, "ATM_CO2_CONST", CS%atm_co2_const, &
                 "Value to send to MARBL as xco2", default=284.317)
  call get_param(param_file, mdl, "ATM_ALT_CO2_CONST", CS%atm_alt_co2_const, &
                 "Value to send to MARBL as xco2_alt_co2", default=284.317)
  call get_param(param_file, mdl, "NDEP_SCALE_FACTOR", CS%ndep_scale_factor, &
                 "Scale factor applied to nitrogen deposition terms", default=1e5)

  ! (2) Read marbl settings file and call put_setting()

  ! (2a) only master task opens file
  if (is_root_PE()) then
     ! read the marbl_in into buffer
     open(unit=marbl_settings_in, file=CS%marbl_settings_file, iostat=read_error)
     if (read_error .ne. 0) then
        write(log_message, '(A, I0, 2A)') "IO ERROR ", read_error, &
              "opening namelist file : ", trim(CS%marbl_settings_file)
        call MOM_error(FATAL, log_message)
     end if
  end if

  ! (2b) master task reads file and broadcasts line-by-line
  marbl_in_line = ''
  do
    ! i. Read next line on master, iostat value out
    !    (Exit loop if read is not successful; either read error or end of file)
    if (is_root_PE()) read(marbl_settings_in, "(A)", iostat=read_error) marbl_in_line(1)
    call broadcast(read_error, root_PE())
    if (read_error .ne. 0) exit

    ! ii. Broadcast line just read in on root PE to all tasks
    call broadcast(marbl_in_line, 256, root_PE())

    ! iii. All tasks call put_setting (TODO: openMP blocks?)
    call MARBL_instances%put_setting(marbl_in_line(1))
  end do
  ! iv. (TEMPORARY) don't set tracer restoring
  call MARBL_instances%put_setting("tracer_restore_vars(1) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(2) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(3) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(4) = ''")
  call MARBL_instances%put_setting("tracer_restore_vars(5) = ''")

  ! (2c) we should always reach the EOF to capture the entire file...
  if (.not. is_iostat_end(read_error)) then
     write(log_message, '(3A, I0)') "IO ERROR reading ", trim(CS%marbl_settings_file), ": ", read_error
     call MOM_error(FATAL, log_message)
  else
     if (is_root_PE()) then
       write(log_message, '(3A)') "Read '", trim(CS%marbl_settings_file), "' until EOF."
       call MOM_error(NOTE, log_message)
     end if
  end if
  if (is_root_PE()) close(marbl_settings_in)

  ! (3) call marbl%init()
  ! TODO: the units in gcm_delta_z, gcm_zw, and gcm_zt are wrong, but we want to strip these values
  !       out of init anyway because MOM updates them every time step / every column
  call MARBL_instances%init(&
                            gcm_num_levels = nz, &
                            gcm_num_PAR_subcols = 1, &
                            gcm_num_elements_surface_flux = 1, & ! FIXME: change to number of grid cells on MPI task
                            gcm_delta_z = GV%sInterface(2:nz+1) - GV%sInterface(1:nz), &
                            gcm_zw = GV%sInterface(2:nz+1), &
                            gcm_zt = GV%sLayer, &
                            lgcm_has_global_ops = .true. &
                           )
  if (MARBL_instances%StatusLog%labort_marbl) &
    call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%init", "configure_MARBL_tracers")
  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()

  ! (4) Initialize forcing fields
  !     i. store all surface forcing indices
  CS%u10_sqr_ind = -1
  CS%sss_ind = -1
  CS%sst_ind = -1
  CS%ifrac_ind = -1
  CS%dust_dep_ind = -1
  CS%fe_dep_ind = -1
  CS%nox_flux_ind = -1
  CS%nhy_flux_ind = -1
  CS%atmpress_ind = -1
  CS%xco2_ind = -1
  CS%xco2_alt_ind = -1
  do m=1,size(MARBL_instances%surface_flux_forcings)
    select case (trim(MARBL_instances%surface_flux_forcings(m)%metadata%varname))
      case('u10_sqr')
        CS%u10_sqr_ind = m
      case('sss')
        CS%sss_ind = m
      case('sst')
        CS%sst_ind = m
      case('Ice Fraction')
        CS%ifrac_ind = m
      case('Dust Flux')
        CS%dust_dep_ind = m
      case('Iron Flux')
        CS%fe_dep_ind = m
      case('NOx Flux')
        CS%nox_flux_ind = m
      case('NHy Flux')
        CS%nhy_flux_ind = m
      case('Atmospheric Pressure')
        CS%atmpress_ind = m
      case('xco2')
        CS%xco2_ind = m
      case('xco2_alt_co2')
        CS%xco2_alt_ind = m
      case DEFAULT
        write(log_message, "(A,1X,A)") trim(MARBL_instances%surface_flux_forcings(m)%metadata%varname), &
                                   'is not a valid surface flux forcing field name.'
        call MOM_error(FATAL, log_message)
    end select
  end do

  !     ii. store all surface forcing indices
  CS%dustflux_ind = -1
  CS%PAR_col_frac_ind = -1
  CS%surf_shortwave_ind = -1
  CS%potemp_ind = -1
  CS%salinity_ind = -1
  CS%pressure_ind = -1
  CS%fesedflux_ind = -1
  CS%o2_scalef_ind = -1
  CS%remin_scalef_ind = -1
  do m=1,size(MARBL_instances%interior_tendency_forcings)
    ! i. Check to see if this is a tracer restoring field or timescale
    ! ii. If not, should be one of the following:
    select case (trim(MARBL_instances%interior_tendency_forcings(m)%metadata%varname))
      case('Dust Flux')
        CS%dustflux_ind = m
      case('PAR Column Fraction')
        CS%PAR_col_frac_ind = m
      case('Surface Shortwave')
        CS%surf_shortwave_ind = m
      case('Potential Temperature')
        CS%potemp_ind = m
      case('Salinity')
        CS%salinity_ind = m
      case('Pressure')
        CS%pressure_ind = m
      case('Iron Sediment Flux')
        CS%fesedflux_ind = m
      case('O2 Consumption Scale Factor')
        CS%o2_scalef_ind = m
      case('Particulate Remin Scale Factor')
        CS%remin_scalef_ind = m
      case DEFAULT
        write(log_message, "(A,1X,A)") trim(MARBL_instances%interior_tendency_forcings(m)%metadata%varname), &
                                   'is not a valid interior tendency forcing field name.'
        call MOM_error(FATAL, log_message)
    end select
  end do
end subroutine configure_MARBL_tracers

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
function register_MARBL_tracers(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(MARBL_tracers_CS),     pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure.

! Local variables
! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MARBL_tracers" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=128) :: desc_name ! The variable's descriptor.
  character(len=48)  :: units ! The variable's units.
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_MARBL_tracers
  integer :: isd, ied, jsd, jed, nz, m
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_MARBL_tracers called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  call configure_MARBL_tracers(GV, param_file, CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  ! ** Input directory
  ! TODO: just use DIN_LOC_ROOT
  call get_param(param_file, mdl, "CESM_INPUTDIR", inputdir, default="/glade/work/mlevy/cesm_inputdata")
  ! ** Tracer initial conditions
  call get_param(param_file, mdl, "MARBL_TRACERS_IC_FILE", CS%IC_file, &
                 "The file in which the MARBL tracers initial values can be found.", &
                 default="ecosys_jan_IC_omip_MOM_tx0.66v1_c190925.nc")
  if (scan(CS%IC_file,'/') == 0) then
    ! Add the directory if CS%IC_file is not already a complete path.
    CS%IC_file = trim(slasher(inputdir))//trim(CS%IC_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_IC_FILE", CS%IC_file)
  endif
  ! ** FESEDFLUX
  call get_param(param_file, mdl, "MARBL_FESEDFLUX_FILE", CS%fesedflux_file, &
                 "The file in which the iron sediment flux forcing field can be found.", &
                 default="fesedflux_total_reduce_oxic_tx0.66v1.c201204.nc")
  if (scan(CS%fesedflux_file,'/') == 0) then
    ! Add the directory if CS%fesedflux_file is not already a complete path.
    CS%fesedflux_file = trim(slasher(inputdir))//trim(CS%fesedflux_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_FESEDFLUX_FILE", CS%fesedflux_file)
  endif
  ! ** FEVENTFLUX
  call get_param(param_file, mdl, "MARBL_FEVENTFLUX_FILE", CS%feventflux_file, &
                 "The file in which the iron vent flux forcing field can be found.", &
                 default="feventflux_5gmol_tx0.66v1.c201204.nc")
  if (scan(CS%feventflux_file,'/') == 0) then
    ! Add the directory if CS%feventflux_file is not already a complete path.
    CS%feventflux_file = trim(slasher(inputdir))//trim(CS%feventflux_file)
    call log_param(param_file, mdl, "INPUTDIR/MARBL_TRACERS_FEVENTFLUX_FILE", CS%feventflux_file)
  endif
  ! ** Scale factor for FESEDFLUX
  call get_param(param_file, mdl, "MARBL_FESEDFLUX_SCALE_FACTOR", CS%fesedflux_scale_factor, &
                 "Conversion factor between FESEDFLUX file and MARBL units (umol / m^2 / d -> nmol / cm^2 / s)", &
                 default=1000. / 10000. / 86400.)

  CS%ntr = size(MARBL_instances%tracer_metadata)
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))
  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr)) ; CS%tr(:,:,:,:) = 0.0

  do m = 1, CS%ntr
    write(var_name(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(desc_name(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%long_name)
    write(units(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%units)
    CS%tr_desc(m) = var_desc(trim(var_name), trim(units), trim(desc_name), caller=mdl)

    ! This is needed to force the compiler not to do a copy in the registration
    ! calls.  Curses on the designers and implementers of Fortran90.
    tr_ptr => CS%tr(:,:,:,m)
    call query_vardesc(CS%tr_desc(m), name=var_name, &
                       caller="register_MARBL_tracers")
    ! Register the tracer for horizontal advection, diffusion, and restarts.
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, units = units, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = aof_set_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_MARBL_tracers")
  enddo

  ! Set up memory for saved state
  call setup_saved_state(MARBL_instances%surface_flux_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, &
                         CS%surface_flux_saved_state)
  call setup_saved_state(MARBL_instances%interior_tendency_saved_state, HI, GV, restart_CS, CS%tracers_may_reinit, &
                         CS%interior_tendency_saved_state)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_MARBL_tracers = .true.
end function register_MARBL_tracers

!> This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_MARBL_tracers(restart, day, G, GV, US, h, diag, OBC, CS, sponge_CSp)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(MARBL_tracers_CS),                pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_MARBL_tracers.
  type(sponge_CS),                    pointer    :: sponge_CSp    !< A pointer to the control structure
                                                                  !! for the sponges, if they are in use.

! Local variables
  character(len=200) :: log_message
  character(len=48) :: name       ! A variable's name in a NetCDF file.
  character(len=100) :: longname   ! The long name of that variable.
  character(len=48) :: units      ! The units of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                  ! years m3 s-1 or years kg s-1.
  character(len=48) :: tracer_name
  logical :: OK
  logical :: fesedflux_has_edges, fesedflux_use_missing
  real    :: fesedflux_missing
  integer :: i, j, k, kbot, m, diag_size

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Allocate memory for surface tracer fluxes
  allocate(CS%STF(SZI_(G), SZJ_(G), CS%ntr))
  CS%STF(:,:,:) = 0.

  ! Register diagnostics returned from MARBL (surface flux first, then interior tendency)
  call register_MARBL_diags(MARBL_instances%surface_flux_diags, diag, day, G, CS%surface_flux_diags)
  call register_MARBL_diags(MARBL_instances%interior_tendency_diags, diag, day, G, CS%interior_tendency_diags)

  ! Register per-tracer diagnostics computed from MARBL surface flux / interior tendency values
  allocate(CS%id_surface_flux_out(CS%ntr))
  allocate(CS%interior_tendency_out(CS%ntr))
  allocate(CS%interior_tendency_out_zint(CS%ntr))
  allocate(CS%interior_tendency_out_zint_100m(CS%ntr))
  do m=1,CS%ntr
    write(name, "(2A)") "STF_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Surface Flux"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%id_surface_flux_out(m) = register_diag_field("ocean_model", &
                                                    trim(name), &
                                                    diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
                                                    day, &
                                                    trim(longname), &
                                                    trim(units))

    write(name, "(2A)") "J_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Source Sink Term"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), "/s"
    CS%interior_tendency_out(m)%id = register_diag_field("ocean_model", &
                                                         trim(name), &
                                                         diag%axesTL, & ! T=> tracer grid? L => layer center
                                                         day, &
                                                         trim(longname), &
                                                         trim(units))
    if (CS%interior_tendency_out(m)%id > 0) then
      allocate(CS%interior_tendency_out(m)%field_3d(SZI_(G),SZJ_(G), SZK_(G)))
      CS%interior_tendency_out(m)%field_3d(:,:,:) = 0.
    end if

    write(name, "(2A)") "Jint_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), " Source Sink Term Vertical Integral"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%interior_tendency_out_zint(m)%id = register_diag_field("ocean_model", &
                                                              trim(name), &
                                                              diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                              day, &
                                                              trim(longname), &
                                                              trim(units))
    if (CS%interior_tendency_out_zint(m)%id > 0) then
      allocate(CS%interior_tendency_out_zint(m)%field_2d(SZI_(G),SZJ_(G)))
      CS%interior_tendency_out_zint(m)%field_2d(:,:) = 0.
    end if

    write(name, "(2A)") "Jint_100m_", trim(MARBL_instances%tracer_metadata(m)%short_name)
    write(longname, "(2A)") trim(MARBL_instances%tracer_metadata(m)%long_name), &
                            " Source Sink Term Vertical Integral, 0-100m"
    write(units, "(2A)") trim(MARBL_instances%tracer_metadata(m)%units), " m/s"
    CS%interior_tendency_out_zint_100m(m)%id = register_diag_field("ocean_model", &
                                                                   trim(name), &
                                                                   diag%axesT1, &
                                                                   day, &
                                                                   trim(longname), &
                                                                   trim(units))
    if (CS%interior_tendency_out_zint_100m(m)%id > 0) then
      allocate(CS%interior_tendency_out_zint_100m(m)%field_2d(SZI_(G),SZJ_(G)))
      CS%interior_tendency_out_zint_100m(m)%field_2d(:,:) = 0.
    end if

  end do

  ! Initialize tracers from file (unless they were initialized by restart file)
  ! Also save indices of tracers that have river fluxes
  CS%tracer_inds%no3_ind = 0
  CS%tracer_inds%po4_ind = 0
  CS%tracer_inds%don_ind = 0
  CS%tracer_inds%donr_ind = 0
  CS%tracer_inds%dop_ind = 0
  CS%tracer_inds%dopr_ind = 0
  CS%tracer_inds%sio3_ind = 0
  CS%tracer_inds%fe_ind = 0
  CS%tracer_inds%doc_ind = 0
  CS%tracer_inds%docr_ind = 0
  CS%tracer_inds%alk_ind = 0
  CS%tracer_inds%alk_alt_co2_ind = 0
  CS%tracer_inds%dic_ind = 0
  CS%tracer_inds%dic_alt_co2_ind = 0
  if (.not.restart) then
    do m= 1, CS%ntr
      write(name(:),'(A)') trim(MARBL_instances%tracer_metadata(m)%short_name)
      if (trim(name) == "NO3") then
        CS%tracer_inds%no3_ind = m
      elseif (trim(name) == "PO4") then
         CS%tracer_inds%po4_ind = m
      elseif (trim(name) == "DON") then
         CS%tracer_inds%don_ind = m
      elseif (trim(name) == "DONr") then
         CS%tracer_inds%donr_ind = m
      elseif (trim(name) == "DOP") then
         CS%tracer_inds%dop_ind = m
      elseif (trim(name) == "DOPr") then
         CS%tracer_inds%dopr_ind = m
      elseif (trim(name) == "SiO3") then
         CS%tracer_inds%sio3_ind = m
        elseif (trim(name) == "Fe") then
          CS%tracer_inds%fe_ind = m
       elseif (trim(name) == "DOC") then
          CS%tracer_inds%doc_ind = m
       elseif (trim(name) == "DOCr") then
          CS%tracer_inds%docr_ind = m
       elseif (trim(name) == "ALK") then
          CS%tracer_inds%alk_ind = m
       elseif (trim(name) == "ALK_ALT_CO2") then
          CS%tracer_inds%alk_alt_co2_ind = m
       elseif (trim(name) == "DIC") then
          CS%tracer_inds%dic_ind = m
       elseif (trim(name) == "DIC_ALT_CO2") then
          CS%tracer_inds%dic_alt_co2_ind = m
       end if
      OK = tracer_Z_init(CS%tr(:,:,:,m), h, CS%IC_file, name, G, GV, US, -1e34)
      if (.not.OK) call MOM_error(FATAL,"initialize_MARBL_tracers: "//&
                                  "Unable to read "//trim(name)//" from "//&
                                  trim(CS%IC_file)//".")
    end do
  end if

  ! Log indices for each tracer to ensure we set them all correctly
  write(log_message, "(A,I0)") "NO3 index: ", CS%tracer_inds%no3_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "PO4 index: ", CS%tracer_inds%po4_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DON index: ", CS%tracer_inds%don_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DONr index: ", CS%tracer_inds%donr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOP index: ", CS%tracer_inds%dop_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOPr index: ", CS%tracer_inds%dopr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "SiO3 index: ", CS%tracer_inds%sio3_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "Fe index: ", CS%tracer_inds%fe_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOC index: ", CS%tracer_inds%doc_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DOCr index: ", CS%tracer_inds%docr_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "ALK index: ", CS%tracer_inds%alk_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "ALK_ALT_CO2 index: ", CS%tracer_inds%alk_alt_co2_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DIC index: ", CS%tracer_inds%dic_ind
  call MOM_error(NOTE, log_message)
  write(log_message, "(A,I0)") "DIC_ALT_CO2 index: ", CS%tracer_inds%dic_alt_co2_ind
  call MOM_error(NOTE, log_message)

  ! Read initial fesedflux and feventflux fields
  ! (1) get vertical dimension
  !     -- comes from fesedflux_file, assume same dimension in feventflux
  !        (maybe these fields should be combined?)
  !     -- note: read_Z_edges treats depth as positive UP => 0 at surface, negative at depth
  fesedflux_use_missing = .false.
  call read_Z_edges(CS%fesedflux_file, "FESEDFLUXIN", CS%fesedflux_z, CS%fesedflux_nz, &
                    fesedflux_has_edges, fesedflux_use_missing, fesedflux_missing, &
                    scale=US%m_to_Z)

  ! (2) Allocate memory for fesedflux and feventflux
  allocate(CS%fesedflux_in(SZI_(G), SZJ_(G), CS%fesedflux_nz))
  allocate(CS%feventflux_in(SZI_(G), SZJ_(G), CS%fesedflux_nz))
  allocate(CS%fesedflux_dz(SZI_(G), SZJ_(G), CS%fesedflux_nz-1))

  ! (3) Read data
  !     TODO: Add US term to scale
  call MOM_read_data(CS%fesedflux_file, "FESEDFLUXIN", CS%fesedflux_in(:,:,:), G%Domain, &
                     scale=CS%fesedflux_scale_factor)
  call MOM_read_data(CS%feventflux_file, "FESEDFLUXIN", CS%feventflux_in(:,:,:), G%Domain, &
                     scale=CS%fesedflux_scale_factor)

  ! (4) Relocate values that are below ocean bottom to layer that intersects bathymetry
  !     Remember, fesedflux_z = 0 at surface and is < 0 below surface
  do k=CS%fesedflux_nz-1, 1, -1
    kbot = k + 1 ! level k is between z(k) and z(k+1)
    do j=G%jsc, G%jec
      do i=G%isc, G%iec
        ! Also figure out layer thickness while we're here
        CS%fesedflux_dz(i,j,k) = CS%fesedflux_z(k) - CS%fesedflux_z(kbot)
        ! If top interface is below ocean bottom, bring flux from bottom layer up
        ! and set thickness of level below to 0
        if (-CS%fesedflux_z(k) > G%bathyT(i,j)) then
          CS%fesedflux_in(i,j,k) = CS%fesedflux_in(i,j,k) + CS%fesedflux_in(i,j,kbot)
          CS%fesedflux_in(i,j,kbot) = 0.
          CS%feventflux_in(i,j,k) = CS%feventflux_in(i,j,k) + CS%feventflux_in(i,j,kbot)
          CS%feventflux_in(i,j,kbot) = 0.
          CS%fesedflux_dz(i,j,k) = 0.
        else if (-CS%fesedflux_z(kbot) > G%bathyT(i,j)) then
          ! Otherwise, if lower interface is below bathymetry move interface to ocean bottom
          CS%fesedflux_dz(i,j,k) = G%bathyT(i,j) + CS%fesedflux_z(k)
        end if
      end do
    end do
  end do

end subroutine initialize_MARBL_tracers

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
subroutine register_MARBL_diags(MARBL_diags, diag, day, G, id_diags)

  type(marbl_diagnostics_type), intent(in)    :: MARBL_diags !< MARBL diagnostics from MARBL_instances
  type(time_type), target,      intent(in)    :: day  !< Time of the start of the run.
  type(diag_ctrl), target,      intent(in)    :: diag !< Structure used to regulate diagnostic output.
  !integer, allocatable,         intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index number
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(temp_marbl_diag), allocatable, intent(inout) :: id_diags(:) !< allocatable array storing diagnostic index
                                                                   !! number and buffer space for collecting diags
                                                                   !! from all columns

  integer :: m, diag_size

  diag_size = size(MARBL_diags%diags)
  allocate(id_diags(diag_size))
  do m = 1, diag_size
    id_diags(m)%id = -1
    if (trim(MARBL_diags%diags(m)%vertical_grid) .eq. "none") then ! 2D field
      id_diags(m)%id = register_diag_field("ocean_model", &
        trim(MARBL_diags%diags(m)%short_name), &
        diag%axesT1, & ! T => tracer grid? 1 => no vertical grid
        day, &
        trim(MARBL_diags%diags(m)%long_name), &
        trim(MARBL_diags%diags(m)%units))
      allocate(id_diags(m)%field_2d(SZI_(G),SZJ_(G)))
      id_diags(m)%field_2d(:,:) = 0.
    else ! 3D field
      id_diags(m)%id = register_diag_field("ocean_model", &
        trim(MARBL_diags%diags(m)%short_name), &
        diag%axesTL, & ! T=> tracer grid? L => layer center
        day, &
        trim(MARBL_diags%diags(m)%long_name), &
        trim(MARBL_diags%diags(m)%units))
      allocate(id_diags(m)%field_3d(SZI_(G),SZJ_(G), SZK_(G)))
      id_diags(m)%field_3d(:,:,:) = 0.
    end if
  end do

end subroutine register_MARBL_diags

!> This subroutine allocates memory for saved state fields and registers them in the restart files
subroutine setup_saved_state(MARBL_saved_state, HI, GV, restart_CS, tracers_may_reinit, local_saved_state)

  type(marbl_saved_state_type),                  intent(in)    :: MARBL_saved_state !< MARBL saved state from
                                                                                    !! MARBL_instances
  type(hor_index_type),                          intent(in)    :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),                       intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MOM_restart_CS), pointer,                 intent(in)    :: restart_CS !< control structure to add saved state
                                                                             !! to restarts
  logical,                                       intent(in)    :: tracers_may_reinit  !< used to determine mandatory
                                                                                      !! flag in restart
  type(saved_state_for_MARBL_type), allocatable, intent(inout) :: local_saved_state(:) !< allocatable array for local
                                                                                       !! saved state

  integer :: num_fields, m
  character(len=200) :: log_message, varname

  num_fields = MARBL_saved_state%saved_state_cnt
  allocate(local_saved_state(num_fields))

  do m=1,num_fields
    write(varname, "(2A)") "MARBL_", trim(MARBL_saved_state%state(m)%short_name)
    select case (MARBL_saved_state%state(m)%rank)
      case (2)
        allocate(local_saved_state(m)%field_2d(SZI_(HI),SZJ_(HI)))
        local_saved_state(m)%field_2d(:,:) = 0.
        call register_restart_field(local_saved_state(m)%field_2d, varname, .not.tracers_may_reinit, restart_CS)
      case (3)
        if (trim(MARBL_saved_state%state(m)%vertical_grid).eq."layer_avg") then
          allocate(local_saved_state(m)%field_3d(SZI_(HI),SZJ_(HI), SZK_(GV)))
          local_saved_state(m)%field_3d(:,:,:) = 0.
          call register_restart_field(local_saved_state(m)%field_3d, varname, .not.tracers_may_reinit, restart_CS)
        else
          write(log_message, "(3A, I0, A)") "'", trim(MARBL_saved_state%state(m)%vertical_grid), &
                "' is an invalid vertical grid for saved state (ind = ", m, ")"
          call MOM_error(FATAL, log_message)
        end if
      case DEFAULT
        write(log_message, "(I0, A, I0, A)") MARBL_saved_state%state(m)%rank, &
              " is an invalid rank for saved state (ind = ", m, ")"
        call MOM_error(FATAL, log_message)
    end select
    local_saved_state(m)%short_name = trim(MARBL_saved_state%state(m)%short_name)
    write(local_saved_state(m)%file_varname, "(2A)") "MARBL_", trim(local_saved_state(m)%short_name)
    local_saved_state(m)%units = trim(MARBL_saved_state%state(m)%units)
  end do

end subroutine setup_saved_state

!> This subroutine applies diapycnal diffusion and any other column
!! tracer physics or chemistry to the tracers from this file.
subroutine MARBL_tracers_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(MARBL_tracers_CS),     pointer :: CS   !< The control structure returned by a previous
                                              !! call to register_MARBL_tracers.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [m]

! Local variables
  character(len=256) :: log_message
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: sfc_val  ! The surface value for the tracers.
  real :: Isecs_per_year  ! The number of seconds in a year.
  real :: year            ! The time in years.
  integer :: secs, days   ! Integer components of the time type.
  real, dimension(0:GV%ke) :: zi  ! z-coordinate interface depth
  real, dimension(GV%ke) :: zc, dz  ! z-coordinate layer center depth and cell thickness
  integer :: i, j, k, is, ie, js, je, nz, m, kmt ! TODO: kmt might be temporary?
  real :: kg_m2_s_conversion, ndep_conversion
  logical :: set_kmt ! Temporary variable; true if kmt has been set based on dz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  kg_m2_s_conversion = US%RZ_T_to_kg_m2s
  ndep_conversion = (US%m_to_L)**2 * US%s_to_T

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  ! (1) Compute surface fluxes
  ! FIXME: MARBL can handle computing surface fluxes for all columns simultaneously
  !        I was just thinking going column-by-column at first might be easier
  do i=is,ie
    do j=js,je
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Load proper column data
      !     * surface flux forcings
      !       These fields are getting the correct data
      !       TODO: if top layer is vanishly thin, do we actually want (e.g.) top 5m average temp / salinity?
      !             How does MOM pass SST and SSS to GFDL coupler? (look in core.F90?)
      if (CS%sss_ind > 0) MARBL_instances%surface_flux_forcings(CS%sss_ind)%field_0d(1) = tv%S(i,j,1)
      if (CS%sst_ind > 0) MARBL_instances%surface_flux_forcings(CS%sst_ind)%field_0d(1) = tv%T(i,j,1)
      if (CS%ifrac_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%ifrac_ind)%field_0d(1) = fluxes%MARBL_forcing%ice_fraction(i,j)
      ! MARBL wants u10_sqr in (cm/s)^2
      if (CS%u10_sqr_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%u10_sqr_ind)%field_0d(1) = fluxes%MARBL_forcing%u10_sqr(i,j) * &
                                                                            (US%L_t_to_m_s * cm_per_m)**2
      ! mct_driver/ocn_cap_methods:93 -- ice_ocean_boundary%p(i,j) comes from coupler
      ! We may need a new ice_ocean_boundary%p_atm because %p includes ice in GFDL driver
      if (CS%atmpress_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%atmpress_ind)%field_0d(1) = fluxes%p_surf_full(i,j) * &
                                                                             ((US%R_to_kg_m3 * US%L_T_to_m_s**2) * &
                                                                              atm_per_Pa)

      !       These are okay, but need option to come in from coupler
      if (CS%xco2_ind > 0) MARBL_instances%surface_flux_forcings(CS%xco2_ind)%field_0d(1) = CS%atm_co2_const
      if (CS%xco2_alt_ind > 0) MARBL_instances%surface_flux_forcings(CS%xco2_alt_ind)%field_0d(1) = CS%atm_alt_co2_const

      !       These are okay, but need option to read in from file
      if (CS%dust_dep_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%dust_dep_ind)%field_0d(1) = fluxes%MARBL_forcing%dust_flux(i,j) * &
                                                                             (kg_m2_s_conversion * g_per_kg * &
                                                                              m_per_cm**2)
      if (CS%fe_dep_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%fe_dep_ind)%field_0d(1) = fluxes%MARBL_forcing%iron_flux(i,j) * &
                                                                           (kg_m2_s_conversion * g_per_kg * m_per_cm**2)

      !       These are read from /glade/work/mlevy/cesm_inputdata/ndep_ocn_1850_w_nhx_emis_MOM_tx0.66v1_c210222.nc
      if (CS%nox_flux_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%nox_flux_ind)%field_0d(1) = fluxes%MARBL_forcing%noy_dep(i,j) * &
                                                                             (ndep_conversion * CS%ndep_scale_factor)
      if (CS%nhy_flux_ind > 0) &
        MARBL_instances%surface_flux_forcings(CS%nhy_flux_ind)%field_0d(1) = fluxes%MARBL_forcing%nhx_dep(i,j) * &
                                                                             (ndep_conversion * CS%ndep_scale_factor)

      !     * tracers at surface
      !       TODO: average over some shallow depth (e.g. 5m)
      do m=1,CS%ntr
        MARBL_instances%tracers_at_surface(1,m) = CS%tr(i,j,1,m)
      end do

      !     * surface flux saved state
      do m=1,size(MARBL_instances%surface_flux_saved_state%state)
      !       (currently only 2D fields are saved from surface_flux_compute())
        MARBL_instances%surface_flux_saved_state%state(m)%field_2d(1) = CS%surface_flux_saved_state(m)%field_2d(i,j)
      end do

      ! iii. Compute surface fluxes in MARBL
      call MARBL_instances%surface_flux_compute()
      if (MARBL_instances%StatusLog%labort_marbl) then
        call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%surface_flux_compute()", &
                                                       "MARBL_tracers_column_physics")
        call print_marbl_log(MARBL_instances%StatusLog)
      end if

      ! iv. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(MARBL_instances%surface_flux_saved_state%state)
        CS%surface_flux_saved_state(m)%field_2d(i,j) = MARBL_instances%surface_flux_saved_state%state(m)%field_2d(1)
      end do

      !     * diagnostics
      do m=1,size(MARBL_instances%surface_flux_diags%diags)
        ! All diags are 2D coming from surface
        CS%surface_flux_diags(m)%field_2d(i,j) = real(MARBL_instances%surface_flux_diags%diags(m)%field_2d(1))
      end do

      !     * Surface tracer flux
      CS%STF(i,j,:) = MARBL_instances%surface_fluxes(1,:) * m_per_cm
    end do
  end do
  ! Add River Fluxes to STF
  if (CS%tracer_inds%no3_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%no3_ind) = CS%STF(:,:,CS%tracer_inds%no3_ind) + fluxes%MARBL_forcing%no3_riv_flux(:,:)
  if (CS%tracer_inds%po4_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%po4_ind) = CS%STF(:,:,CS%tracer_inds%po4_ind) + fluxes%MARBL_forcing%po4_riv_flux(:,:)
  if (CS%tracer_inds%don_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%don_ind) = CS%STF(:,:,CS%tracer_inds%don_ind) + fluxes%MARBL_forcing%don_riv_flux(:,:)
  if (CS%tracer_inds%donr_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%donr_ind) = CS%STF(:,:,CS%tracer_inds%donr_ind) + fluxes%MARBL_forcing%donr_riv_flux(:,:)
  if (CS%tracer_inds%dop_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%dop_ind) = CS%STF(:,:,CS%tracer_inds%dop_ind) + fluxes%MARBL_forcing%dop_riv_flux(:,:)
  if (CS%tracer_inds%dopr_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%dopr_ind) = CS%STF(:,:,CS%tracer_inds%dopr_ind) + fluxes%MARBL_forcing%dopr_riv_flux(:,:)
  if (CS%tracer_inds%sio3_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%sio3_ind) = CS%STF(:,:,CS%tracer_inds%sio3_ind) + fluxes%MARBL_forcing%sio3_riv_flux(:,:)
  if (CS%tracer_inds%fe_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%fe_ind) = CS%STF(:,:,CS%tracer_inds%fe_ind) + fluxes%MARBL_forcing%fe_riv_flux(:,:)
  if (CS%tracer_inds%doc_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%doc_ind) = CS%STF(:,:,CS%tracer_inds%doc_ind) + fluxes%MARBL_forcing%doc_riv_flux(:,:)
  if (CS%tracer_inds%docr_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%docr_ind) = CS%STF(:,:,CS%tracer_inds%docr_ind) + fluxes%MARBL_forcing%docr_riv_flux(:,:)
  if (CS%tracer_inds%alk_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%alk_ind) = CS%STF(:,:,CS%tracer_inds%alk_ind) + fluxes%MARBL_forcing%alk_riv_flux(:,:)
  if (CS%tracer_inds%alk_alt_co2_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%alk_alt_co2_ind) = CS%STF(:,:,CS%tracer_inds%alk_alt_co2_ind) + &
                                                 fluxes%MARBL_forcing%alk_alt_co2_riv_flux(:,:)
  if (CS%tracer_inds%dic_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%dic_ind) = CS%STF(:,:,CS%tracer_inds%dic_ind) + fluxes%MARBL_forcing%dic_riv_flux(:,:)
  if (CS%tracer_inds%dic_alt_co2_ind > 0) &
    CS%STF(:,:,CS%tracer_inds%dic_alt_co2_ind) = CS%STF(:,:,CS%tracer_inds%dic_alt_co2_ind) + &
                                                 fluxes%MARBL_forcing%dic_alt_co2_riv_flux(:,:)

  ! (2) Post surface fluxes and their diagnostics (currently all 2D)
  do m=1,CS%ntr
    if (CS%id_surface_flux_out(m) > 0) &
      call post_data(CS%id_surface_flux_out(m), CS%STF(:,:,m), CS%diag)
  end do
  do m=1,size(CS%surface_flux_diags)
    if (CS%surface_flux_diags(m)%id > 0) &
      call post_data(CS%surface_flux_diags(m)%id, CS%surface_flux_diags(m)%field_2d(:,:), CS%diag)
  end do

  ! (3) Apply surface fluxes via vertical diffusion
  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ;do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m) , dt, fluxes, h_work, &
          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV, sfc_flux=CS%STF(:,:,m))
    enddo
  else
    do m=1,CS%ntr
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV, sfc_flux=CS%STF(:,:,m))
    enddo
  endif

  ! (4) Compute interior tendencies
  do i=is,ie
    do j=js,je
      ! i. only want ocean points in this loop
      if (G%mask2dT(i,j) == 0) cycle

      ! ii. Set up vertical domain
      !     TODO: for z*, we may want to set kmt to last layer thicker than (say) 1cm or 1mm
      !           need to update MARBL to handle vanishing layers better
      MARBL_instances%domain%kmt = GV%ke
      set_kmt=.false.
      ! Calculate depth of interface by building up thicknesses from the bottom (top interface is always 0)
      ! MARBL wants this to be positive-down
      zi(GV%ke) = G%bathyT(i,j)
      do k = GV%ke, 1, -1
        ! TODO: if we move this above vertical mixing, use h_old
        dz(k) = h_new(i,j,k)*GV%H_to_Z ! cell thickness
        zc(k) = zi(k) - 0.5 * dz(k)
        zi(k-1) = zi(k) - dz(k)
        ! TODO: better way of handling vanishing layers, e.g. average sediment over bottom 5m
        !       (this is a z* workaround for layers that have vanished at bottom of ocean)
        if ((.not. set_kmt) .and. (dz(k) > 0.01)) then
          set_kmt = .true.
          MARBL_instances%domain%kmt = k
        end if
      enddo

      ! mks -> cgs
      ! zw(1:nz) is bottom cell depth so no element of zw = 0, it is assumed to be top layer depth
      MARBL_instances%domain%zw(:) = zi(1:GV%ke) * cm_per_m
      MARBL_instances%domain%zt(:) = zc(:) * cm_per_m
      MARBL_instances%domain%delta_z(:) = dz(:) * cm_per_m

      ! iii. Load proper column data
      !      * Forcing Fields
      !       These fields are getting the correct data
      if (CS%potemp_ind > 0) MARBL_instances%interior_tendency_forcings(CS%potemp_ind)%field_1d(1,:) = tv%T(i,j,:)
      if (CS%salinity_ind > 0) MARBL_instances%interior_tendency_forcings(CS%salinity_ind)%field_1d(1,:) = tv%S(i,j,:)

      !       This are okay, but need option to read in from file
      !       (Same as dust_dep_ind for surface_flux_forcings)
      if (CS%dustflux_ind > 0) &
        MARBL_instances%interior_tendency_forcings(CS%dustflux_ind)%field_0d(1) = &
            fluxes%MARBL_forcing%dust_flux(i,j) * (kg_m2_s_conversion * g_per_kg * m_per_cm**2)

      !        TODO: Support PAR (currently just using single subcolumn)
      !              (Look for Pen_sw_bnd?)
      if (CS%PAR_col_frac_ind > 0) then
        ! second index is num_subcols, not depth
        MARBL_instances%interior_tendency_forcings(CS%PAR_col_frac_ind)%field_1d(1,:) = 0
        MARBL_instances%interior_tendency_forcings(CS%PAR_col_frac_ind)%field_1d(1,1) = 1
      end if
      if (CS%surf_shortwave_ind > 0) then
        ! second index is num_subcols, not depth
        MARBL_instances%interior_tendency_forcings(CS%surf_shortwave_ind)%field_1d(1,:) = 0
        MARBL_instances%interior_tendency_forcings(CS%surf_shortwave_ind)%field_1d(1,1) = fluxes%sw(i,j)
      end if

      !        TODO: In POP, pressure comes from a function in state_mod.F90; I don't see a similar function here
      !              This formulation is from Levitus 1994, and I think it belongs in MOM_EOS.F90?
      !              Converts depth [m] -> pressure [bars]
      !        NOTE: Andrew recommends using GV%H_to_Pa
      if (CS%pressure_ind > 0) MARBL_instances%interior_tendency_forcings(CS%pressure_ind)%field_1d(1,:) = &
          0.0598088*(exp(-0.025*zc(:)) - 1) + 0.100766*zc(:) + 2.28405e-7*(zc(:)**2)

      if (CS%fesedflux_ind > 0) then
        kmt = MARBL_instances%domain%kmt
        MARBL_instances%interior_tendency_forcings(CS%fesedflux_ind)%field_1d(1,:) = 0.
        call reintegrate_column(CS%fesedflux_nz, &
                                CS%fesedflux_dz(i,j,:) * sum(dz(:)) / (G%bathyT(i,j)), &
                                CS%fesedflux_in(i,j,:) + CS%feventflux_in(i,j,:), &
                                kmt, &
                                dz(1:kmt), &
                                0., &
                                MARBL_instances%interior_tendency_forcings(CS%fesedflux_ind)%field_1d(1,1:kmt))
      end if

      !        TODO: and ability to read these fields from file
      !              also, add constant values to CS
      if (CS%o2_scalef_ind > 0) MARBL_instances%interior_tendency_forcings(CS%o2_scalef_ind)%field_1d(1,:) = 1
      if (CS%remin_scalef_ind > 0) MARBL_instances%interior_tendency_forcings(CS%remin_scalef_ind)%field_1d(1,:) = 1

      !      * Column Tracers
      !        NOTE: POP averages previous two timesteps, should we do that too?
      do m=1,CS%ntr
        MARBL_instances%tracers(m, :) = CS%tr(i,j,:,m)
      end do

      !     * interior tendency saved state
      !       (currently only 3D fields are saved from interior_tendency_compute())
      do m=1,size(MARBL_instances%interior_tendency_saved_state%state)
        MARBL_instances%interior_tendency_saved_state%state(m)%field_3d(:,1) = &
            CS%interior_tendency_saved_state(m)%field_3d(i,j,:)
      end do

      ! iv. Compute interior tendencies in MARBL
      call MARBL_instances%interior_tendency_compute()
      if (MARBL_instances%StatusLog%labort_marbl) then
        call MARBL_instances%StatusLog%log_error_trace("MARBL_instances%interior_tendency_compute()", &
                                                       "MARBL_tracers_column_physics")
        call print_marbl_log(MARBL_instances%StatusLog, G, i, j)
      end if

      ! v. Apply tendencies immediately
      !    First pass - Euler step; if stability issues, we can do something different (subcycle?)
      !    Also, set kmt because we need to avoid using MARBL values below that level
      !    Current plan: treat all values below kmt as equal to the value of the kmt level
      kmt = MARBL_instances%domain%kmt
      do k=1,GV%ke
        if (k <= kmt) then
          CS%tr(i,j,k,:) = CS%tr(i,j,k,:) + G%mask2dT(i,j)*dt*MARBL_instances%interior_tendencies(:, k)
        else
          CS%tr(i,j,k,:) = CS%tr(i,j,k,:) + G%mask2dT(i,j)*dt*MARBL_instances%interior_tendencies(:, kmt)
        end if
      end do

      ! vi. Copy output that MOM6 needs to hold on to
      !     * saved state
      do m=1,size(MARBL_instances%interior_tendency_saved_state%state)
        CS%interior_tendency_saved_state(m)%field_3d(i,j,:) = 0.
        CS%interior_tendency_saved_state(m)%field_3d(i,j,1:kmt) = &
            MARBL_instances%interior_tendency_saved_state%state(m)%field_3d(1:kmt,1)
      end do

      !     * diagnostics
      do m=1,size(MARBL_instances%interior_tendency_diags%diags)
        if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
          CS%interior_tendency_diags(m)%field_2d(i,j) = &
              real(MARBL_instances%interior_tendency_diags%diags(m)%field_2d(1))
        else
          CS%interior_tendency_diags(m)%field_3d(i,j,1:kmt) = &
              real(MARBL_instances%interior_tendency_diags%diags(m)%field_3d(1:kmt,1))
          if (kmt < GV%ke) then
            CS%interior_tendency_diags(m)%field_3d(i,j,kmt+1:GV%ke) = CS%interior_tendency_diags(m)%field_3d(i,j,kmt)
          end if
        end if
      end do
      !     * tendency values themselves (and vertical integrals of them)
      do m=1,CS%ntr
        if (allocated(CS%interior_tendency_out(m)%field_3d)) then
          CS%interior_tendency_out(m)%field_3d(i,j,1:kmt) = G%mask2dT(i,j)*MARBL_instances%interior_tendencies(m,1:kmt)
          if (kmt < GV%ke) then
            CS%interior_tendency_out(m)%field_3d(i,j,kmt+1:GV%ke) = 0.
          end if
        end if

        if (allocated(CS%interior_tendency_out_zint(m)%field_2d)) then
          CS%interior_tendency_out_zint(m)%field_2d(i,j) = 0.
          do k=1,kmt
            CS%interior_tendency_out_zint(m)%field_2d(i,j) = CS%interior_tendency_out_zint(m)%field_2d(i,j) + &
                                                             dz(k) * MARBL_instances%interior_tendencies(m,k)
          end do
          CS%interior_tendency_out_zint(m)%field_2d(i,j) = G%mask2dT(i,j)*CS%interior_tendency_out_zint(m)%field_2d(i,j)
        end if

        if (allocated(CS%interior_tendency_out_zint_100m(m)%field_2d)) then
          CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = 0.
          do k=1,kmt
            if (zi(k) < 100.) then
              CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
                  CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) + dz(k) * MARBL_instances%interior_tendencies(m,k)
            else if (zi(k-1) < 100.) then
              CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
                  CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) + dz(k) * &
                                                                        ((100. - zi(k-1)) / (zi(k) - zi(k-1))) * &
                                                                        MARBL_instances%interior_tendencies(m,k)
            else
              exit
            end if
          end do
          CS%interior_tendency_out_zint_100m(m)%field_2d(i,j) = &
              G%mask2dT(i,j)*CS%interior_tendency_out_zint_100m(m)%field_2d(i,j)
        end if
      end do

    end do
  end do

  ! (5) Post diagnostics from our buffer
  !     i. Interior tendency diagnostics (mix of 2D and 3D)
  !     ii. Interior tendencies themselves
  do m=1,size(CS%interior_tendency_diags)
    if (CS%interior_tendency_diags(m)%id > 0) then
      if (allocated(CS%interior_tendency_diags(m)%field_2d)) then
        call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_2d(:,:), CS%diag)
      else if (allocated(CS%interior_tendency_diags(m)%field_3d)) then
        call post_data(CS%interior_tendency_diags(m)%id, CS%interior_tendency_diags(m)%field_3d(:,:,:), CS%diag)
      else
        write(log_message, "(A, I0, A, I0, A)") "Diagnostic number ", m, " post id ", &
                                                CS%interior_tendency_diags(m)%id," did not allocate 2D or 3D array"
        call MOM_error(FATAL, log_message)
      end if
    end if
  end do
  do m=1,CS%ntr
    if (allocated(CS%interior_tendency_out(m)%field_3d)) &
      call post_data(CS%interior_tendency_out(m)%id, CS%interior_tendency_out(m)%field_3d(:,:,:), CS%diag)
    if (allocated(CS%interior_tendency_out_zint(m)%field_2d)) &
      call post_data(CS%interior_tendency_out_zint(m)%id, CS%interior_tendency_out_zint(m)%field_2d(:,:), CS%diag)
    if (allocated(CS%interior_tendency_out_zint_100m(m)%field_2d)) &
      call post_data(CS%interior_tendency_out_zint_100m(m)%id, CS%interior_tendency_out_zint_100m(m)%field_2d(:,:), &
                     CS%diag)
  end do


end subroutine MARBL_tracers_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function MARBL_tracer_stock(h, stocks, G, GV, CS, names, units, stock_index)
  real, dimension(NIMEM_,NJMEM_,NKMEM_), intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(:),                 intent(out)   :: stocks !< the mass-weighted integrated amount of
                                                            !! each tracer, in kg times concentration units [kg conc].
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  type(MARBL_tracers_CS),                pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to register_MARBL_tracers.
  character(len=*), dimension(:),     intent(out)   :: names !< the names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units !< the units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< the coded index of a specific stock
                                                                   !! being sought.
  integer                                           :: MARBL_tracer_stock   !< Return value: the number of stocks
                                                                   !! calculated here.

! Local variables
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  MARBL_tracer_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="MARBL_tracer_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = 0.0
    do k=1,nz ; do j=js,je ; do i=is,ie
      stocks(m) = stocks(m) + CS%tr(i,j,k,m) * &
                             (G%mask2dT(i,j) * G%areaT(i,j) * h(i,j,k))
    enddo ; enddo ; enddo
    stocks(m) = GV%H_to_kg_m2 * stocks(m)
  enddo
  MARBL_tracer_stock = CS%ntr

end function MARBL_tracer_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine MARBL_tracers_surface_state(state, h, G, CS)
  type(ocean_grid_type),  intent(in)    :: G  !< The ocean's grid structure.
  type(surface),          intent(inout) :: state !< A structure containing fields that
                                              !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                          intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(MARBL_tracers_CS),    pointer       :: CS !< The control structure returned by a previous
                                              !! call to register_MARBL_tracers.

  ! This particular tracer package does not report anything back to the coupler.
  ! The code that is here is just a rough guide for packages that would.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call coupler_type_set_data(CS%tr(:,:,1,m), CS%ind_tr(m), ind_csurf, &
                   state%tr_fields, idim=(/isd, is, ie, ied/), &
                   jdim=(/jsd, js, je, jed/) )
    enddo
  endif

end subroutine MARBL_tracers_surface_state

!> Clean up any allocated memory after the run.
subroutine MARBL_tracers_end(CS)
  type(MARBL_tracers_CS), pointer, intent(inout) :: CS !< The control structure returned by a previous
                                                       !! call to register_MARBL_tracers.
  integer :: m

  call print_marbl_log(MARBL_instances%StatusLog)
  call MARBL_instances%StatusLog%erase()
  call MARBL_instances%shutdown()
  ! TODO: print MARBL timers to stdout as well

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine MARBL_tracers_end

! TODO: some log messages come from a specific grid point, and this routine
!       needs to include the location in the preamble
!> This subroutine writes the contents of the MARBL log using MOM_error(NOTE, ...).
subroutine print_marbl_log(log_to_print, G, i, j)

  use marbl_logging, only : marbl_status_log_entry_type
  use marbl_logging, only : marbl_log_type
  use MOM_coms,      only : PE_here

  class(marbl_log_type),           intent(in) :: log_to_print  !< MARBL log to include in MOM6 logfile
  type(ocean_grid_type), optional, intent(in) :: G             !< The ocean's grid structure
  integer,               optional, intent(in) :: i             !< i of (i,j) index of column providing the log
  integer,               optional, intent(in) :: j             !< j of (i,j) index of column providing the log

  character(len=*), parameter :: subname = 'MARBL_tracers:print_marbl_log'
  character(len=256)          :: message_prefix, message_location, log_message
  type(marbl_status_log_entry_type), pointer :: tmp
  integer :: msg_lev, elem_old

  ! elem_old is used to keep track of whether all messages are coming from the same point
  elem_old = -1
  write(message_prefix, "(A,I0,A)") '(Task ', PE_here(), ')'

  tmp => log_to_print%FullLog
  do while (associated(tmp))
    ! 1) Do I need to write this message? Yes, if all tasks should write this
    !    or if I am master_task
    if ((.not. tmp%lonly_master_writes) .or. is_root_PE()) then
      ! 2) Print message location? (only if ElementInd changed and is positive; requires G)
      if ((present(G)) .and. (tmp%ElementInd .ne. elem_old)) then
        if (tmp%ElementInd .gt. 0) then
          if (present(i) .and. present(j)) then
            write(message_location, "(A,F8.3,A,F7.3,A,I0)") &
                 'Message from (lon, lat) (', G%geoLonT(i,j), ', ', &
                 G%geoLatT(i,j), '). Level: ', tmp%ElementInd
          else
            write(message_location, "(A)") "Grid cell responsible for message is unknown"
          !   i_loc = marbl_col_to_pop_i(tmp%ElementInd, iblock)
          !   j_loc = marbl_col_to_pop_j(tmp%ElementInd, iblock)
          !   write(message_location, "(A,F8.3,A,F7.3,A,I0,A,I0,A)") &
          !        'Message from (lon, lat) (', TLOND(i_loc, j_loc, iblock), &
          !        ", ", TLATD(i_loc, j_loc, iblock), '), which is global (i,j) (', &
          !        this_block%i_glob(i_loc), ', ', this_block%j_glob(j_loc), ')'
          end if ! i,j present
          ! master task does not need prefix
          if (is_root_PE()) then
            write(log_message, "(A)") trim(message_location)
            msg_lev = NOTE
          else
            write(log_message, "(A,1X,A)") trim(message_prefix), trim(message_location)
            msg_lev = WARNING
          end if ! print message prefix?
          call MOM_error(msg_lev, log_message, all_print=.true.)
        end if   ! ElementInd > 0
        elem_old = tmp%ElementInd
      end if     ! ElementInd /= elem_old

      ! 3) Write message from the log
      ! master task does not need prefix
      if (is_root_PE()) then
        write(log_message, "(A)") trim(tmp%LogMessage)
        msg_lev = NOTE
      else
        write(log_message, "(A,1X,A)") trim(message_prefix), trim(tmp%LogMessage)
        msg_lev = WARNING
      end if     ! print message prefix?
      call MOM_error(msg_lev, log_message, all_print=.true.)
    end if       ! write the message?
    tmp => tmp%next
  end do

  if (log_to_print%labort_marbl) then
    call MOM_error(WARNING, 'ERROR reported from MARBL library', all_print=.true.)
    call MOM_error(FATAL, 'Stopping in ' // subname)
  end if

end subroutine print_marbl_log

! TODO: fix the comment below
!> \namespace MARBL_tracers
!!
!!    This file contains an example of the code that is needed to set
!!  up and use a set (in this case two) of dynamically passive tracers
!!  for diagnostic purposes.  The tracers here are dye tracers which
!!  are set to 1 within the geographical region specified. The depth
!!  which a tracer is set is determined by calculating the depth from
!!  the seafloor upwards through the column.

end module MARBL_tracers
