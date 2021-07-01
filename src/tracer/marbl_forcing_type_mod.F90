!> This module provides a common datatype to provide forcing for MARBL tracers
!! regardless of driver
module marbl_forcing_type_mod

!! This module exists to house code used by multiple drivers in config_src/
!! for passing forcing fields to MARBL
!! (This comment can go in the wiki on the NCAR fork?)

use MOM_diag_mediator,        only : safe_alloc_ptr, diag_ctrl, register_diag_field, post_data
use MOM_time_manager,         only : time_type
use MOM_error_handler,        only : MOM_error, WARNING
use MOM_file_parser,          only : get_param, log_param, param_file_type
use MOM_grid,                 only : ocean_grid_type
use MOM_unit_scaling,         only : unit_scale_type
use MOM_interpolate,          only : init_external_field, time_interp_external
use tracer_forcing_utils_mod, only : forcing_timeseries_dataset
use tracer_forcing_utils_mod, only : forcing_timeseries_set_time_type_vars
use tracer_forcing_utils_mod, only : map_model_time_to_forcing_time

implicit none ; private

#include <MOM_memory.h>

!> Data type used to store diagnostic index returned from register_diag_field()
!! For the forcing fields that can be written via post_data()
type, private :: marbl_forcing_diag_ids
  integer :: atm_fine_dust   !< Atmospheric fine dust component of dust_flux
  integer :: atm_coarse_dust !< Atmospheric coarse dust component of dust_flux
  integer :: atm_bc          !< Atmospheric black carbon component of iron_flux
  integer :: ice_dust        !< Sea-ice dust component of dust_flux
  integer :: ice_bc          !< Sea-ice black carbon component of iron_flux
  ! River fluxes
  integer :: no3_riv_flux          !< NO3 riverine flux
  integer :: po4_riv_flux          !< PO4 riverine flux
  integer :: don_riv_flux          !< DON riverine flux
  integer :: donr_riv_flux         !< DONr riverine flux
  integer :: dop_riv_flux          !< DOP riverine flux
  integer :: dopr_riv_flux         !< DOPr riverine flux
  integer :: sio3_riv_flux         !< SiO3 riverine flux
  integer :: fe_riv_flux           !< Fe riverine flux
  integer :: doc_riv_flux          !< DOC riverine flux
  integer :: docr_riv_flux         !< DOCr riverine flux
  integer :: alk_riv_flux          !< ALK riverine flux
  integer :: alk_alt_co2_riv_flux  !< ALK (alternate CO2) riverine flux
  integer :: dic_riv_flux          !< DIC riverine flux
  integer :: dic_alt_co2_riv_flux  !< DIC (alternate CO2) riverine flux
end type marbl_forcing_diag_ids

!> Contains pointers to the forcing fields needed to drive MARBL
type, public :: marbl_forcing_type
  real, pointer, dimension(:,:) :: noy_dep => NULL() !< NOy Deposition [R Z T-1 ~> kgN m-2 s-1]
  real, pointer, dimension(:,:) :: nhx_dep => NULL() !< NHx Deposition [R Z T-1 ~> kgN m-2 s-1]
  real, pointer, dimension(:,:) :: dust_flux => NULL() !< Flux of dust into the ocean [m2 m-2]
  real, pointer, dimension(:,:) :: iron_flux => NULL() !< Flux of dust into the ocean [m2 m-2]
  real, pointer, dimension(:,:) :: ice_fraction => NULL() !< Fraction of ocean cell under seaice [m2 m-2]
  real, pointer, dimension(:,:) :: u10_sqr => NULL() !< 10m wind speed squared [L2 T-2 ~> m2 s-2]

  real, pointer, dimension(:,:) :: no3_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: po4_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: sio3_riv_flux => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: fe_riv_flux   => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: alk_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: doc_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: docr_riv_flux => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: don_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: donr_riv_flux => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: dop_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: dopr_riv_flux => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: dic_riv_flux  => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: alk_alt_co2_riv_flux => NULL() !< [mmol / m^2 / s]
  real, pointer, dimension(:,:) :: dic_alt_co2_riv_flux => NULL() !< [mmol / m^2 / s]
end type marbl_forcing_type

!> Control structure for this module
type, public :: marbl_forcing_CS
  logical :: read_ndep                      !< If true, use nitrogen deposition supplied from an input file.
                                            !! This is temporary, we will always read NDEP
  character(len=200) :: ndep_file           !< If read_ndep, then this is the file from which to read
  type(forcing_timeseries_dataset) :: riv_flux_dataset !< File and time axis information for river fluxes
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.

  real    :: dust_ratio_to_fe_bioavail_frac !< TODO: Add description
  real    :: fe_bioavail_frac_offset        !< TODO: Add description
  real    :: atm_fe_to_bc_ratio             !< TODO: Add description
  real    :: atm_bc_fe_bioavail_frac        !< TODO: Add description
  real    :: seaice_fe_to_bc_ratio          !< TODO: Add description
  real    :: seaice_bc_fe_bioavail_frac     !< TODO: Add description
  real    :: iron_frac_in_atm_fine_dust     !< Fraction of fine dust from the atmosphere that is iron
  real    :: iron_frac_in_atm_coarse_dust   !< Fraction of coarse dust from the atmosphere that is iron
  real    :: iron_frac_in_seaice_dust       !< Fraction of dust from the sea ice that is iron

  type(marbl_forcing_diag_ids) :: diag_ids  !< used for registering and posting some MARBL forcing fields as diagnostics


  integer :: id_noydep   = -1     !< id number for time_interp_external.
  integer :: id_nhxdep   = -1     !< id number for time_interp_external.
  integer :: id_din_riv  = -1     !< id number for time_interp_external.
  integer :: id_don_riv  = -1     !< id number for time_interp_external.
  integer :: id_dip_riv  = -1     !< id number for time_interp_external.
  integer :: id_dop_riv  = -1     !< id number for time_interp_external.
  integer :: id_dsi_riv  = -1     !< id number for time_interp_external.
  integer :: id_dfe_riv  = -1     !< id number for time_interp_external.
  integer :: id_dic_riv  = -1     !< id number for time_interp_external.
  integer :: id_alk_riv  = -1     !< id number for time_interp_external.
  integer :: id_doc_riv  = -1     !< id number for time_interp_external.

end type marbl_forcing_CS

!> Contains pointers to IOB fields that are used to compute MARBL forcing
type, public :: marbl_ice_ocean_boundary_type
  real, pointer, dimension(:,:) :: atm_fine_dust_flux   => NULL() !< Fine dust flux from atmosphere [kg/m^2/s]
  real, pointer, dimension(:,:) :: atm_coarse_dust_flux => NULL() !< Coarse dust flux from atmosphere [kg/m^2/s]
  real, pointer, dimension(:,:) :: seaice_dust_flux     => NULL() !< Dust flux from seaice [kg/m^2/s]
  real, pointer, dimension(:,:) :: atm_bc_flux          => NULL() !< Black carbon flux from atmosphere [kg/m^2/s]
  real, pointer, dimension(:,:) :: seaice_bc_flux       => NULL() !< Black carbon flux from seaice [kg/m^2/s]
  real, pointer, dimension(:,:) :: ice_fraction         => NULL() !< Fraction of ocn covered with ice
  real, pointer, dimension(:,:) :: u10_sqr              => NULL() !< 10m wind speed squared (m^2/s^2)
end type marbl_ice_ocean_boundary_type

public :: marbl_forcing_init
public :: marbl_forcing_type_init
public :: marbl_iob_allocate
public :: convert_marbl_IOB_to_forcings

contains

  subroutine marbl_forcing_init(G, param_file, diag, day, inputdir, CS)
    type(ocean_grid_type),           intent(in)    :: G           !< The ocean's grid structure
    type(param_file_type),           intent(in)    :: param_file  !< A structure to parse for run-time parameters
    type(diag_ctrl), target,         intent(in)    :: diag        !< Structure used to regulate diagnostic output.
    type(time_type), target,         intent(in)    :: day         !< Time of the start of the run.
    character(len=*),                intent(in)    :: inputdir    !< Directory containing input files
    type(marbl_forcing_CS), pointer, intent(inout) :: CS          !< A pointer that is set to point to control
                                                                  !! structure for MARBL forcing

    character(len=40)  :: mdl = "MOM_forcing_type"  ! This module's name.
    integer :: riv_flux_file_start_year
    integer :: riv_flux_file_end_year
    integer :: riv_flux_file_data_ref_year
    integer :: riv_flux_file_model_ref_year
    integer :: riv_flux_forcing_year

    if (associated(CS)) then
      call MOM_error(WARNING, "marbl_forcing_init called with an associated "// &
                              "control structure.")
      return
    endif

    allocate(CS)
    CS%diag => diag

    call get_param(param_file, mdl, "DUST_RATIO_TO_FE_BIOAVAIL_FRAC", CS%dust_ratio_to_fe_bioavail_frac, &
    "TODO: Add description", default=1./170.)
    call get_param(param_file, mdl, "FE_BIOAVAIL_FRAC_OFFSET", CS%fe_bioavail_frac_offset, &
        "TODO: Add description", default=0.01)
    call get_param(param_file, mdl, "ATM_FE_TO_BC_RATIO", CS%atm_fe_to_bc_ratio, &
        "TODO: Add description", default=1.)
    call get_param(param_file, mdl, "ATM_BC_FE_BIOAVAIL_FRAC", CS%atm_bc_fe_bioavail_frac, &
        "TODO: Add description", default=0.06)
    call get_param(param_file, mdl, "SEAICE_FE_TO_BC_RATIO", CS%seaice_fe_to_bc_ratio, &
        "TODO: Add description", default=1.)
    call get_param(param_file, mdl, "SEAICE_BC_FE_BIOAVAIL_FRAC", CS%seaice_bc_fe_bioavail_frac, &
        "TODO: Add description", default=0.06)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_FINE_DUST", CS%iron_frac_in_atm_fine_dust, &
        "Fraction of fine dust from the atmosphere that is iron", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_ATM_COARSE_DUST", CS%iron_frac_in_atm_coarse_dust, &
        "Fraction of coarse dust from the atmosphere that is iron", default=0.035)
    call get_param(param_file, mdl, "IRON_FRAC_IN_SEAICE_DUST", CS%iron_frac_in_seaice_dust, &
        "Fraction of dust from sea ice that is iron", default=0.035)

    call get_param(param_file, mdl, "READ_NDEP", CS%read_ndep, &
        "If true, use nitrogen deposition supplied from "//&
        "an input file", default=.true.)
    if (CS%read_ndep) then
      ! TODO: we only want to read this variable in when running with MARBL
      call get_param(param_file, mdl, "NDEP_FILE", CS%ndep_file, &
            "The file in which the nitrogen deposition is found in "//&
            "variables NOy_deposition and NHx_deposition.", &
            default='ndep_ocn_1850_w_nhx_emis_MOM_tx0.66v1_c210222.nc')
      if (scan(CS%ndep_file,'/') == 0) then
        ! CS%ndep_file = trim(inputdir) // trim(CS%ndep_file)
        CS%ndep_file = trim('/glade/work/mlevy/cesm_inputdata/') // trim(CS%ndep_file)
        call log_param(param_file, mdl, "INPUTDIR/NDEP_FILE", CS%ndep_file)
      end if
      CS%id_noydep = init_external_field(CS%ndep_file, 'NDEP_NOy_month', domain=G%Domain%mpp_domain)
      CS%id_nhxdep = init_external_field(CS%ndep_file, 'NDEP_NHx_month', domain=G%Domain%mpp_domain)
    end if

    ! ** River fluxes
    call get_param(param_file, mdl, "RIV_FLUX_FILE", CS%riv_flux_dataset%file_name, &
                   "The file in which the river fluxes can be found", &
                   default="riv_nut.gnews_gnm.rx1_to_tx0.66v1_nnsm_e1000r300_190315.20210405.nc")
    ! call get_param(param_file, mdl, "RIV_FLUX_OFFSET_YEAR", CS%riv)
    if (scan(CS%riv_flux_dataset%file_name,'/') == 0) then
      ! CS%riv_flux_dataset%file_name = trim(inputdir) // trim(CS%riv_flux_dataset%file_name)
      CS%riv_flux_dataset%file_name = trim('/glade/work/mlevy/cesm_inputdata/') // trim(CS%riv_flux_dataset%file_name)
      call log_param(param_file, mdl, "INPUTDIR/RIV_FLUX_FILE", CS%riv_flux_dataset%file_name)
    end if
    call get_param(param_file, mdl, "RIV_FLUX_L_TIME_VARYING", CS%riv_flux_dataset%l_time_varying, &
                   ".true. for time-varying forcing, .false. for static forcing", default=.false.)
    if (CS%riv_flux_dataset%l_time_varying) then
      call get_param(param_file, mdl, "RIV_FLUX_FILE_START_YEAR", riv_flux_file_start_year, &
                     "Time coordinate of earliest date in RIV_FLUX_FILE", default=1900)
      call get_param(param_file, mdl, "RIV_FLUX_FILE_END_YEAR", riv_flux_file_end_year, &
                     "Time coordinate of earliest date in RIV_FLUX_FILE", default=1999)
      call get_param(param_file, mdl, "RIV_FLUX_FILE_DATA_REF_YEAR", riv_flux_file_data_ref_year, &
                     "Time coordinate of latest date in RIV_FLUX_FILE", default=1900)
      call get_param(param_file, mdl, "RIV_FLUX_FILE_MODEL_REF_YEAR", riv_flux_file_model_ref_year, &
                     "Time coordinate of latest date in RIV_FLUX_FILE",  default=1900)
    else
      call get_param(param_file, mdl, "RIV_FLUX_FORCING_YEAR", riv_flux_forcing_year, &
                     "Year from RIV_FLUX_FILE to use for forcing",  default=1900)
    end if
    call forcing_timeseries_set_time_type_vars(riv_flux_file_start_year, &
                                               riv_flux_file_end_year, &
                                               riv_flux_file_data_ref_year, &
                                               riv_flux_file_model_ref_year, &
                                               riv_flux_forcing_year, &
                                               CS%riv_flux_dataset)

    CS%id_din_riv = init_external_field(CS%riv_flux_dataset%file_name, 'din_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_don_riv = init_external_field(CS%riv_flux_dataset%file_name, 'don_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_dip_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dip_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_dop_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dop_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_dsi_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dsi_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_dfe_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dfe_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_dic_riv = init_external_field(CS%riv_flux_dataset%file_name, 'dic_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_alk_riv = init_external_field(CS%riv_flux_dataset%file_name, 'alk_riv_flux', domain=G%Domain%mpp_domain)
    CS%id_doc_riv = init_external_field(CS%riv_flux_dataset%file_name, 'doc_riv_flux', domain=G%Domain%mpp_domain)

    ! Register diagnostic fields for outputing forcing values
    CS%diag_ids%atm_fine_dust = register_diag_field("ocean_model", &
                                                    "ATM_FINE_DUST_FLUX_CPL", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "ATM_FINE_DUST_FLUX from cpl", &
                                                    "kg/m^2/s")
    CS%diag_ids%atm_coarse_dust = register_diag_field("ocean_model", &
                                                      "ATM_COARSE_DUST_FLUX_CPL", &
                                                      CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                      day, &
                                                      "ATM_COARSE_DUST_FLUX from cpl", &
                                                      "kg/m^2/s")
    CS%diag_ids%atm_bc = register_diag_field("ocean_model", &
                                             "ATM_BLACK_CARBON_FLUX_CPL", &
                                             CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                             day, &
                                             "ATM_BLACK_CARBON_FLUX from cpl", &
                                             "kg/m^2/s")

    CS%diag_ids%ice_dust = register_diag_field("ocean_model", &
                                               "SEAICE_DUST_FLUX_CPL", &
                                               CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                               day, &
                                               "SEAICE_DUST_FLUX from cpl", &
                                               "kg/m^2/s")
    CS%diag_ids%ice_bc = register_diag_field("ocean_model", &
                                             "SEAICE_BLACK_CARBON_FLUX_CPL", &
                                             CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                             day, &
                                             "SEAICE_BLACK_CARBON_FLUX from cpl", &
                                             "kg/m^2/s")

    CS%diag_ids%no3_riv_flux = register_diag_field("ocean_model", &
                                                   "NO3_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Nitrate Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%po4_riv_flux = register_diag_field("ocean_model", &
                                                   "PO4_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Phosphate Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%don_riv_flux = register_diag_field("ocean_model", &
                                                   "DON_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Nitrogen Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%donr_riv_flux = register_diag_field("ocean_model", &
                                                    "DONR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DON Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%dop_riv_flux = register_diag_field("ocean_model", &
                                                   "DOP_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Phosphorus Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%dopr_riv_flux = register_diag_field("ocean_model", &
                                                    "DOPR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DOP Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%sio3_riv_flux = register_diag_field("ocean_model", &
                                                    "SiO3_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Dissolved Inorganic Silicate Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%fe_riv_flux = register_diag_field("ocean_model", &
                                                  "Fe_RIV_FLUX", &
                                                  CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                  day, &
                                                  "Dissolved Inorganic Iron Riverine Flux", &
                                                  "mmol/m^3 m/s")
    CS%diag_ids%doc_riv_flux = register_diag_field("ocean_model", &
                                                   "DOC_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Organic Carbon Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%docr_riv_flux = register_diag_field("ocean_model", &
                                                    "DOCR_RIV_FLUX", &
                                                    CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                    day, &
                                                    "Refractory DOC Riverine Flux", &
                                                    "mmol/m^3 m/s")
    CS%diag_ids%alk_riv_flux = register_diag_field("ocean_model", &
                                                   "ALK_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Alkalinity Riverine Flux", &
                                                   "meq/m^3 m/s")
    CS%diag_ids%alk_alt_co2_riv_flux = register_diag_field("ocean_model", &
                                                          "ALK_ALT_CO2_RIV_FLUX", &
                                                          CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                          day, &
                                                          "Alkalinity Riverine Flux, Alternative CO2", &
                                                          "meq/m^3 m/s")
    CS%diag_ids%dic_riv_flux = register_diag_field("ocean_model", &
                                                   "DIC_RIV_FLUX", &
                                                   CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                   day, &
                                                   "Dissolved Inorganic Carbon Riverine Flux", &
                                                   "mmol/m^3 m/s")
    CS%diag_ids%dic_alt_co2_riv_flux = register_diag_field("ocean_model", &
                                                          "DIC_ALT_CO2_RIV_FLUX", &
                                                          CS%diag%axesT1, & ! T=> tracer grid? 1 => no vertical grid
                                                          day, &
                                                          "Dissolved Inorganic Carbon Riverine Flux, Alternative CO2", &
                                                          "mmol/m^3 m/s")

  end subroutine marbl_forcing_init

  subroutine marbl_forcing_type_init(isd,ied,jsd,jed,MARBL_forcing)
    integer,                           intent(in)    :: isd            !< start of i-indices for current block
    integer,                           intent(in)    :: ied            !< end of i-indices for current block
    integer,                           intent(in)    :: jsd            !< start of j-indices for current block
    integer,                           intent(in)    :: jed            !< end of j-indices for current block
    type(marbl_forcing_type), pointer, intent(inout) :: MARBL_forcing  !< MARBL-specific forcing fields

    if (associated(MARBL_forcing)) then
      call MOM_error(WARNING, "marbl_forcing_type_init called with an associated "// &
                              "marbl forcing structure.")
      return
    endif

    allocate(MARBL_forcing)

    ! Fields from coupler
    call safe_alloc_ptr(MARBL_forcing%noy_dep,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%nhx_dep,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%dust_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%iron_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%ice_fraction,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%u10_sqr,isd,ied,jsd,jed)

    ! Fields read from file
    call safe_alloc_ptr(MARBL_forcing%no3_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%po4_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%sio3_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%fe_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%alk_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%doc_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%docr_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%don_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%donr_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%dop_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%dopr_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%dic_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%alk_alt_co2_riv_flux,isd,ied,jsd,jed)
    call safe_alloc_ptr(MARBL_forcing%dic_alt_co2_riv_flux,isd,ied,jsd,jed)

  end subroutine marbl_forcing_type_init

  subroutine marbl_iob_allocate(isc, iec, jsc, jec, MARBL_IOB)

    integer,                                      intent(in) :: isc           !< start of i-indices for current block
    integer,                                      intent(in) :: iec           !< end of i-indices for current block
    integer,                                      intent(in) :: jsc           !< start of j-indices for current block
    integer,                                      intent(in) :: jec           !< end of j-indices for current block
    type(marbl_ice_ocean_boundary_type), pointer, intent(inout) :: MARBL_IOB  !< MARBL-specific ice-ocean boundary type

    if (associated(MARBL_IOB)) then
      call MOM_error(WARNING, "marbl_iob_allocate called with an associated "// &
                              "ice-ocean boundary structure.")
      return
    endif

    allocate(MARBL_IOB)
    allocate(MARBL_IOB% atm_fine_dust_flux (isc:iec,jsc:jec),  &
             MARBL_IOB% atm_coarse_dust_flux (isc:iec,jsc:jec),&
             MARBL_IOB% seaice_dust_flux (isc:iec,jsc:jec),    &
             MARBL_IOB% atm_bc_flux (isc:iec,jsc:jec),         &
             MARBL_IOB% seaice_bc_flux (isc:iec,jsc:jec),      &
             MARBL_IOB% ice_fraction (isc:iec,jsc:jec),        &
             MARBL_IOB% u10_sqr (isc:iec,jsc:jec))

    MARBL_IOB%atm_fine_dust_flux(:,:)   = 0.0
    MARBL_IOB%atm_coarse_dust_flux(:,:) = 0.0
    MARBL_IOB%seaice_dust_flux(:,:)     = 0.0
    MARBL_IOB%atm_bc_flux(:,:)          = 0.0
    MARBL_IOB%seaice_bc_flux(:,:)       = 0.0
    MARBL_IOB%ice_fraction              = 0.0
    MARBL_IOB%u10_sqr                   = 0.0

  end subroutine marbl_iob_allocate

  subroutine convert_marbl_IOB_to_forcings(MARBL_IOB, Time, G, US, i0, j0, MARBL_forcing, CS)

    type(marbl_ice_ocean_boundary_type), pointer, intent(in)    :: MARBL_IOB      !< MARBL-specific ice-ocean boundary
                                                                                  !! type
    type(time_type),                              intent(in)    :: Time           !< The time of the fluxes, used for
                                                                                  !! interpolating the salinity to the
                                                                                  !! right time, when it is being
                                                                                  !! restored.
    type(ocean_grid_type),                        intent(in)    :: G              !< The ocean's grid structure
    type(unit_scale_type),                        intent(in)    :: US             !< A dimensional unit scaling type
    integer,                                      intent(in)    :: i0             !< i index offset
    integer,                                      intent(in)    :: j0             !< j index offset
    type(marbl_forcing_type),                     intent(inout) :: MARBL_forcing  !< MARBL-specific forcing fields
    type(marbl_forcing_CS), pointer,              intent(inout) :: CS             !< A pointer that is set to point to
                                                                                  !! control structure for MARBL forcing

    ! These are Fortran parameters in POP
    real, parameter :: DONriv_refract = 0.1
    real, parameter :: DOCriv_refract = 0.2
    real, parameter :: DOPriv_refract = 0.025

    real, dimension(SZI_(G),SZJ_(G)) :: time_varying_data  !< The field read in from forcing file with time dimension
    type(time_type) :: Time_riv_flux  !< For reading river flux fields, we use a modified version of Time
    integer :: i, j, is, ie, js, je
    real :: ndep_conversion          !< Combination of unit conversion factors for rescaling
                                     !! nitrogen deposition [g(N) m-2 s-1 ~> mol L-2 T-2]
    real :: kg_m2_s_conversion       !< A combination of unit conversion factors for rescaling
                                     !! mass fluxes [R Z s m2 kg-1 T-1 ~> 1].

    is   = G%isc   ; ie   = G%iec    ; js   = G%jsc   ; je   = G%jec
    ndep_conversion = (1/14.) * ((US%L_to_m)**2 * US%T_to_s)
    kg_m2_s_conversion = US%kg_m2s_to_RZ_T

    ! Post fields from coupler to diagnostics
    ! TODO: units from diag register are incorrect; we should be converting these in the cap, I think
    if (CS%diag_ids%atm_fine_dust > 0) &
      call post_data(CS%diag_ids%atm_fine_dust, &
                     kg_m2_s_conversion * MARBL_IOB%atm_fine_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_coarse_dust > 0) &
      call post_data(CS%diag_ids%atm_coarse_dust, &
                     kg_m2_s_conversion * MARBL_IOB%atm_coarse_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%atm_bc > 0) &
      call post_data(CS%diag_ids%atm_bc, kg_m2_s_conversion * MARBL_IOB%atm_bc_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_dust > 0) &
      call post_data(CS%diag_ids%ice_dust, kg_m2_s_conversion * MARBL_IOB%seaice_dust_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))
    if (CS%diag_ids%ice_bc > 0) &
      call post_data(CS%diag_ids%ice_bc, kg_m2_s_conversion * MARBL_IOB%seaice_bc_flux(is-i0:ie-i0,js-j0:je-j0), &
                     CS%diag, mask=G%mask2dT(is:ie,js:je))

    do j=js,je ; do i=is,ie
      if (associated(MARBL_IOB%atm_fine_dust_flux)) then
        ! TODO: MARBL wants g/cm^2/s; we should convert to RZ_T in ocn_cap_methods then back to MARBL units here
        MARBL_forcing%dust_flux(i,j) = (G%mask2dT(i,j) * kg_m2_s_conversion) * &
                                       (MARBL_IOB%atm_fine_dust_flux(i-i0,j-j0) + &
                                        MARBL_IOB%atm_coarse_dust_flux(i-i0,j-j0) + &
                                        MARBL_IOB%seaice_dust_flux(i-i0,j-j0))
      end if

      if (associated(MARBL_IOB%atm_bc_flux)) then
        ! Contribution of atmospheric dust to iron flux
        MARBL_forcing%iron_flux(i,j) = (CS%fe_bioavail_frac_offset * &
                                 (CS%iron_frac_in_atm_fine_dust * MARBL_IOB%atm_fine_dust_flux(i-i0,j-j0) + &
                                  CS%iron_frac_in_atm_coarse_dust * MARBL_IOB%atm_coarse_dust_flux(i-i0,j-j0)))
        ! Contribution of atmospheric black carbon to iron flux
        MARBL_forcing%iron_flux(i,j) = MARBL_forcing%iron_flux(i,j) + (CS%atm_bc_fe_bioavail_frac * &
                                 (CS%atm_fe_to_bc_ratio * MARBL_IOB%atm_bc_flux(i-i0,j-j0)))
        ! Contribution of seaice dust to iron flux
        MARBL_forcing%iron_flux(i,j) = MARBL_forcing%iron_flux(i,j) + (CS%fe_bioavail_frac_offset * &
                                 (CS%iron_frac_in_seaice_dust * MARBL_IOB%seaice_dust_flux(i-i0,j-j0)))
        ! Contribution of seaice black carbon to iron flux
        MARBL_forcing%iron_flux(i,j) = MARBL_forcing%iron_flux(i,j) + (CS%seaice_bc_fe_bioavail_frac * &
                                 (CS%seaice_fe_to_bc_ratio * MARBL_IOB%seaice_bc_flux(i-i0,j-j0)))
        ! Unit conversion
        MARBL_forcing%iron_flux(i,j) = (G%mask2dT(i,j) * kg_m2_s_conversion) * MARBL_forcing%iron_flux(i,j)
      end if

      if (associated(MARBL_IOB%ice_fraction)) then
        MARBL_forcing%ice_fraction(i,j) = G%mask2dT(i,j) * MARBL_IOB%ice_fraction(i-i0,j-j0)
      end if

      if (associated(MARBL_IOB%u10_sqr)) then
        MARBL_forcing%u10_sqr(i,j) = G%mask2dT(i,j) * US%m_s_to_L_T**2 * MARBL_IOB%u10_sqr(i-i0,j-j0)
      end if
    enddo; enddo

    if (CS%read_ndep) then
      call time_interp_external(CS%id_noydep,Time,time_varying_data)
      MARBL_forcing%noy_dep = ndep_conversion * time_varying_data
      call time_interp_external(CS%id_nhxdep,Time,time_varying_data)
      MARBL_forcing%nhx_dep = ndep_conversion * time_varying_data
      do j=js,je ; do i=is,ie
        MARBL_forcing%noy_dep(i,j) = G%mask2dT(i,j) * MARBL_forcing%noy_dep(i,j)
        MARBL_forcing%nhx_dep(i,j) = G%mask2dT(i,j) * MARBL_forcing%nhx_dep(i,j)
      enddo; enddo
    else
      ! This is temporary while I test the file we created
      MARBL_forcing%noy_dep(:,:) = 0.
      MARBL_forcing%nhx_dep(:,:) = 0.
    endif

    ! River fluxes
    MARBL_forcing%alk_riv_flux(:,:) = 0.
    MARBL_forcing%alk_alt_co2_riv_flux(:,:) = 0.

    ! DIN river flux affects NO3, ALK, and ALK_ALT_CO2
    Time_riv_flux = map_model_time_to_forcing_time(Time, CS%riv_flux_dataset)

    call time_interp_external(CS%id_din_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%no3_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)
    MARBL_forcing%alk_riv_flux(:,:) = MARBL_forcing%alk_riv_flux(:,:) - time_varying_data(:,:)
    MARBL_forcing%alk_alt_co2_riv_flux(:,:) = MARBL_forcing%alk_alt_co2_riv_flux(:,:) - time_varying_data(:,:)

    call time_interp_external(CS%id_dip_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%po4_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

    call time_interp_external(CS%id_don_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%don_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DONriv_refract) * &
                                      time_varying_data(:,:)
    MARBL_forcing%donr_riv_flux(:,:) = G%mask2dT(:,:) * DONriv_refract * &
                                       time_varying_data(:,:)

    call time_interp_external(CS%id_dop_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%dop_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DOPriv_refract) * &
                                      time_varying_data(:,:)
    MARBL_forcing%dopr_riv_flux(:,:) = G%mask2dT(:,:) * DOPriv_refract * &
                                       time_varying_data(:,:)

    call time_interp_external(CS%id_dsi_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%sio3_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

    call time_interp_external(CS%id_dfe_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%fe_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

    call time_interp_external(CS%id_dic_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%dic_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)
    MARBL_forcing%dic_alt_co2_riv_flux(:,:) = G%mask2dT(:,:) * time_varying_data(:,:)

    call time_interp_external(CS%id_alk_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%alk_riv_flux(:,:) = MARBL_forcing%alk_riv_flux(:,:) + time_varying_data(:,:)
    MARBL_forcing%alk_alt_co2_riv_flux(:,:) = MARBL_forcing%alk_alt_co2_riv_flux(:,:) + time_varying_data(:,:)

    call time_interp_external(CS%id_doc_riv,Time_riv_flux,time_varying_data)
    MARBL_forcing%doc_riv_flux(:,:) = G%mask2dT(:,:) * (1. - DOCriv_refract) * time_varying_data(:,:)
    MARBL_forcing%docr_riv_flux(:,:) = G%mask2dT(:,:) * DOCriv_refract * time_varying_data(:,:)

    ! Post to diags
    if (CS%diag_ids%no3_riv_flux > 0) &
      call post_data(CS%diag_ids%no3_riv_flux, MARBL_forcing%no3_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%po4_riv_flux > 0) &
      call post_data(CS%diag_ids%po4_riv_flux, MARBL_forcing%po4_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%don_riv_flux > 0) &
      call post_data(CS%diag_ids%don_riv_flux, MARBL_forcing%don_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%donr_riv_flux > 0) &
      call post_data(CS%diag_ids%donr_riv_flux, MARBL_forcing%donr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dop_riv_flux > 0) &
      call post_data(CS%diag_ids%dop_riv_flux, MARBL_forcing%dop_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dopr_riv_flux > 0) &
      call post_data(CS%diag_ids%dopr_riv_flux, MARBL_forcing%dopr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%sio3_riv_flux > 0) &
      call post_data(CS%diag_ids%sio3_riv_flux, MARBL_forcing%sio3_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%fe_riv_flux > 0) &
      call post_data(CS%diag_ids%fe_riv_flux, MARBL_forcing%fe_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%doc_riv_flux > 0) &
      call post_data(CS%diag_ids%doc_riv_flux, MARBL_forcing%doc_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%docr_riv_flux > 0) &
      call post_data(CS%diag_ids%docr_riv_flux, MARBL_forcing%docr_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%alk_riv_flux > 0) &
      call post_data(CS%diag_ids%alk_riv_flux, MARBL_forcing%alk_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%alk_alt_co2_riv_flux > 0) &
      call post_data(CS%diag_ids%alk_alt_co2_riv_flux, MARBL_forcing%alk_alt_co2_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dic_riv_flux > 0) &
      call post_data(CS%diag_ids%dic_riv_flux, MARBL_forcing%dic_riv_flux(:,:), CS%diag)
    if (CS%diag_ids%dic_alt_co2_riv_flux > 0) &
      call post_data(CS%diag_ids%dic_alt_co2_riv_flux, MARBL_forcing%dic_alt_co2_riv_flux(:,:), CS%diag)

  end subroutine convert_marbl_IOB_to_forcings

end module marbl_forcing_type_mod
