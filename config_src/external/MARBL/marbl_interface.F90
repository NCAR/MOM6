!> A non-functioning template of the MARBL interface
module marbl_interface

    use MOM_error_handler,            only : MOM_error, FATAL
    use marbl_logging,                only : marbl_log_type
    use marbl_interface_public_types, only : marbl_forcing_fields_type
    use marbl_interface_public_types, only : marbl_tracer_metadata_type
    use marbl_interface_public_types, only : marbl_saved_state_type
    use marbl_interface_public_types, only : marbl_diagnostics_type
    use marbl_interface_public_types, only : marbl_domain_type
    implicit none
    private ! Only want marbl_interface_class to be public, not supporting functions

    !> A non-functioning template of the MARBL_interface class
    !!
    !> All variables are dummy representations of actual members of the real marbl_interface_class
    !! that are used in the MARBL tracer routines.
    type, public :: marbl_interface_class
        type(marbl_log_type) :: StatusLog  !< dummy log
        type(marbl_forcing_fields_type), allocatable :: surface_flux_forcings(:)  !< dummy forcing array
        type(marbl_forcing_fields_type), allocatable :: interior_tendency_forcings(:)  !< dummy forcing array
        type(marbl_tracer_metadata_type), allocatable :: tracer_metadata(:)  !< dummy metadata array
        type(marbl_domain_type) :: domain  !< dummy domain
        type(marbl_saved_state_type) :: surface_flux_saved_state  !< dummy saved state
        type(marbl_saved_state_type) :: interior_tendency_saved_state  !< dummy saved state
        type(marbl_diagnostics_type) :: interior_tendency_diags  !< dummy diagnostics
        type(marbl_diagnostics_type) :: surface_flux_diags  !< dummy diagnostics
        real, allocatable :: tracers(:,:)  !< dummy tracer array
        real, allocatable :: tracers_at_surface(:,:)  !< dummy tracer surface array
        real, allocatable :: surface_fluxes(:,:)  !< dummy fluxes
        real, allocatable :: interior_tendencies(:,:)  !< dummy tendencies
       contains
        procedure, public  :: init  !< dummy routine
        procedure, public  :: surface_flux_compute  !< dummy surface flux routine
        procedure, public  :: interior_tendency_compute  !< dummy interior tendency routine
        procedure, public  :: shutdown  !< dummy shutdown routine
        procedure, public :: put_setting  !< dummy put_setting routine
    end type marbl_interface_class

    !> Error message that appears if the dummy interface is called
    character(len=*), parameter :: error_msg = "MOM6 built the MARBL stubs rather than the full library"

contains

    !> Dummy version of MARBL's put_setting() function
    subroutine put_setting(self, str_in)
        class(marbl_interface_class), intent(in) :: self
        character(len=*),             intent(in) :: str_in

        call MOM_error(FATAL, error_msg)
    end subroutine put_setting

    !> Dummy version of MARBL's init() function
    subroutine init(self,                  &
        gcm_num_levels,                    &
        gcm_num_PAR_subcols,               &
        gcm_num_elements_surface_flux,     &
        gcm_delta_z,                       &
        gcm_zw,                            &
        gcm_zt,                            &
        lgcm_has_global_ops)

        class(marbl_interface_class), intent(inout) :: self
        integer,                      intent(in)    :: gcm_num_levels
        integer,                      intent(in)    :: gcm_num_PAR_subcols
        integer,                      intent(in)    :: gcm_num_elements_surface_flux
        real,                         intent(in)    :: gcm_delta_z(gcm_num_levels)
        real,                         intent(in)    :: gcm_zw(gcm_num_levels)
        real,                         intent(in)    :: gcm_zt(gcm_num_levels)
        logical,                      intent(in)    :: lgcm_has_global_ops

        call MOM_error(FATAL, error_msg)
    end subroutine init

    !> Dummy version of MARBL's surface_flux_compute() function
    subroutine surface_flux_compute(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine surface_flux_compute

    !> Dummy version of MARBL's interior_tendency_compute() function
    subroutine interior_tendency_compute(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine interior_tendency_compute

    !> Dummy version of MARBL's shutdown() function
    subroutine shutdown(self)

        class(marbl_interface_class), intent(inout) :: self

        call MOM_error(FATAL, error_msg)

    end subroutine shutdown

end module marbl_interface