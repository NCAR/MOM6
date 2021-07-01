!> A non-functioning template of the public structures provided through MARBL interface
module marbl_interface_public_types

    implicit none
    private ! Only want a few types to be public

    !> A non-functioning template of MARBL diagnostic type
    type :: marbl_single_diagnostic_type
        character(len=0)                  :: long_name  !< dummy name
        character(len=0)                  :: short_name  !< dummy name
        character(len=0)                  :: units  !< dummy units
        character(len=0)                  :: vertical_grid  !< dummy grid
        logical                           :: compute_now  !< dummy flag
        logical                           :: ltruncated_vertical_extent  !< dummy flag
        integer                           :: ref_depth  !< dummy depth
        real, allocatable, dimension(:)   :: field_2d  !< dummy field
        real, allocatable, dimension(:,:) :: field_3d  !< dummy field
    end type marbl_single_diagnostic_type

    !> A non-functioning template of MARBL diagnostic type
    type, public :: marbl_diagnostics_type
        type(marbl_single_diagnostic_type), dimension(:), pointer :: diags => NULL()  !< dummy point
    end type marbl_diagnostics_type

    !> A non-functioning template of MARBL saved state type
    type :: marbl_single_saved_state_type
        integer :: rank  !< dummy rank
        character(len=0) :: short_name  !< dummy name
        character(len=0) :: units  !< dummy units
        character(len=0) :: vertical_grid  !< dummy grid
        real, allocatable :: field_2d(:)  !< dummy field
        real, allocatable :: field_3d(:,:)  !< dummy field
    end type marbl_single_saved_state_type

    !> A non-functioning template of MARBL saved state type
    type, public :: marbl_saved_state_type
        integer :: saved_state_cnt  !< dummy counter
        type(marbl_single_saved_state_type), dimension(:), pointer :: state => NULL()  !< dummy pointer
    end type marbl_saved_state_type

    !> A non-functioning template of MARBL forcing metadata type
    type :: marbl_forcing_fields_metadata_type
        character(len=0) :: varname  !< dummy name
    end type marbl_forcing_fields_metadata_type

    !> A non-functioning template of MARBL forcing type
    type, public :: marbl_forcing_fields_type
        type(marbl_forcing_fields_metadata_type) :: metadata  !< dummy metadata
        real, pointer :: field_0d(:)   => NULL()  !< dummy pointer
        real, pointer :: field_1d(:,:) => NULL()  !< dummy pointer
    end type marbl_forcing_fields_type

    !> A non-functioning template of MARBL tracer metadata type
    type, public :: marbl_tracer_metadata_type
        character(len=0) :: short_name  !< dummy name
        character(len=0) :: long_name  !< dummy name
        character(len=0) :: units  !< dummy units
    end type marbl_tracer_metadata_type

    !> A non-functioning template of MARBL domain type
    type, public :: marbl_domain_type
        integer :: kmt  !< dummy index
        real, allocatable :: zt(:)  !< dummy depths
        real, allocatable :: zw(:)  !< dummy depths
        real, allocatable :: delta_z(:)  !< dummy thicknesses
    end type marbl_domain_type

end module marbl_interface_public_types