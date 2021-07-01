!> This module provides a common framework for interacting with time_interp_external.
!! Its main goal is to allow for offsets between the time axis in the forcing file
!! and the model run itself. It also provides a mechanism for using specific dates
!! as earliest / latest dates to read from file (so use data_start until Time + m2d_offset
!! catches up, and use data_end once Time + m2doffset exceeds it)
module tracer_forcing_utils_mod

use MOM_time_manager, only : time_type, real_to_time, operator(+), operator(<), operator(>)

!> Data type used to store information about forcing datasets that are time series
!! E.g. how do we align the data in the model with the time axis in the file?
type, public :: forcing_timeseries_dataset
    character(len=200) :: file_name  !< name of file containing river flux forcing
    logical :: l_time_varying        !< .true. => forcing is dependent on model time, .false. => static forcing
    ! logical :: l_FMS_modulo        !< .true. => let FMS handle determining time level to read (e.g. for climatologies)
    type(time_type) :: data_forcing  !< convert data_forcing_year to time type
    type(time_type) :: data_start    !< convert data_start_year to time type
    type(time_type) :: data_end      !< convert data_end_year to time type
    type(time_type) :: m2d_offset    !< add to model time to get data time
end type forcing_timeseries_dataset

public :: forcing_timeseries_set_time_type_vars
public :: map_model_time_to_forcing_time

contains

    !> Set time_type variables in forcing_timeseries_dataset type based on integer input
    !! TODO: make this part of forcing_timeseries_dataset class if OO is okay in MOM6?
    subroutine forcing_timeseries_set_time_type_vars(data_start_year, data_end_year, data_ref_year, &
                                                     model_ref_year, data_forcing_year, forcing_dataset)

        integer,                          intent(in)    :: data_start_year    !< first year of data to read
                                                                              !! (this is ignored for static forcing)
        integer,                          intent(in)    :: data_end_year      !< last year of data to read
                                                                              !! (this is ignored for static forcing)
        integer,                          intent(in)    :: data_ref_year      !< for time-varying forcing, align
                                                                              !! data_ref_year in file with
                                                                              !! model_ref_year in model
        integer,                          intent(in)    :: model_ref_year     !< for time-varying forcing, align
                                                                              !! data_ref_year in file with
                                                                              !! model_ref_year in model
        integer,                          intent(in)    :: data_forcing_year  !< for static forcing, read file at this
                                                                              !! date (this is ignored for time-varying
                                                                              !! forcing)
        type(forcing_timeseries_dataset), intent(inout) :: forcing_dataset    !< information about forcing file

        if (forcing_dataset%l_time_varying) then
            forcing_dataset%data_start = real_to_time(year_to_sec(data_start_year))
            forcing_dataset%data_end = real_to_time(year_to_sec(data_end_year))
            forcing_dataset%m2d_offset = real_to_time(year_to_sec(data_ref_year - model_ref_year))
        else
            forcing_dataset%data_forcing = real_to_time(year_to_sec(data_forcing_year))
        end if

    end subroutine forcing_timeseries_set_time_type_vars

    !> If necessary, apply an offset to convert from model time to forcing time and then
    !! ensure result is within acceptable bounds
    function map_model_time_to_forcing_time(Time, forcing_dataset)

        type(time_type),                  intent(in)  :: Time             !< Model time
        type(forcing_timeseries_dataset), intent(in)  :: forcing_dataset  !< information about forcing file
        type(time_type) :: map_model_time_to_forcing_time                 !< time to read forcing file

        if (forcing_dataset%l_time_varying) then
            map_model_time_to_forcing_time = Time + forcing_dataset%m2d_offset
            ! If Time + offset is not between data_start and data_end, use whichever of those values is closer
            if (map_model_time_to_forcing_time < forcing_dataset%data_start) &
                map_model_time_to_forcing_time = forcing_dataset%data_start
            if (map_model_time_to_forcing_time > forcing_dataset%data_end) &
                map_model_time_to_forcing_time = forcing_dataset%data_end
        else
            map_model_time_to_forcing_time = forcing_dataset%data_forcing
        end if

    end function map_model_time_to_forcing_time

    !> real_to_time converts from seconds since 0000-01-01 to time_type so we need to convert from years -> seconds
    function year_to_sec(year)

        integer, intent(in) :: year
        real :: year_to_sec

        year_to_sec = 86400. * 365. * real(year)

    end function year_to_sec

end module tracer_forcing_utils_mod