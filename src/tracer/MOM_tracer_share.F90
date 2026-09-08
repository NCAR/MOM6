! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This module contains miscellaneous utilities usuable by tracer modules.
module MOM_tracer_share

use MOM_coms,            only : broadcast
use MOM_error_handler,   only : is_root_PE, MOM_error, FATAL
use MOM_io,              only : MOM_infra_file, MOM_field

implicit none ; private

public :: MOM_tracer_read_lines
public :: MOM_IO_handles_find_name

contains


!> This subroutine reads lines from filename into file_contents.
subroutine MOM_tracer_read_lines(filename, file_contents)
  character(len=*), intent(in)               :: filename          !< Name of file being read.
  character(len=:), intent(out), allocatable :: file_contents(:)  !< Variable where file contents are written.

  ! local variables
  character(len=256) :: fileline      ! line from filename
  integer            :: max_linelen   ! maximum line length
  integer            :: read_iter     ! count lines on read_iter=1, store them on read_iter=2
  integer            :: line_ind      ! which line in file currently being read
  integer            :: readunit      ! I/O unit for reading filename
  integer            :: ios           ! status from I/O call
  integer            :: line_cnt      ! number of lines in file

  ! read contents of filename on root PE
  if (is_root_PE()) then
    max_linelen = 0
    do read_iter = 1, 2
      line_ind = 1
      open(newunit=readunit, file=filename, form='formatted', status='old', iostat=ios)
      if (ios /= 0) call MOM_error(FATAL, "error opening " // trim(filename))
      do
        read(readunit,'(a)', iostat=ios) fileline
        if (ios /= 0) then
          if (is_iostat_end(ios)) exit ! done reading lines from file
          call MOM_error(FATAL, "error reading from " // trim(filename))
        endif
        if (read_iter == 1) then
          max_linelen = max(max_linelen, len_trim(fileline))
        else
          file_contents(line_ind) = trim(fileline)
        endif
        line_ind = line_ind + 1
      enddo
      close(readunit)
      line_cnt = line_ind - 1
      if (line_cnt == 0) call MOM_error(FATAL, trim(filename) // " appears to be empty")
      if (read_iter == 1) allocate(character(len=max_linelen) :: file_contents(line_cnt))
    enddo
    close(readunit)
  endif

  ! broadcast results to other PEs
  call broadcast(max_linelen)
  call broadcast(line_cnt)
  if (.not. is_root_PE()) allocate(character(len=max_linelen) :: file_contents(line_cnt))
  call broadcast(file_contents, max_linelen)

end subroutine MOM_tracer_read_lines

!> This function returns the first index of the IO_handle whose
!! corresponding file has a variable with varname=name.
!! If no file has a variable with varname=name 0 is returned.
function MOM_IO_handles_find_name(IO_handles, name) result(file_ind)

  type(MOM_infra_file), intent(inout) :: IO_handles(:) !< handles to files being searched
  character(len=*), intent(in)        :: name          !< name being searched for in IO_handles
  integer                             :: file_ind      !< returned index into IO_handles

  ! local variables
  type(MOM_field), allocatable :: fields(:)
  character(len=80) :: file_varname
  integer :: nvar, var_ind
  logical :: var_found

  var_found = .false.
  do file_ind = 1, size(IO_handles)
    call IO_handles(file_ind)%get_file_info(nvar=nvar)
    allocate(fields(nvar))
    call IO_handles(file_ind)%get_file_fields(fields)
    do var_ind = 1, nvar
      call IO_handles(file_ind)%get_field_atts(fields(var_ind), name=file_varname)
      if (file_varname == name) then
        var_found = .true.
        exit
      endif
    enddo
    deallocate(fields)
    if (var_found) exit
  enddo
  if (.not. var_found) file_ind = 0

end function MOM_IO_handles_find_name

end module MOM_tracer_share
