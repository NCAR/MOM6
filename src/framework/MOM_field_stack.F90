!> This module contains utilities for maintaning stacks of model fields.
!! The current implementation handles 3 dimensional fields.
module MOM_field_stack

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_error, FATAL, MOM_get_verbosity, MOM_mesg

implicit none ; private

public field_stack_init, field_stack_end
public field_stack_push, field_stack_pop
public field_stack_peek, field_stack_drop

!> The field_stack type
type, public :: field_stack_type ; private
  real, dimension(:,:,:,:), pointer :: field => NULL() !< field values stored in stack
  integer :: top = 0                                   !< index of top of stack, zero for an empty stack
  character(:), allocatable :: name                    !< stack name to include in informational messages
end type field_stack_type

contains

!> Initialize a field_stack object
subroutine field_stack_init(field_stack_obj, is, ie, js, je, nz, name)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being initialized
  integer,                   intent(in) :: is              !< The start index to allocate for the 1st dimension
  integer,                   intent(in) :: ie              !< The end index to allocate for the 1st dimension
  integer,                   intent(in) :: js              !< The start index to allocate for the 2nd dimension
  integer,                   intent(in) :: je              !< The end index to allocate for the 2nd dimension
  integer,                   intent(in) :: nz              !< The size to allocate for the 3rd dimension
  character(len=*),          intent(in) :: name            !< stack name to include in informational messages

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_init"

  if (associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has already been initialized")

  allocate(field_stack_obj%field(is:ie, js:je, nz, 0))
  field_stack_obj%top = 0
  field_stack_obj%name = trim(name)

  call write_stack_info_msg(sub_name, field_stack_obj, "init")

end subroutine field_stack_init


!> Deallocate memory associated with field_stack object
subroutine field_stack_end(field_stack_obj, info_msg)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being operated on
  character(len=*),          intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_end"

  if (.not. associated(field_stack_obj%field)) return

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  deallocate(field_stack_obj%field)
  field_stack_obj%field => NULL()
  field_stack_obj%top = 0
  deallocate(field_stack_obj%name)

end subroutine field_stack_end


!> Increase memory available in stack, preserving current stack values
subroutine field_stack_grow(field_stack_obj, info_msg)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being operated on
  character(len=*),          intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_grow"
  real, dimension(:,:,:,:), pointer :: field_ptr_tmp !< temporary pointer
  integer :: is, ie, js, je, nz, stack_size !< bounds of stack's field

  if (.not. associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has not been initialized")

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  field_ptr_tmp => field_stack_obj%field

  is = lbound(field_stack_obj%field, 1) ; ie = ubound(field_stack_obj%field, 1)
  js = lbound(field_stack_obj%field, 2) ; je = ubound(field_stack_obj%field, 2)
  nz = size(field_stack_obj%field, 3)
  stack_size = size(field_stack_obj%field, 4)

  stack_size = stack_size + 1

  allocate(field_stack_obj%field(is:ie, js:je, nz, stack_size))

  if (field_stack_obj%top > 0) then
    field_stack_obj%field(:,:,:,1:field_stack_obj%top) = field_ptr_tmp(:,:,:,1:field_stack_obj%top)
  endif

  deallocate(field_ptr_tmp)

end subroutine field_stack_grow


!> Push field values to stack
subroutine field_stack_push(field_stack_obj, field, info_msg)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being operated on
  real, dimension(:,:,:),    intent(in) :: field           !< field values to be pushed to stack
  character(len=*),          intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_push"
  integer :: i, j, k
  integer :: i_offset, j_offset

  if (.not. associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has not been initialized")

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  if (field_stack_obj%top == size(field_stack_obj%field, 4)) &
      call field_stack_grow(field_stack_obj, info_msg)

  ! verify that field has same shape as stack field
  if ((size(field, 1) /= size(field_stack_obj%field, 1)) .or. &
      (size(field, 2) /= size(field_stack_obj%field, 2)) .or. &
      (size(field, 3) /= size(field_stack_obj%field, 3))) then
    call MOM_error(FATAL, &
        sub_name // ": field argument shape doesn't match stack field shape")
  endif

  field_stack_obj%top = field_stack_obj%top + 1

  i_offset = lbound(field_stack_obj%field, 1) - 1
  j_offset = lbound(field_stack_obj%field, 2) - 1
  do k=1,size(field,3) ; do j=1,size(field,2) ; do i=1,size(field,1)
    field_stack_obj%field(i+i_offset,j+j_offset,k,field_stack_obj%top) = field(i,j,k)
  enddo ; enddo ; enddo

end subroutine field_stack_push


!> Pop field values from stack
subroutine field_stack_pop(field_stack_obj, field, info_msg)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being operated on
  real, dimension(:,:,:),   intent(out) :: field           !< where popped field values are stored
  character(len=*),          intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_pop"
  integer :: i, j, k
  integer :: i_offset, j_offset

  if (.not. associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has not been initialized")

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  if (field_stack_obj%top == 0) call MOM_error(FATAL, &
      sub_name // ": attempting to pop from empty stack")

  ! verify that field has same shape as stack field
  if ((size(field, 1) /= size(field_stack_obj%field, 1)) .or. &
      (size(field, 2) /= size(field_stack_obj%field, 2)) .or. &
      (size(field, 3) /= size(field_stack_obj%field, 3))) then
    call MOM_error(FATAL, &
        sub_name // ": field argument shape doesn't match stack field shape")
  endif

  i_offset = lbound(field_stack_obj%field, 1) - 1
  j_offset = lbound(field_stack_obj%field, 2) - 1
  do k=1,size(field,3) ; do j=1,size(field,2) ; do i=1,size(field,1)
    field(i,j,k) = field_stack_obj%field(i+i_offset,j+j_offset,k,field_stack_obj%top)
  enddo ; enddo ; enddo

  field_stack_obj%top = field_stack_obj%top - 1

end subroutine field_stack_pop


!> Associate pointer with top item on stack
subroutine field_stack_peek(field_stack_obj, field_ptr, info_msg)
  type(field_stack_type), intent(in) :: field_stack_obj !< field_stack object being inspected
  real, dimension(:,:,:), pointer    :: field_ptr       !< pointer to be associated with top of stack
  character(len=*),       intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_peek"

  if (.not. associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has not been initialized")

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  if (field_stack_obj%top == 0) call MOM_error(FATAL, &
      sub_name // ": attempting to peek at empty stack")

  field_ptr => field_stack_obj%field(:,:,:,field_stack_obj%top)

end subroutine field_stack_peek


!> Drop top field values from stack
subroutine field_stack_drop(field_stack_obj, info_msg)
  type(field_stack_type), intent(inout) :: field_stack_obj !< field_stack object being operated on
  character(len=*),          intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(*), parameter :: sub_name = "field_stack_drop"
  integer :: i, j, k
  integer :: i_offset, j_offset

  if (.not. associated(field_stack_obj%field)) call MOM_error(FATAL, &
      sub_name // ": field_stack_obj has not been initialized")

  call write_stack_info_msg(sub_name, field_stack_obj, info_msg)

  if (field_stack_obj%top == 0) call MOM_error(FATAL, &
      sub_name // ": attempting to drop from empty stack")

  field_stack_obj%top = field_stack_obj%top - 1

end subroutine field_stack_drop


!> Write informational message using MOM_mesg
subroutine write_stack_info_msg(caller, field_stack_obj, info_msg)
  character(len=*),       intent(in) :: caller          !< name of subroutine writing message
  type(field_stack_type), intent(in) :: field_stack_obj !< field_stack object being referred to
  character(len=*),       intent(in) :: info_msg        !< informational message to write

  ! Local variables
  character(2) :: top_i2 !< character version of stack top value
  character(2) :: size_i2 !< character version of stack size value

  ! Only construct full informational message if verbosity is sufficiently high
  if (MOM_get_verbosity() >= 4) then
    write(top_i2, '(I2)') field_stack_obj%top
    write(size_i2, '(I2)') size(field_stack_obj%field, 4)
    call MOM_mesg( &
        caller // ": " // trim(field_stack_obj%name) // " " // trim(info_msg) // ", top=" &
        // top_i2 // ", size=" // size_i2, 4)
  endif

end subroutine write_stack_info_msg

end module MOM_field_stack
