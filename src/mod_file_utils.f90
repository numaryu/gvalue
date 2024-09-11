module mod_file_utils
  implicit none
  private
  public :: get_unused_unit
  public :: unit_stdin, unit_stdout, unit_stderr

  integer, parameter :: unit_stdin = 5
  integer, parameter :: unit_stdout = 6
  integer, parameter :: unit_stderr = 0

contains

  subroutine get_unused_unit(unit)
    implicit none
    !> A new unit not associated with an open file
    integer, intent (out) :: unit
    logical :: od
    unit = 50
    do
       inquire (unit=unit, opened=od)
       if (.not.od) return
       unit = unit + 1
    end do
  end subroutine get_unused_unit
  
end module mod_file_utils
