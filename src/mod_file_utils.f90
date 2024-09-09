module mod_file_utils
  implicit none
  private
  public :: get_unused_unit

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
