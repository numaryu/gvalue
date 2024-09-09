module mod_vars
  implicit none
  private

  public :: ydeg, ydeg_total
  public :: spow
  
  public :: init_vars, finish_vars
  
  integer, parameter :: ndeg = 6
  
  real, allocatable :: ydeg(:,:), ydeg_total(:)
  real, allocatable :: spow(:)

  logical :: initialized = .false.
  logical :: debug = .false.
  
contains

  subroutine init_vars
    use mod_grid, only: negrid => negrid_compat
    if (initialized) return
    if (.not.allocated(ydeg)) allocate(ydeg(negrid,ndeg))
    if (.not.allocated(ydeg_total)) allocate(ydeg_total(negrid))
    if (.not.allocated(spow)) allocate(spow(negrid))
    initialized = .true.
  end subroutine init_vars

  subroutine finish_vars
    if (.not.initialized) return
    if (allocated(ydeg)) deallocate(ydeg)
    if (allocated(ydeg_total)) deallocate(ydeg_total)
    if (allocated(spow)) deallocate(spow)
    initialized = .false.
  end subroutine finish_vars
  
end module mod_vars
