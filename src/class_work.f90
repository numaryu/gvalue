module class_work
  use class_orbital, only: orbital
  use class_grid, only: grid
  implicit none
  private

  public :: work

  type work
     private
     integer :: number
     type(subwork), pointer :: worker(:) => null()
   contains
     final :: destroy
     procedure :: execute
  end type work

  type subwork
     character (len=100) :: name
     type(orbital) :: medium
     type(grid) :: egrid
  end type subwork

  ! declare constructor
  interface work
     module procedure init_work
  end interface work

contains

  type(work) function init_work()
    integer :: iarg, l
    character (len=100) :: name
    logical :: ex
    init_work%number = command_argument_count()
    if (.not.associated(init_work%worker)) allocate(init_work%worker(init_work%number))

    ! parse arguments
    do iarg = 1, init_work%number
       call get_command_argument(iarg, name)
       l = len_trim(name)
       if (l > 3) then
          if (name(l-2:l) == ".in") then
             name = name(1:l-3)
          end if
       end if
       init_work%worker(iarg)%name = trim(name)

       inquire(file=trim(init_work%worker(iarg)%name)//'.in', exist = ex)
       if (.not.ex) stop 'input file does not exist!'
    end do

  end function init_work

  subroutine destroy(self)
    type(work), intent(in out) :: self
    if (associated(self%worker)) nullify(self%worker)
  end subroutine destroy

  subroutine execute(self)
    class(work) :: self
    integer :: iwork
    do iwork = 1, self%number
       write(6,'(/1x,a,i0,a,a/)') '=====> WORK ', iwork, ': ',self%worker(iwork)%name

       call init_orbital(self%worker(iwork))
       call init_grid(self%worker(iwork))

       call self%worker(iwork)%medium%init_orbital_vars(self%worker(iwork)%egrid%number)
       call self%worker(iwork)%medium%calculate_stopping_power(self%worker(iwork)%egrid)
       call self%worker(iwork)%medium%calculate_degradation(self%worker(iwork)%egrid)
       call self%worker(iwork)%medium%calculate_yield(self%worker(iwork)%egrid)

       call self%worker(iwork)%medium%print_results(self%worker(iwork)%name, self%worker(iwork)%egrid)

       write(6,'(/1x,a,i0,a/)') '=====> WORK ', iwork, ': done'
    end do
  contains

    subroutine init_orbital(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
      character (len=100) :: file_orbital
      integer :: norbital
      real :: number_density
      namelist /param_orbital/ file_orbital
      namelist /orbital_dim/ norbital
      namelist /param_medium/ number_density

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, param_orbital)
      close(unit)
      write(unit_stdout, param_orbital)

      call get_unused_unit(unit)
      open(unit, file=trim(file_orbital))
      read(unit, orbital_dim)
      read(unit, param_medium)
      close(unit)
      write(unit_stdout, orbital_dim)
      write(unit_stdout, param_medium)

      worker%medium = orbital(norbital, file_orbital)

      worker%medium%number_density = number_density
    end subroutine init_orbital

    subroutine init_grid(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
      ! number of division in energy grid
      ! nediv means T of i+nediv is a half of T of i
      integer :: nediv
      ! egrid_max, egrid_min are the maximum, minimum energy
      real :: egrid_max, egrid_min

      namelist /grid_param/ nediv, egrid_max, egrid_min

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, grid_param)
      close(unit)
      write(unit_stdout, grid_param)

      worker%egrid = grid(nediv, egrid_max, egrid_min)
    end subroutine init_grid

  end subroutine execute

end module class_work
