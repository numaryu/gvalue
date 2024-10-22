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
     procedure :: print_results
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

       write(6,'(/1x,a,i0,a/)') '=====> WORK ', iwork, ': done'
    end do
  contains

    subroutine init_orbital(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
      character (len=100) :: file_orbital
      integer :: norbital, ngeneration
      real :: number_density
      logical :: ex
      namelist /param_orbital/ file_orbital, ngeneration
      namelist /orbital_dim/ norbital
      namelist /param_medium/ number_density

      ! default values
      norbital = 1
      ngeneration = 6

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, param_orbital)
      close(unit)
      write(unit_stdout, param_orbital)

      inquire(file=trim(file_orbital), exist = ex)
      if (.not.ex) stop 'orbital file does not exist!'

      call get_unused_unit(unit)
      open(unit, file=trim(file_orbital))
      read(unit, orbital_dim)
      read(unit, param_medium)
      close(unit)
      write(unit_stdout, orbital_dim)
      write(unit_stdout, param_medium)

      worker%medium = orbital(norbital, ngeneration, file_orbital)

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

  subroutine print_results(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use class_grid, only: grid
    class(work) :: self
    integer :: iwork
    integer :: io, ie
    integer :: unit
    character (len=100) :: file

    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'.dat'

       call get_unused_unit(unit)
       open(unit, file=trim(file))

       write(unit,'("#")')
       write(unit,'("#",a)') ' Yield:'
       write(unit,'("#",4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
       do io = 1, self%worker(iwork)%medium%number
          write(unit,'("#",15x,i5,3(8x,f12.4))') io, &
               self%worker(iwork)%medium%yield_ionize(io), &
               self%worker(iwork)%medium%yield_singlet(io), &
               self%worker(iwork)%medium%yield_triplet(io)
       end do
       if (self%worker(iwork)%medium%number > 1) then
          write(unit,'("#",a20,3(8x,f12.4))') 'sum', &
               sum(self%worker(iwork)%medium%yield_ionize), &
               sum(self%worker(iwork)%medium%yield_singlet), &
               sum(self%worker(iwork)%medium%yield_triplet)
       end if

       write(unit,'("#")')
       write(unit,'("#",a)') ' G-value:'
       write(unit,'("#",4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
       do io = 1, self%worker(iwork)%medium%number
          write(unit,'("#",15x,i5,3(8x,f12.4))') io, &
               self%worker(iwork)%medium%gvalue_ionize(io), &
               self%worker(iwork)%medium%gvalue_singlet(io), &
               self%worker(iwork)%medium%gvalue_triplet(io)
       end do
       if (self%worker(iwork)%medium%number > 1) then
          write(unit,'("#",a20,3(8x,f12.4))') 'sum', &
               sum(self%worker(iwork)%medium%gvalue_ionize), &
               sum(self%worker(iwork)%medium%gvalue_singlet), &
               sum(self%worker(iwork)%medium%gvalue_triplet)
       end if
       write(unit,'("#")')

       write(unit,'("#",11(1x,a20))') 'Energy', 'Stopping Power', 'Degradation (sum)', &
            'Total Cross Sec. i', 'Total Cross Sec. s', 'Total Cross Sec. t', &
            'Platzman i', 'Platzman s', 'Platzman t', 'Mean Free Path', 'Range'

       do ie = 1, self%worker(iwork)%egrid%number
          write(unit,'(1x,11(1x,e20.12))') self%worker(iwork)%egrid%val(ie), &
               self%worker(iwork)%medium%stop_power(ie), &
               sum(self%worker(iwork)%medium%degradation_gen(:, ie)), &
               self%worker(iwork)%medium%total_cross_section_ionize(ie), &
               self%worker(iwork)%medium%total_cross_section_singlet(ie), &
               self%worker(iwork)%medium%total_cross_section_triplet(ie), &
               sum(self%worker(iwork)%medium%platzman_ionize_orbital(:, ie)), &
               sum(self%worker(iwork)%medium%platzman_singlet_orbital(:, ie)), &
               sum(self%worker(iwork)%medium%platzman_triplet_orbital(:, ie)), &
               self%worker(iwork)%medium%mean_free_path(ie), &
               self%worker(iwork)%medium%range(ie)
       end do
       close(unit)

       unit = unit_stdout

       write(unit,*)
       write(unit,'(1x,a)') ' Yield:'
       write(unit,'(1x,4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
       do io = 1, self%worker(iwork)%medium%number
          write(unit,'(1x,15x,i5,3(8x,f12.4))') io, &
               self%worker(iwork)%medium%yield_ionize(io), &
               self%worker(iwork)%medium%yield_singlet(io), &
               self%worker(iwork)%medium%yield_triplet(io)
       end do
       if (self%worker(iwork)%medium%number > 1) then
          write(unit,'(1x,a20,3(8x,f12.4))') 'sum', &
               sum(self%worker(iwork)%medium%yield_ionize), &
               sum(self%worker(iwork)%medium%yield_singlet), &
               sum(self%worker(iwork)%medium%yield_triplet)
       end if

       write(unit,*)
       write(unit,'(1x,a)') ' G-value:'
       write(unit,'(1x,4(a20))') 'Orbital', 'ionize', 'singlet', 'triplet'
       do io = 1, self%worker(iwork)%medium%number
          write(unit,'(1x,15x,i5,3(8x,f12.4))') io, &
               self%worker(iwork)%medium%gvalue_ionize(io), &
               self%worker(iwork)%medium%gvalue_singlet(io), &
               self%worker(iwork)%medium%gvalue_triplet(io)
       end do
       if (self%worker(iwork)%medium%number > 1) then
          write(unit,'(1x,a20,3(8x,f12.4))') 'sum', &
               sum(self%worker(iwork)%medium%gvalue_ionize), &
               sum(self%worker(iwork)%medium%gvalue_singlet), &
               sum(self%worker(iwork)%medium%gvalue_triplet)
       end if
    end do

    call print_results_degradation(self)

  end subroutine print_results

  subroutine print_results_degradation(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use class_grid, only: grid
    class(work) :: self
    integer :: iwork
    integer :: ie, igen, ngen
    integer :: unit
    character (len=100) :: file

    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'_degradation.dat'

       ngen = size(self%worker(iwork)%medium%degradation_gen, 1)

       call get_unused_unit(unit)
       open(unit, file=trim(file))

       write(unit,'("#",2(1x,a20))', advance='no') 'Energy', 'Degradation (sum)'
       do igen = 1, ngen
          write(unit,'(1x, a16, 1x, i3)', advance='no') 'Degradation ', igen
       end do
       write(unit,*)

       do ie = 1, self%worker(iwork)%egrid%number
          write(unit,'(1x,2(1x,e20.12))', advance='no') &
               self%worker(iwork)%egrid%val(ie), &
               sum(self%worker(iwork)%medium%degradation_gen(:, ie))
          do igen = 1, ngen
             write(unit,'(1x,e20.12)', advance='no') &
                  self%worker(iwork)%medium%degradation_gen(igen, ie)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

  end subroutine print_results_degradation

end module class_work
