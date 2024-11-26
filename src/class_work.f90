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
     integer :: nmedia, ngeneration
     type(orbital), pointer :: medium(:) => null()
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
    integer :: iwork, imedia
    do iwork = 1, self%number
       write(6,'(/1x,a,i0,a,a/)') '=====> WORK ', iwork, ': ',self%worker(iwork)%name

       call init_subwork(self%worker(iwork))

       do imedia = 1, self%worker(iwork)%nmedia
          call self%worker(iwork)%medium(imedia)%init_orbital_vars(self%worker(iwork)%egrid%number, &
               self%worker(iwork)%ngeneration)
          call self%worker(iwork)%medium(imedia)%calculate_stopping_power(self%worker(iwork)%egrid)
          call self%worker(iwork)%medium(imedia)%calculate_degradation(self%worker(iwork)%egrid, &
               self%worker(iwork)%ngeneration)
          call self%worker(iwork)%medium(imedia)%calculate_yield(self%worker(iwork)%egrid)
       end do

       write(6,'(/1x,a,i0,a/)') '=====> WORK ', iwork, ': done'
    end do
  contains

    subroutine init_subwork(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
      integer :: imedia
      integer :: ngeneration, nmedia
      namelist /param_calc/ ngeneration, nmedia

      ! default values
      ngeneration = 6
      nmedia = 1

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, param_calc)
      close(unit)
      write(unit_stdout, param_calc)

      worker%ngeneration = ngeneration
      worker%nmedia = nmedia
      if (.not.associated(worker%medium)) allocate(worker%medium(nmedia))

      call init_orbital(worker)
      call init_grid(worker)

    end subroutine init_subwork

    subroutine init_orbital(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
      integer, parameter :: nmedia_max = 10
      integer :: imedia
      character (len=100) :: file_medium(worker%nmedia)
      real :: number_density(worker%nmedia)
      integer :: norbital
      character (len=100) :: name
      logical :: ex
      namelist /param_medium/ file_medium, number_density
      namelist /param_orbital/ norbital, name

      ! default values
      number_density = 1.
      norbital = 1
      file_medium = ''

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, param_medium)
      close(unit)
      write(unit_stdout, param_medium)

      do imedia = 1, worker%nmedia
         inquire(file=trim(file_medium(imedia)), exist = ex)
         if (.not.ex) stop 'orbital file does not exist!'

         call get_unused_unit(unit)
         open(unit, file=trim(file_medium(imedia)))
         read(unit, param_orbital)
         close(unit)
         write(unit_stdout, param_orbital)

         worker%medium(imedia) = orbital(norbital, name, &
              file_medium(imedia), number_density(imedia))
      end do

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

      namelist /param_grid/ nediv, egrid_max, egrid_min

      ! default values
      nediv = 40
      egrid_max = 1.e5
      egrid_min = 1.e0

      call get_unused_unit(unit)
      open(unit, file=trim(worker%name)//'.in')
      read(unit, param_grid)
      close(unit)
      write(unit_stdout, param_grid)

      worker%egrid = grid(nediv, egrid_max, egrid_min)
    end subroutine init_grid

  end subroutine execute

  subroutine print_results(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use class_grid, only: grid
    class(work) :: self
    integer :: iwork, imedia
    integer :: io, ie
    integer :: unit
    character (len=100) :: file

    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'.dat'

       call get_unused_unit(unit)
       open(unit, file=trim(file))

       do imedia = 1, self%worker(iwork)%nmedia

          write(unit,'("#")')
          write(unit,'("#",1x 2a)') 'Media: ', trim(self%worker(iwork)%medium(imedia)%name)
          write(unit,'("#")')
          call print_yield(unit, imedia)
          call print_gvalue(unit, imedia)

          write(unit,'("#")')
          write(unit,'("#",13(1x,a20))') 'Energy', 'Stopping Power', 'Degradation (sum)', &
               'Total Cross Sec. i', 'Total Cross Sec. s', 'Total Cross Sec. t', &
               'Platzman i', 'Platzman s', 'Platzman t', 'Mean Free Path', 'Range i', 'Range s', 'Range t'

          do ie = 1, self%worker(iwork)%egrid%number
             write(unit,'(1x,13(1x,e20.12))') self%worker(iwork)%egrid%val(ie), &
                  self%worker(iwork)%medium(1)%stop_power(ie), &
                  sum(self%worker(iwork)%medium(1)%degradation_gen(:, ie)), &
                  self%worker(iwork)%medium(1)%total_cross_section_ionize(ie), &
                  self%worker(iwork)%medium(1)%total_cross_section_singlet(ie), &
                  self%worker(iwork)%medium(1)%total_cross_section_triplet(ie), &
                  sum(self%worker(iwork)%medium(1)%platzman_ionize_orbital(:, ie)), &
                  sum(self%worker(iwork)%medium(1)%platzman_singlet_orbital(:, ie)), &
                  sum(self%worker(iwork)%medium(1)%platzman_triplet_orbital(:, ie)), &
                  self%worker(iwork)%medium(1)%mean_free_path(ie), &
                  self%worker(iwork)%medium(1)%range_ionize(ie), &
                  self%worker(iwork)%medium(1)%range_singlet(ie), &
                  self%worker(iwork)%medium(1)%range_triplet(ie)
          end do

          write(unit,*)
       end do

       close(unit)

       unit = unit_stdout

       do imedia = 1, self%worker(iwork)%nmedia
          call print_yield(unit, imedia)
          call print_gvalue(unit, imedia)
       end do
    end do

    call print_results_degradation(self)
    call print_results_gvalue(self)

  contains

    subroutine print_yield(unit, imedia)
      integer, intent(in) :: unit
      integer, intent(in) :: imedia

      write(unit,'("#")')
      write(unit,'("#",a)') ' Yield:'
      write(unit,'("#",7(a20))') 'Orbital', 'ionize', 'singlet', 'triplet', &
           'Energy i', 'Energy s', 'Energy t'

      do io = 1, self%worker(iwork)%medium(imedia)%number
         write(unit,'("#",15x,i5,6(8x,f12.4))') io, &
              self%worker(iwork)%medium(imedia)%yield_ionize(io), &
              self%worker(iwork)%medium(imedia)%yield_singlet(io), &
              self%worker(iwork)%medium(imedia)%yield_triplet(io), &
              self%worker(iwork)%medium(imedia)%energy_ionize(io), &
              self%worker(iwork)%medium(imedia)%energy_singlet(io), &
              self%worker(iwork)%medium(imedia)%energy_triplet(io)
      end do

      if (self%worker(iwork)%medium(imedia)%number > 1) then
         write(unit,'("#",a20,3(8x,f12.4))') 'sum', &
              sum(self%worker(iwork)%medium(imedia)%yield_ionize), &
              sum(self%worker(iwork)%medium(imedia)%yield_singlet), &
              sum(self%worker(iwork)%medium(imedia)%yield_triplet)
      end if

      write(unit,*)

    end subroutine print_yield

    subroutine print_gvalue(unit, imedia)
      integer, intent(in) :: unit
      integer, intent(in) :: imedia

      write(unit,'("#")')
      write(unit,'("#",a)') ' G-value:'
      write(unit,'("#",7(a20))') 'Orbital', 'ionize', 'singlet', 'triplet', &
           'Energy i', 'Energy s', 'Energy t'

      do io = 1, self%worker(iwork)%medium(imedia)%number
         write(unit,'("#",15x,i5,6(8x,f12.4))') io, &
              self%worker(iwork)%medium(imedia)%gvalue_ionize(io), &
              self%worker(iwork)%medium(imedia)%gvalue_singlet(io), &
              self%worker(iwork)%medium(imedia)%gvalue_triplet(io), &
              self%worker(iwork)%medium(imedia)%energy_ionize(io), &
              self%worker(iwork)%medium(imedia)%energy_singlet(io), &
              self%worker(iwork)%medium(imedia)%energy_triplet(io)
      end do

      if (self%worker(iwork)%medium(imedia)%number > 1) then
         write(unit,'("#",a20,3(8x,f12.4))') 'sum', &
              sum(self%worker(iwork)%medium(imedia)%gvalue_ionize), &
              sum(self%worker(iwork)%medium(imedia)%gvalue_singlet), &
              sum(self%worker(iwork)%medium(imedia)%gvalue_triplet)
      end if

      write(unit,*)

    end subroutine print_gvalue

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

       ngen = size(self%worker(iwork)%medium(1)%degradation_gen, 1)

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
               sum(self%worker(iwork)%medium(1)%degradation_gen(:, ie))
          do igen = 1, ngen
             write(unit,'(1x,e20.12)', advance='no') &
                  self%worker(iwork)%medium(1)%degradation_gen(igen, ie)
          end do
          write(unit,*)
       end do
       close(unit)
    end do

  end subroutine print_results_degradation

  subroutine print_results_gvalue(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    class(work) :: self
    integer :: iwork
    character (len=100) :: file
    integer :: unit
    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'_gvalue.dat'
       call get_unused_unit(unit)
       open(unit, file=trim(file))
       write(unit,'("#",3(1x,a20))') 'ionize', 'singlet', 'triplet'
       write(unit,'(1x,3(1x,e20.12))') &
            sum(self%worker(iwork)%medium(1)%gvalue_ionize), &
            sum(self%worker(iwork)%medium(1)%gvalue_singlet), &
            sum(self%worker(iwork)%medium(1)%gvalue_triplet)
       close(unit)
    end do
  end subroutine print_results_gvalue

end module class_work
