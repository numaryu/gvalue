module class_work
  use class_medium, only: medium
  use class_grid, only: grid
  use class_mixture, only: mixture
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
     type(medium), pointer :: medium(:) => null()
     type(grid) :: egrid
     type(mixture) :: mediamix
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
    integer :: iwork, imedia, igen
    real, allocatable :: stop_power(:,:), degradation(:,:)
    real, allocatable :: number_density(:)
    real, allocatable :: total_cross_section_total(:,:)

    do iwork = 1, self%number
       write(6,'(/1x,a,i0,a,a/)') '=====> WORK ', iwork, ': ',self%worker(iwork)%name

       call init_subwork(self%worker(iwork))

       if (.not.allocated(stop_power)) allocate( &
            stop_power(self%worker(iwork)%nmedia, self%worker(iwork)%egrid%number))
       if (.not.allocated(degradation)) allocate( &
            degradation(self%worker(iwork)%nmedia, self%worker(iwork)%egrid%number))
       if (.not.allocated(number_density)) allocate( &
            number_density(self%worker(iwork)%nmedia))
       if (.not.allocated(total_cross_section_total)) allocate( &
            total_cross_section_total(self%worker(iwork)%nmedia,self%worker(iwork)%egrid%number))

       do imedia = 1, self%worker(iwork)%nmedia
          call self%worker(iwork)%medium(imedia)%init_orbital_vars(self%worker(iwork)%egrid%number, &
               self%worker(iwork)%ngeneration)
          call self%worker(iwork)%medium(imedia)%calculate_stopping_power(self%worker(iwork)%egrid)

          stop_power(imedia,:) = self%worker(iwork)%medium(imedia)%stop_power
          number_density(imedia) = self%worker(iwork)%medium(imedia)%number_density
       end do

       call calculate_mixture(stop_power, self%worker(iwork)%mediamix%stop_power_mixture, &
            op = 'sum')

       do igen = 1, self%worker(iwork)%ngeneration
          do imedia = 1, self%worker(iwork)%nmedia
             call self%worker(iwork)%medium(imedia)%calculate_degradation(self%worker(iwork)%egrid, &
                  self%worker(iwork)%mediamix, igen)
             degradation(imedia,:) = sum(self%worker(iwork)%medium(imedia)%degradation_gen(:,:), dim=1)
          end do
       end do

       call calculate_mixture(degradation, self%worker(iwork)%mediamix%degradation_mixture, &
            ratio = number_density/sum(number_density))

       do imedia = 1, self%worker(iwork)%nmedia
          call self%worker(iwork)%medium(imedia)%calculate_yield(self%worker(iwork)%egrid, &
               self%worker(iwork)%mediamix)
          total_cross_section_total(imedia,:) = self%worker(iwork)%medium(imedia)%total_cross_section_total
       end do

       call calculate_mixture(total_cross_section_total, &
            self%worker(iwork)%mediamix%mean_free_path_mixture, &
            ratio = number_density)
       where (self%worker(iwork)%mediamix%mean_free_path_mixture /= 0.)
          self%worker(iwork)%mediamix%mean_free_path_mixture = &
               1./self%worker(iwork)%mediamix%mean_free_path_mixture
       end where

       if (allocated(stop_power)) deallocate(stop_power)
       if (allocated(degradation)) deallocate(degradation)
       if (allocated(number_density)) deallocate(number_density)
       if (allocated(total_cross_section_total)) deallocate(total_cross_section_total)

       write(6,'(/1x,a,i0,a/)') '=====> WORK ', iwork, ': done'
    end do
  contains

    subroutine init_subwork(worker)
      use mod_file_utils, only: get_unused_unit, unit_stdin, unit_stdout
      type(subwork), intent(in out) :: worker
      integer :: unit
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

      call init_medium(worker)
      call init_grid(worker)
      call init_mixture(worker)

    end subroutine init_subwork

    subroutine init_medium(worker)
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

         worker%medium(imedia) = medium(norbital, name, &
              file_medium(imedia), number_density(imedia))
      end do

    end subroutine init_medium

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

    subroutine init_mixture(worker)
      type(subwork) , intent(in out) :: worker
      worker%mediamix = mixture(worker%egrid%number)
    end subroutine init_mixture

    subroutine calculate_mixture(val, mixed, ratio, op)
      real, intent(in) :: val(:,:)
      real, intent(out) :: mixed(:)
      real, intent(in), optional :: ratio(:)
      character (len=*), optional :: op
      real :: ratio2(size(val,1), size(val,2))
      if (present(ratio)) then
         ratio2 = spread(ratio, 2, size(val,2))
      else if (present(op)) then
         if (op == 'sum') then
            ratio2 = 1
         else if (op == 'mean') then
            ratio2 = 1./size(val,1)
         else
            stop 'invalid option'
         end if
      end if
      mixed(:) = sum(val*ratio2, dim=1)
    end subroutine calculate_mixture

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
          write(unit,'("#",14(1x,a20))') 'Energy', &
               'Stopping Power', 'Stopping Power (mix)', 'Degradation (sum)', &
               'Total Cross Sec. i', 'Total Cross Sec. s', 'Total Cross Sec. t', &
               'Platzman i', 'Platzman s', 'Platzman t', 'Mean Free Path (mix)', &
               'Range i', 'Range s', 'Range t'

          do ie = 1, self%worker(iwork)%egrid%number
             write(unit,'(1x,14(1x,e20.12))') self%worker(iwork)%egrid%val(ie), &
                  self%worker(iwork)%medium(imedia)%stop_power(ie), &
                  self%worker(iwork)%mediamix%stop_power_mixture(ie), &
                  sum(self%worker(iwork)%medium(imedia)%degradation_gen(:, ie)), &
                  self%worker(iwork)%medium(imedia)%total_cross_section_ionize(ie), &
                  self%worker(iwork)%medium(imedia)%total_cross_section_singlet(ie), &
                  self%worker(iwork)%medium(imedia)%total_cross_section_triplet(ie), &
                  sum(self%worker(iwork)%medium(imedia)%platzman_ionize_orbital(:, ie)), &
                  sum(self%worker(iwork)%medium(imedia)%platzman_singlet_orbital(:, ie)), &
                  sum(self%worker(iwork)%medium(imedia)%platzman_triplet_orbital(:, ie)), &
                  self%worker(iwork)%mediamix%mean_free_path_mixture(ie), &
                  self%worker(iwork)%medium(imedia)%range_ionize(ie), &
                  self%worker(iwork)%medium(imedia)%range_singlet(ie), &
                  self%worker(iwork)%medium(imedia)%range_triplet(ie)
          end do

          write(unit,*)
          write(unit,*)
       end do

       close(unit)

       unit = unit_stdout

       do imedia = 1, self%worker(iwork)%nmedia

          write(unit,'("#")')
          write(unit,'("#",1x 2a)') 'Media: ', trim(self%worker(iwork)%medium(imedia)%name)
          write(unit,'("#")')

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

      write(unit,'("#")')

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

      write(unit,'("#")')

    end subroutine print_gvalue

  end subroutine print_results

  subroutine print_results_degradation(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    use class_grid, only: grid
    class(work) :: self
    integer :: iwork, imedia
    integer :: ie, igen, ngen
    integer :: unit
    character (len=100) :: file

    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'_degradation.dat'

       ngen = size(self%worker(iwork)%medium(1)%degradation_gen, 1)

       call get_unused_unit(unit)
       open(unit, file=trim(file))

       do imedia = 1, self%worker(iwork)%nmedia
          write(unit,'("#")')
          write(unit,'("#",1x 2a)') 'Media: ', trim(self%worker(iwork)%medium(imedia)%name)
          write(unit,'("#")')
          write(unit,'("#",2(1x,a20))', advance='no') 'Energy', 'Degradation (sum)'
          do igen = 1, ngen
             write(unit,'(1x, a16, 1x, i3)', advance='no') 'Degradation ', igen
          end do
          write(unit,*)

          do ie = 1, self%worker(iwork)%egrid%number
             write(unit,'(1x,2(1x,e20.12))', advance='no') &
                  self%worker(iwork)%egrid%val(ie), &
                  sum(self%worker(iwork)%medium(imedia)%degradation_gen(:, ie))
             do igen = 1, ngen
                write(unit,'(1x,e20.12)', advance='no') &
                     self%worker(iwork)%medium(imedia)%degradation_gen(igen, ie)
             end do
             write(unit,*)
          end do

          write(unit,*)
          write(unit,*)
       end do

       close(unit)

       file = trim(self%worker(iwork)%name)//'_degradation_mix.dat'

       call get_unused_unit(unit)
       open(unit, file=trim(file))

       write(unit,'("#",2(1x,a20))') 'Energy', 'Degradation (mix)'

       do ie = 1, self%worker(iwork)%egrid%number
          write(unit,'(1x,2(1x,e20.12))') &
               self%worker(iwork)%egrid%val(ie), &
               self%worker(iwork)%mediamix%degradation_mixture(ie)
       end do

       close(unit)
    end do

  end subroutine print_results_degradation

  subroutine print_results_gvalue(self)
    use mod_file_utils, only: get_unused_unit, unit_stdout
    class(work) :: self
    integer :: iwork, imedia
    character (len=100) :: file
    integer :: unit
    do iwork = 1, self%number
       file = trim(self%worker(iwork)%name)//'_gvalue.dat'
       call get_unused_unit(unit)
       open(unit, file=trim(file))

       do imedia = 1, self%worker(iwork)%nmedia
          write(unit,'("#")')
          write(unit,'("#",1x 2a)') 'Media: ', trim(self%worker(iwork)%medium(imedia)%name)
          write(unit,'("#")')
          write(unit,'("#",3(1x,a20))') 'ionize', 'singlet', 'triplet'
          write(unit,'(1x,3(1x,e20.12))') &
               sum(self%worker(iwork)%medium(imedia)%gvalue_ionize), &
               sum(self%worker(iwork)%medium(imedia)%gvalue_singlet), &
               sum(self%worker(iwork)%medium(imedia)%gvalue_triplet)
          write(unit,*)
       end do

       close(unit)
    end do
  end subroutine print_results_gvalue

end module class_work
