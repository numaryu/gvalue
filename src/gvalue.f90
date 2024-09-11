!
!     *** Theory of G-value ***
!
!     Ionization and singlet and triplet excitations.
!
!     References
!     1) S. Sato, K. Okazaki, S. Ohno, Bull. Chem. Soc. Jpn, 47, 2174 (1974).
!     2) K. Okazaki, S. Sato, S. Ohno, Bull. Chem. Soc. Jpn, 48, 1411 (1975).
!
program gvalue
  use mod_orbital, only: init_orbital, finish_orbital
  use mod_orbital, only: medium
  use mod_grid, only: init_grid, finish_grid
  use mod_grid, only: egrid

  call init_orbital
  call init_grid

  call medium%init_orbital_vars(egrid%number)
  call medium%calculate_stopping_power()
  call medium%calculate_degradation()
  call medium%calculate_yield()

  call medium%print_results()
  
  stop

end program gvalue
