!
!     *** Theory of G-value ***
!
!     Ionization and singlet and triplet excitations.
!
!     The original version (S. Sato, Radiation Chemistry 14, 2-19 (1979) [in Japanese])
!     is written in Fortran 77, is rewritten by numaryu into Fortran 2023.
!
!     References
!     1) S. Sato, K. Okazaki, S. Ohno, Bull. Chem. Soc. Jpn, 47, 2174 (1974).
!     2) K. Okazaki, S. Sato, S. Ohno, Bull. Chem. Soc. Jpn, 48, 1411 (1975).
!
program gvalue
  use class_work, only: work
  type(work) :: mywork
  
  mywork = work()
  call mywork%execute()

  stop

end program gvalue
