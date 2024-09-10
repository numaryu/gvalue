program gvalue
  ! use mod_constants, only: bb
  use mod_constants, only: bb => bb_compat
  use mod_orbital, only: init_orbital, finish_orbital
  use mod_orbital, only: medium
  use mod_orbital, only: n => norbital_compat
  use mod_orbital, only: a => energy_ionize_compat, b => energy_singlet_compat, c => energy_triplet_compat
  use mod_orbital, only: d => energy_kinetic_compat
  use mod_orbital, only: e => number_electrons_compat
  use mod_grid, only: init_grid, finish_grid
  use mod_grid, only: egrid
  use mod_grid, only: s => egrid_compat
  use mod_grid, only: t0 => egrid_max_compat
  use mod_grid, only: mesh => mesh_compat
  use mod_vars, only: init_vars, finish_vars
  use mod_vars, only: y => ydeg_total
  !
  !     *** Theory of G-value ***
  !
  !     Ionization and singlet and triplet excitations.
  !     Number of atomic or molecular orbitals n く= 20.
  !
  !     References
  !     1) S. Sato, K. Okazaki, S. Ohno, Bull. Chem. Soc. Jpn, 47, 2174 (1974).
  !     2) K. Okazaki, S. Sato, S. Ohno, Bull. Chem. Soc. Jpn, 48, 1411 (1975).
  !
  dimension st(1000),y1(1000),y2(1000),y3(1000),y4(1000),y5(1000),y6(1000)

  integer ::la(20), lb(20), lc(20)
  integer :: mt

  call init_orbital
  call init_grid
  call init_vars

  !
  !     la(i) ... No. of the lowest electron energy above a(i)
  !     lb(i) ... No. of the lowest electron energy above b(i)
  !     lc(i) ... No. of the lowest electron energy above c(i)
  !
  do i=1, n
     la(i) = mesh(a(i),t0)
     lb(i) = mesh(b(i),t0)
     lc(i) = mesh(c(i),t0)
  end do
  mt = lc(1)
  write(6,*) 'mt = ', mt, minval(c), c(1)
  
  write(6,'(a,5x,a,i5)') "#", "n=",n
  write(6,'(a,5x,a,4(8x,a),11x,a,3(5x,a))') &
       "#", "i", "a(i)","b(i)","c(i)","d(i)","e(i)","la(i)","lb(i)","lc(i)"
  write(6,'(a,i6,4f12.4,e15.4,3i10)') ("#",i,a(i),b(i),c(i),d(i),e(i),la(i),lb(i),lc(i),i=1,n)
  
  !     
  !     Calculation of stopping power
  !     st(i) ... stopping power at T=s(i)
  !
  call medium%init_orbital_vars(egrid%number)
  call medium%calculate_stopping_power()
  do i=1,6
     call medium%calculate_degradation(i)
  end do
  
  diff = 0.
  do i=1,mt
     sst=0.0
     do j=1,n
        sst=sst+sop(j,i)
     end do
     st(i)=sst
     diff = diff + abs(st(i)-medium%stop_power(i))
  end do
  write(6,*) 'st diff',diff
  
  !
  !     Calculation of degradation spectrum
  !
  !     y1(i) ... Degradation spectrum of primary electrons at T=s(i)
  !     y2(i) ... Degradation spectrum of secondary electrons at T=s(i)
  !     y3(i) ... Degradation spectrum of tertiary electrons at T=s(i)
  !     y4(i) ... Degradation spectrum of fourth electrons at T=s(i)
  !     y5(i) ... Degradation spectrum of fifth electrons at T=s(i)
  !     y6(i) ... Degradation spectrum of sixth electrons at T=s(i)
  !     y(i)  ... Total degradation spectrum of electrons at T=s(i)
  !
  !     Upper and lower limits of integral
  !
  do i=1,mt
     y1(i)=1.0/st(i)
  end do

  na=mesh(t0,t0)+1
  t1=(s(na)-a(1))/2.0
  nb=mesh(t1,t0)+1
  t2=(s(nb)-a(1))/2.0
  nc=mesh(t2,t0)+1
  t3=(s(nc)-a(1))/2.0
  nd=mesh(t3,t0)+1
  t4=(s(nd)-a(1))/2.0
  ne=mesh(t4,t0)+1
  t5=(s(ne)-a(1))/2.0
  nf=mesh(t5,t0)+1
  write(6,*) 't0, t1, t2, t3, t4, t5: ', t0, t1, t2, t3, t4, t5
  write(6,*) 'na, nb, nc, nd, ne, nf: ', na, nb, nc, nd, ne, nf
  ! t_i is the upper bound of T2
  !     t5 < t4 < t3 < t2 < t1 < t0
  !     nf > ne > nd > nc > nb > na
  call yi(y1,y2,st,nb,na,lc(1), 2.0,n,t0)
  call yi(y2,y3,st,nc,nb,lc(1), 4.0,n,t0)
  call yi(y3,y4,st,nd,nc,lc(1), 8.0,n,t0)
  call yi(y4,y5,st,ne,nd,lc(1),16.0,n,t0)
  call yi(y5,y6,st,nf,ne,lc(1),32.0,n,t0)

  diff1 = 0.
  diff2 = 0.
  diff3 = 0.
  diff4 = 0.
  diff5 = 0.
  diff6 = 0.
  do i=1, mt
     diff1 = diff1 + abs(y1(i) - medium%degradation_gen_total(1, i))
     diff2 = diff2 + abs(y2(i) - medium%degradation_gen_total(2, i))
     diff3 = diff3 + abs(y3(i) - medium%degradation_gen_total(3, i))
     diff4 = diff4 + abs(y4(i) - medium%degradation_gen_total(4, i))
     diff5 = diff5 + abs(y5(i) - medium%degradation_gen_total(5, i))
     diff6 = diff6 + abs(y6(i) - medium%degradation_gen_total(6, i))
     write(333,*) s(i), y2(i), medium%degradation_gen_total(2, i)
  end do
  write(6,*) 'y1 diff', diff1
  write(6,*) 'y2 diff', diff2
  write(6,*) 'y3 diff', diff3
  write(6,*) 'y4 diff', diff4
  write(6,*) 'y5 diff', diff5
  write(6,*) 'y6 diff', diff6
!  stop
  
  do i=1,mt
     y(i)=y1(i)+y2(i)+y3(i)+y4(i)+y5(i)+y6(i)
  end do

!!$  write(6,'(a,5x,a,10x,a,6(9x,a),10x,a)') &
!!$       & "#","i","s(i)","y1(i)","y2(i)","y3(i)","y4(i)","y5(i)","y6(i)","y(i)"
!!$  write(6,'(1x,i6,8e14.5)') (i,s(i),y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y(i),i=1,mt)

  open(10,file='degradation.dat',status='unknown')
  write(10,'(a,5x,a,10x,a,6(9x,a),10x,a)') &
       & "#","i","s(i)","y1(i)","y2(i)","y3(i)","y4(i)","y5(i)","y6(i)","y(i)"
  write(10,'(1x,i6,8e14.5)') (i,s(i),y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y(i),i=1,mt)
  close(10)

  stop
  !
  !     Calculation of G-value
  !
  do i=1,n
     call yield(la(i),lb(i),lc(i),1,2,t0,i)
  end do

  stop

end program gvalue

function sop(n,k)
  use mod_constants, only: bb => bb_compat
  use mod_orbital, only: a => energy_ionize_compat, b => energy_singlet_compat, c => energy_triplet_compat
  use mod_orbital, only: d => energy_kinetic_compat
  use mod_orbital, only: e => number_electrons_compat
  use mod_grid, only: s => egrid_compat
  !
  !     Function of stopping power
  !       derived from differential cross section
  !

  ux=a(n)
  ezs=b(n)
  ezt=c(n)
  ex=d(n)
  p=s(k)+ux
  r=s(k)+ux-ezs
  u=s(k)+ux-ezt
  eee=e(n)

  if (s(k) >= ux) then

     sop=eee*bb/(p+ex)*(1.0/2.0*alog(p**4/(16.0*ezs**2*r*u)) &
          +2.0-p/(2.0*r)-p/(2.0*u) &
          +2.0*ex/3.0* &
          (2.0/ezs+1.0/u+1.0/r-4.0/p-p/2.0*(1.0/r**2+1.0/u**2)) &
          )

  else if (s(k) >= ezs) then
  
     sop=eee*bb/(p+ex)*(1.0/2.0*alog(s(k)**2*ux**2/(ezs**2*r*u)) &
          +1.0+s(k)/ux-p/(2.0*r)-p/(2.0*u) &
          +ex/3.0* &
          (2.0*p/ux**2-4.0/ux+2.0/r+2.0/u &
          -p/r**2-p/u**2+4.0/ezs-4.0/s(k)) &
          )

  else if (s(k) >= ezt) then
  
     sop=eee*bb/(p+ex)*(1.0/2.0*alog(ux/u) &
          +p/(2.0*ux) &
          -p/(2.0*u) &
          +2.0*ex/3.0* &
          (1.0/u-1.0/ux-p/(2.0*u**2)+p/(2.0*ux**2)) &
          )
  
  else  
     sop = 0.0
  end if

end function sop

subroutine yi(ya,yb,st,na,nb,mt,x0,m,tx)
  use mod_constants, only: bb => bb_compat
  use mod_orbital, only: a => energy_ionize_compat, b => energy_singlet_compat, c => energy_triplet_compat
  use mod_orbital, only: d => energy_kinetic_compat
  use mod_orbital, only: e => number_electrons_compat
  !
  !     Subroutine 2 of calculation for degradation spectrum
  !       Second integral
  !
  common /yt/ ra(20,1000)
  dimension ya(1000),yb(1000),st(1000)

  ! i here is the energy grid index
  ! na is the upper bound of generated electron energy (uppder bound of T_{2})
  do i=1,na-1
     yb(i)=0.0
  end do
     
  if(na.gt.mt) return
     
  ! for each i-shell
  ! na is the upper bound of T2
  ! nb is the upper bound of T1
  do i=1,m
     call ysin(ya,a(i),d(i),i,na,nb,mt,x0,tx)
  end do

  write(6,*) 'na mt', na, mt
  za=0.0
  ! energy grid
  do i=na,mt
!!     za = 0.0
     ! shell
     do j=1,m
        za=za+ra(j,i)*e(j)
     end do
     yb(i)=za/st(i)
  end do

  return
end subroutine yi

subroutine ysin(yc,ux,ex,l,mm,ll,mt,xa,tt0)
  use mod_constants, only: bb => bb_compat
  use mod_orbital, only: a => energy_ionize_compat, b => energy_singlet_compat, c => energy_triplet_compat
  use mod_orbital, only: d => energy_kinetic_compat
  use mod_orbital, only: e => number_electrons_compat
  use mod_grid, only: s => egrid_compat
  !
  !     Subroutine 1 of caluclation for degradation spectrum
  !       First integral
  !
  !     Differential cross section for production of secondary
  !       electron with energy of s(i) in collision between incident
  !       electron with energy of s(j) and molecular electron
  !       in k-th orbital
  !
  !     ux ... Ii
  !     ex ... Ei
  !     mm < ll ... second integration range indices
  !
  !     l ... shell index
  !
  !     mm is the upper bound of T2
  !     ll is the upper bound of T1
  !     mt is the lower bound index (triplet excitation)
  !
  common /yt/ ra(20,1000)
  dimension yc(1000)
  
  !     determines the upper bound index of T2
  do i=mm,mt+1
     ra(l,i-1)=0.0
     tm=(tt0-ux*(xa-1.0))/xa
     k=i
     if (s(i).le.tm) exit
  end do

  write(6,*) 'k, ll, m', k, ll, mt
!!! s(i) = T2
  do i=k,mt
     vx2=0.0
!!! s(j) = T1 
     do j=ll,mt
        if(s(j)-2.0*s(i)-ux > 0.) then
!!!   Arthmetic IF statement (to be deleted in Fortran 2018)
!!!     if s(j)-2.0*s(i)-ux < 0, jump to 4
!!!     if s(j)-2.0*s(i)-ux = 0, jump to 4
!!!     if s(j)-2.0*s(i)-ux > 0, statement 9 is performed
           vx2=vx2+yc(j)*(s(j-1)-s(j))*bb/(s(j)+ux+ex)*( &
                1.0/(s(i)+ux)**2 &
                +4.0/3.0*ex/(s(i)+ux)**3 &
                +4.0/3.0*ex/(s(j)-s(i))**3 &
                +1.0/(s(j)-s(i))**2 &
                )*(s(i-1)-s(i))
        else
           exit
        end if
     end do
     ra(l,i)=vx2
  end do
  
  return
end subroutine ysin

subroutine yield(ni,ns,nt,lx,lw,tx,k)
  use mod_constants, only: bb => bb_compat
  use mod_orbital, only: a => energy_ionize_compat, b => energy_singlet_compat, c => energy_triplet_compat
  use mod_orbital, only: d => energy_kinetic_compat
  use mod_orbital, only: e => number_electrons_compat
  use mod_grid, only: s => egrid_compat
  use mod_vars, only: y => ydeg_total
  !
  !     Subroutine of calculation for G-value
  !
  dimension qi(1000),qs(1000),qt(1000), &
       dxi(1000),dxs(1000),dxt(1000)
  character (len=20) :: yield_file

  ux=a(k)
  eys=b(k)
  eyt=c(k)
  ex=d(k)

  !
  !     Calculation of cross section Q
  !
  !     qi(i) ... Total cross section of ionization at T=s(i)
  !     qs(i) ... Total cross section of singlet excitation at T=s(i)
  !     qt(i) ... Total cross section of triplet excitation at T=s(i)
  !
  !       These were derived from differential cross section
  !
  !     Calculation of T*Y*Q at T
  !
  !     T = s(i) ... electron energy
  !     Y = y(i) ... degradation spectrum
  !     Q = qi(i),qs(i),qt(i) ... total cross section
  !
  !     dxi(i) ... T*Y*Q for ionization at T
  !     dxs(i) ... T*Y*Q for singlet excitation at T
  !     dxt(i) ... T*Y*Q for triplet excitation
  !
  do i=1,ni
     qi(i)=bb/(s(i)+ux+ex)*( &
          1.0/ux-1.0/s(i) &
          +2.0*ex/(3.0*ux**2)-2.0*ex/(3.0*s(i)**2) &
          )
     dxi(i)=s(i)*y(i)*qi(i)

     qs(i)=bb/(s(i)+ux+ex)*( &
          1.0/eys-1.0/ux &
          +1.0/2.0/s(i)-1.0/2.0/(s(i)+ux-eys) &
          +ex/3.0*( &
          2.0/eys**2-2.0/ux**2+1.0/s(i)**2-1.0/(s(i)+ux-eys)**2 &
          ) &
          )
     dxs(i)=s(i)*y(i)*qs(i)

     qt(i)=bb/2.0/(s(i)+ux+ex)*( &
          1.0/s(i)-1.0/(s(i)+ux-eyt) &
          +2.0/3.0*ex*(1.0/s(i)**2-1.0/(s(i)+ux-eyt)**2))
     dxt(i)=s(i)*y(i)*qt(i)
  end do

  if (ni+1 <= ns) then

     do i=ni+1,ns
        qs(i)=bb/(s(i)+ux+ex)*(1.0/eys-1.0/s(i)+1.0/(2.0*ux) &
             -1.0/(2.0*(s(i)+ux-eys))+ex/3.0*(2.0/eys**2 &
             -2.0/s(i)**2+1.0/ux**2 &
             -1.0/(s(i)+ux-eys)**2))
        dxs(i)=s(i)*y(i)*qs(i)
     end do

  end if

  if(ni+1 <= nt) then

     do i=ni+1,nt
        qt(i)=bb/(2.0*(s(i)+ux+ex))*(1.0/ux-1.0/(s(i)+ux-eyt) &
             +2.0/3.0*ex*(1.0/ux**2-1.0/(s(i)+ux-eyt)**2))
        dxt(i)=s(i)*y(i)*qt(i)
     end do

  end if

  !
  !     Calculation of G-value
  !
  !     Integral of T*Y*Q
  !     G = 100 * NS/T0
  !     NS  ... Number of species produced in S-process
  !
  !     fxi ... Number of ionization from k-th orbital
  !     gxi ... G-value of ionization from k-th orbital
  !     fxs ... Number of singlet excitations from k-th orbital
  !     gxs ... G-value of singlet excitations from k-th orbital
  !     fxt ... Number of triplet excitations from k-th orbital
  !     gxt ... G-value of triplet excitations from k-th orbital
  !

  wi=0.0
  do i=lw,ni
     wi=wi+dxi(i)
  end do
  wi=wi+dxi(lx)/2.0
  fxi=1.0/40.0*alog(2.0)*e(k)*wi
  gxi=100.0*fxi/tx

  ws=0.0
  do i=lw,ns
     ws=ws+dxs(i)
  end do
  ws=ws+dxs(lx)/2.0
  fxs=1.0/40.0*alog(2.0)*e(k)*ws
  gxs=100.0*fxs/tx

  wt=0.0
  do i=lw,nt
     wt=wt+dxt(i)
  end do
  wt=wt+dxt(lx)/2.0
  fxt=1.0/40.0*alog(2.0)*e(k)*wt
  gxt=100.0*fxt/tx
  
  write(yield_file,'(a,i0,a)') "yield",k,".dat"

  open(30,file=yield_file,status='unknown')
  write(30,'(a,i0)') "# shell: ",k
  write(30,'(a,6(10x,a))') "#","hi(i)","gi(i)","hs(i)","gs(i)","ht(i)","gt(i)"
  write(30,'(a,6e15.5)') "#",fxi,gxi,fxs,gxs,fxt,gxt
  write(30,'(a,5x,a,10x,a,3(9x,a),3(8x,a))') &
       "#","i","s(i)","qi(i)","qs(i)","qt(i)","dxi(i)","dxs(i)","dxt(i)"
  write(30,'(1x,i6,7e14.5)') (i,s(i),qi(i),qs(i),qt(i),dxi(i),dxs(i),dxt(i),i=1,nt)
  close(30)

  return
end subroutine yield
