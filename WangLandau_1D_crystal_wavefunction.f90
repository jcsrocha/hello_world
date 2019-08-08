!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Definitions of kind value to integer variables
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE integer_kind
  IMPLICIT NONE

  INTEGER, PARAMETER :: K15=selected_int_kind(15)
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
END MODULE integer_kind
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the pseudo-random numbers generator
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE raset1
  USE integer_kind
  IMPLICIT NONE

  REAL(8), SAVE :: u(97), cc, cd, cm
  INTEGER(K4B), SAVE :: i97, j97
END MODULE raset1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the Wang-Landau Algorithm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE WL_variables
  USE integer_kind
  IMPLICIT NONE
  INTEGER(K4B) :: iseed1, iseed2
  INTEGER(K15) :: iE_current, jE_trial, rangeE_min, rangeE_max, &
                & accept_counter, mcs_bin, mcs, nEbins,Np
  REAL(8) :: lnf_min, p, pi, mcs_inv, accept_ratio, epsilon, epsilon_inv, &
           & lnf, delta_E, en_min, en_max, E_current, E_trial, k_current, &
           & k_trial, E0, t0, a , Lx, sigma
  INTEGER(K15),ALLOCATABLE,DIMENSION(:) :: he
  REAL(8),ALLOCATABLE,DIMENSION(:) :: lng
  complex(8),ALLOCATABLE,DIMENSION(:) :: psi,ppsi
  complex(8) :: imag
  LOGICAL :: flat

END MODULE WL_variables

!##################################################################################################!
!Wang-Landau Method for 1D crystal 1 band
!##################################################################################################!
PROGRAM WangLandau_1D_crystal
  USE integer_kind
  USE WL_variables
  IMPLICIT NONE
  INTEGER :: i
  CHARACTER(30) :: lnge

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Parameters and Variables initialization
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  CALL initial()
  call system('rm kxEner.out')
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Wang-Landau Algorithm
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  he = 0; lng = 1.0D0; lnf = 1.0D0 !Initialization of the Variables
  DO WHILE (lnf .GT. lnf_min) 
    CALL sweep()    ! Perform mcs Monte Carlo sweeps (mcs = mcs_bin*nEbins)
    CALL flathist() !Check if histogram is flat
    lng = lng - MINVAL(lng)
    lnf = 0.5D0*lnf
  END DO

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .!
  !Simulation Output (log of the density of states)
  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .!
  OPEN(UNIT=10, FILE='oWL_1D_crystal.out')
  DO i=RangeE_min,RangeE_max
    WRITE(10,*) DBLE(i)*epsilon,lng(i)
  END DO
  CLOSE(10)

  DEALLOCATE(he,lng)

  STOP
!##################################################################################################!
CONTAINS
!------------------------------------------------------------------------------!
!Initialization (initial configuration, constants, table of energy and space Ex)
!------------------------------------------------------------------------------!
SUBROUTINE initial()
USE integer_kind
USE WL_variables
IMPLICIT NONE
  INTEGER(K4B) :: i, j, count, count_rate, count_max
  INTEGER(k4B),DIMENSION(33) :: iseed_aux
  REAL :: seed_aux

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Random Number Generator inicialization
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  DO i=1,33
    DO j=0,i
      CALL SYSTEM_CLOCK(count, count_rate, count_max)
      iseed_aux(i) = count
    END DO
    iseed_aux(i) = iseed_aux(i)*i + iseed_aux(i)*i
  END DO
  CALL random_seed(PUT=iseed_aux)
  DO i=1,10
    CALL random_number(seed_aux)
    iseed1 = NINT(31328.D0*seed_aux)
    CALL random_number(seed_aux)
    iseed2 = NINT(30081.D0*seed_aux)
  END DO
  PRINT*, 'Seeds for the Marsaglia random number generator', iseed1,iseed2
  CALL rmarin(iseed1,iseed2)
  PRINT*, 'The random number generator was initalized.'
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Reading input file
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  OPEN(UNIT=10, FILE='input_WL_1D_crystal.dat')
  READ(10,*) mcs_bin
  READ(10,*) lnf_min
  READ(10,*) p
  READ(10,*) E0, t0, a
  READ(10,*) Lx
  READ(10,*) sigma  
  READ(10,*) epsilon
  CLOSE(10)
  PRINT*, 'Monte Carlo sweep per bin = ', mcs_bin
  PRINT*, 'Minimum value for ln(f) = ', lnf_min
  PRINT*, 'Flatness criteria p = ', p
  PRINT*, 'Energy E = <psi|H|psi>'
  PRINT*, 'psi(x)=A*exp(-sigma*(x-x0)**2+i*kx*(x-x0))'
  PRINT*, 'H=Sum_n  (E0+U(n))| n >< n|-t0 | n +1>< n|-t0 | n -1>< n|'  
  PRINT*, 'E0=', E0
  PRINT*, 'sigma=', sigma
  PRINT*, 't0 = ', t0, 'and a = ', a
  PRINT*, 'Lx = ', Lx, 'and Np = ', Np
  PRINT*, 'Energy discretization with bin size epsilon = ', epsilon

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Defining some variables and allocating memory:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    
  Np=int(Lx/a)
  ALLOCATE(psi(0:Np+1))
  ALLOCATE(ppsi(0:Np+1))

  psi=dcmplx(0.0D0,0.0D0)
  ppsi=psi
  imag=cmplx(0.0,1.0D0)
  epsilon_inv = 1.0D0/epsilon
  pi = 4.0D0*DATAN(1.0D0)  !Defining pi
  
  call psi1D(sigma,0.0D0,0.D0)
  call ppsi1D()
  en_min=DREAL(SUM(CONJG(psi)*ppsi))
  call psi1D(sigma,pi,100.D0)
  call ppsi1D()
  en_max=DREAL(SUM(CONJG(psi)*ppsi))
  
  
  PRINT*, 'Minimum energy = ', en_min
  PRINT*, 'Maximum energy = ', en_max
  rangeE_min = NINT(en_min*epsilon_inv)
  rangeE_max = NINT(en_max*epsilon_inv)
  nEbins = 1 + ABS(rangeE_min) + ABS(rangeE_max)
  mcs = mcs_bin*nEbins       !Monte Carlos sweeps
  PRINT*, 'Monte Carlo sweeps =', mcs
  mcs_inv = 1.0D0/DBLE(mcs)  !inverse of monte carlo sweep
  ALLOCATE(he(rangeE_min:rangeE_max),lng(rangeE_min:rangeE_max))
  PRINT*,'The memory was allocated and the variables was initilizated.'
  
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Finishing the initialization subroutine:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  k_current = pi*ranmar()  !Defining the initial configuration
  CALL energy(E_current, iE_current, k_current)
  PRINT*,'The initialization was done.'

  RETURN
END SUBROUTINE initial

!------------------------------------------------------------------------------!
!Performing mcs Monte Carlo sweeps:
!------------------------------------------------------------------------------!
SUBROUTINE sweep()
  USE integer_kind
  USE WL_variables
  IMPLICIT NONE
  INTEGER(K15) :: k

  accept_counter = 0
  DO k=1,mcs
    k_trial = pi*ranmar()
    CALL energy(E_trial, jE_trial, k_trial)
    CALL wangLandau_step()
  END DO
  accept_ratio = DBLE(accept_counter)*mcs_inv

  RETURN
END SUBROUTINE sweep

!------------------------------------------------------------------------------!
!Energy
!------------------------------------------------------------------------------!
SUBROUTINE energy(E_aux, iE_aux, k_aux)
  USE WL_variables
  IMPLICIT NONE
  INTEGER(K15), INTENT(OUT) :: iE_aux
  REAL(8), INTENT(OUT) :: E_aux
  REAL(8), INTENT(IN) :: k_aux

  call psi1D(sigma,k_aux,100.D0*ranmar())
  call ppsi1D()

  E_aux=DREAL(SUM(CONJG(psi)*ppsi))
  iE_aux = NINT(E_aux*epsilon_inv)

  RETURN
END SUBROUTINE energy

!------------------------------------------------------------------------------!
!Wave function
!------------------------------------------------------------------------------!
SUBROUTINE psi1D(sigma1,kx,x0)
  USE WL_variables
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: kx,x0,sigma1
  REAL(8) :: xi,xi02,norm
  
  DO i=1,Np
   xi = DBLE(i)*a-x0
   xi02 = xi*xi
   psi(i) = cmplx(dcos(kx*xi),dsin(kx*xi))*dexp(-sigma1*xi02)
  END DO
  norm=DREAL(sum(CONJG(psi)*psi))
  psi=psi/sqrt(norm)

  RETURN
END SUBROUTINE psi1D

!------------------------------------------------------------------------------!
!Wave function bracked  
!------------------------------------------------------------------------------!
SUBROUTINE ppsi1D()
  USE WL_variables
  IMPLICIT NONE

  ppsi=cmplx(0.0D0,0.0D0)
  DO i=1,Np
   ppsi(i) = 2.0D0*(E0 + 0.01D0*DBLE(i))*psi(i) - t0*psi(i+1) - t0*psi(i-1)
  END DO
  
  RETURN
END SUBROUTINE ppsi1D

!------------------------------------------------------------------------------!
!Wang-Landau Algorithm
!------------------------------------------------------------------------------!
SUBROUTINE wangLandau_step()
  USE WL_variables
  IMPLICIT NONE
  REAL(8) :: aux1, aux2

  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  !Wang Landau sweep
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  IF (lng(jE_trial) < lng(iE_current)) THEN
    iE_current = jE_trial
    E_current = E_trial
    accept_counter = accept_counter + 1
  ELSE
    aux1 = ranmar()
    aux2 = DEXP(lng(iE_current) - lng(jE_trial))
    IF (aux1 < aux2) THEN
      iE_current = jE_trial
      E_current = E_trial
      accept_counter = accept_counter + 1
    END IF
  END IF
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  !Updating the histogram and the density of states:
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
   lng(iE_current) = lng(iE_current) + lnf
   he(iE_current) = he(iE_current) + 1

  RETURN
END SUBROUTINE wangLandau_step

!------------------------------------------------------------------------------!
!Verifica se o histograma esta flat
!------------------------------------------------------------------------------!
SUBROUTINE flathist()
  USE integer_kind
  USE WL_variables
  IMPLICIT NONE
  INTEGER(K15) :: i, kk
  REAL(8) :: lngE_med(rangeE_min:rangeE_max),media,mini,aux,p2
  LOGICAL :: flat_core

    !--------------------------------------------------!
    !check if the walker myid is flat
    !--------------------------------------------------!
    media=0.0D0; mini = 1.0D9; kk = 0
    flat = .FALSE.
    DO i = rangeE_min, rangeE_max
      kk = kk + 1
      aux = DBLE(he(i))
      media = media + aux
      IF (aux .LT. mini) mini = aux
    END DO
    media=media/DBLE(kk)
    p2=mini/media
    IF (p2 > p) flat = .TRUE.
    PRINT*,'flatness',p2, flat, 'lnf =', lnf

  RETURN
END SUBROUTINE flathist

!------------------------------------------------------------------------------!
!Pseudorandom numbers generator
!------------------------------------------------------------------------------!
SUBROUTINE rmarin(ij, kl)
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The onL changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
USE integer_kind
USE raset1
IMPLICIT NONE
INTEGER(K4B), INTENT(IN) :: ij, kl
INTEGER(K4B) :: i, j, k, l, ii, jj, m
REAL(8) :: s, t

IF( ij < 0  .OR.  ij > 31328  .OR. kl < 0  .OR.  kl > 30081 ) THEN
  WRITE(*, '(A)') ' The first random number seed must have a value ',  &
               'between 0 AND 31328'
  WRITE(*, '(A)') ' The second seed must have a value between 0 and 30081'
  STOP
END IF

i = MOD(ij/177, 177) + 2
j = MOD(ij, 177) + 2
k = MOD(kl/169, 178) + 1
l = MOD(kl, 169)
DO ii = 1, 97
  s = 0.0D0
  t = 0.5D0
  DO jj = 1, 24
    m = MOD(MOD(i*j, 179)*k, 179)
    i = j
    j = k
    k = m
    l = MOD(53*l + 1, 169)
    IF (MOD(l*m, 64) >= 32) THEN
      s = s + t
    END IF
    t = 0.5D0*t
  END DO
  u(ii) = s
END DO
cc = 362436.0D0/16777216.0D0
cd = 7654321.0D0/16777216.0D0
cm = 16777213.0D0/16777216.0D0
i97 = 97
j97 = 33

RETURN
END SUBROUTINE rmarin

!------------------------------------------------------------------------------!
! This is the random number generator proposed by George Marsaglia
! in Florida State University Report: FSU-SCRI-87-50
!------------------------------------------------------------------------------!
FUNCTION ranmar() RESULT(fn_val)
USE raset1
IMPLICIT NONE
REAL(8) :: fn_val
! Local variable
REAL(8):: uni

  uni = u(i97) - u(j97)
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
  u(i97) = uni
  i97 = i97 - 1
  IF(i97 == 0) i97 = 97
  j97 = j97 - 1
  IF(j97 == 0) j97 = 97
  cc = cc - cd
  IF( cc < 0.0D0 ) cc = cc + cm
  uni = uni - cc
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
!  IF( uni == 0.0D0 ) uni = 2.0D-38
  fn_val = uni

  RETURN
END FUNCTION ranmar

!------------------------------------------------------------------------------!
!The same generator as above but for a vector of n random numbers 
!------------------------------------------------------------------------------!
SUBROUTINE vranmar(n, fn_val)
USE integer_kind
USE raset1
IMPLICIT NONE
  INTEGER(K4B), INTENT(IN) :: n
  REAL(8), DIMENSION(1:n),INTENT(OUT) :: fn_val

  REAL(8) :: uni
  INTEGER(K4B) :: i

  DO i = 1, n
    uni = u(i97) - u(j97)
    IF( uni < 0.0D0 ) uni = uni + 1.0D0
    u(i97) = uni
    i97 = i97 - 1
    IF(i97 == 0) i97 = 97
    j97 = j97 - 1
    IF(j97 == 0) j97 = 97
    cc = cc - cd
    IF( cc < 0.0D0 ) cc = cc + cm
    uni = uni - cc
    IF( uni < 0.0D0 ) uni = uni + 1.0D0
!    IF( uni == 0.0D0 ) uni = 2.0D-38
    fn_val(i) = uni
  END DO

  RETURN
END SUBROUTINE vranmar

!------------------------------------------------------------------------------!
!Initial random number generator (numerical recipes)
!------------------------------------------------------------------------------!
DOUBLE PRECISION FUNCTION ran_init(idum)
  USE integer_kind
  IMPLICIT NONE
  INTEGER(K4B), INTENT(INOUT)  :: idum
  INTEGER(K4B), PARAMETER      :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL(8), SAVE                :: am
  INTEGER(K4B), SAVE           :: ix=-1,iy=-1,k

  IF (idum <= 0 .OR. iy < 0) THEN
    am = nearest(1.0,-1.0)/IM
    iy = ior(ieor(888889999,abs(idum)),1)
    ix = ieor(777755555,abs(idum))
    idum = abs(idum) + 1
  END IF
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy - k*IQ) - IR*k
  IF (iy < 0) iy = iy + IM
  ran_init = am*ior(iand(IM,ieor(ix,iy)),1)

  RETURN
END FUNCTION ran_init
!------------------------------------------------------------------------------!

END PROGRAM WangLandau_1D_crystal
!##################################################################################################!
