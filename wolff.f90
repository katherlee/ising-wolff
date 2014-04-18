include 'mkl_vsl.f90'
program isingWolff
  use MKL_VSL_TYPE
  use MKL_VSL
  implicit none
!  integer, parameter :: nsteps = 10000
  real(kind=8), dimension(10000) :: ms, cor
  real(kind=8) :: mean, var
  integer :: i
!  real*8 :: T
  integer :: brng, seed  
  integer(kind=4) :: errcode
  type (VSL_STREAM_STATE) :: stream
! initialize VSL RNG stream
!$OMP PARALLEL SHARED(brng) PRIVATE(seed, errcode, stream, mean, var, T)
  brng = VSL_BRNG_MT2203
  seed = 16703
!$seed = mod(OMP_get_thread_num()*16703, 997)
  errcode = vslnewstream(stream, brng, seed)
  call vsl_test(errcode)
!$OMP DO
!  do i=1,100
!     T = 2.0d0 + i*0.01d0
!     call wolff(T,100,nsteps,ms,stream)
!     call meanAndVar(ms, nsteps, mean, var)
!     write (*,*) T, mean, var
!  end do
!$OMP END DO
  call wolff(2.d0/log(sqrt(2.d0)+1.d0), 100, 10000, ms, stream)
  call autocor(ms, 10000, cor)

  do i=1,10000
     print *, cor(i)
  end do

!-----------------------!
! Destroying RNG stream !
!-----------------------!
  errcode = vsldeletestream(stream)
  call vsl_test(errcode)
!$OMP END PARALLEL
end program isingWolff

subroutine wolff(T, L, nsteps, ms, stream)
! mkl initialization
  use MKL_VSL_TYPE
  use MKL_VSL

  implicit none

!----------------------!
! Variable declaration !
!----------------------!

! Input parameters
  real(kind=8), intent(in) :: T
  integer, intent(in) :: L
  integer, intent(in) :: nsteps
  real(kind = 8), intent(out), dimension(nsteps) :: ms

! VSL RNG parameters
  type (VSL_STREAM_STATE), intent(inout) :: stream
  integer(kind = 4) ::  errcode
  real(kind = 8), dimension(1) :: rExp
  integer(kind = 4), dimension(2) :: rInt
  

! Total sites
  integer :: N
! Heating steps
  integer, parameter :: heating = 0
! Spin map
  integer(kind=4), dimension(L,L) :: map
  integer(kind=4), dimension(L,L) :: flag
  integer(kind=4), dimension(L) :: iup
  integer(kind=4), dimension(L) :: idown
! Spin stack
  integer, dimension(L*L,2) :: stack
  integer :: top

! Cluster array
  integer, dimension(L*L) :: cluster
  integer :: nCluster

! intermediate variables
  integer :: i,j,k
  integer :: x, y, x2, y2
  integer :: s

!----------------!
! Initialization !
!----------------!

  N = L*L

! initialize VSL RNG stream
!  brng = VSL_BRNG_MT2203
!  seed = 16703
  
!  errcode = vslnewstream(stream, brng, seed)
!  call vsl_test(errcode)
  

! initialize the spin map
  errcode = virnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N, map, 0, 2)
  call vsl_test(errcode)
  map = map*2 - 1

! initialize the neighbor indices
  do i = 1, L
     iup(i) = mod(i, L) + 1
     idown(i) = mod(i+L-2, L) +1
  end do

!-----------!
! Iteration !
!-----------!
  do i = 1, nsteps+heating
     ! Initialize flags
     flag = 0
     ! Initialize stack
     stack = 0
     top = 0
     ! Pick a random site
     errcode = virnguniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 2, rInt, 1, L+1)
     call vsl_test(errcode)
     top = top + 1
     stack(top,:) = rInt
     flag(rInt(1),rInt(2)) = 1
     s = map(rInt(1), rInt(2))
     ! Build a cluster
     clst: do while(top > 0)
!        print *, T, i, top
        x = stack(top,1)
        y = stack(top,2)
        top = top-1
        neighbor: do k = 1,4
           !examine the neighbors
           x2 = x
           y2 = y
           if (k==1)then
              x2 = iup(x)
           elseif (k==2) then
              x2 = idown(x)
           elseif (k==3) then
              y2 = iup(y)
           else !k==4
              y2 = idown(y)
           endif
           ifcan: if( flag(x2,y2) == 0 .AND. map(x2,y2) == s) then
              errcode = vdrngexponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF, &
                   stream, 1, rExp, 0.d0, 1.d0)
              call vsl_test(errcode)
!              print *,T, rExp(1) , 2.d0/T
              ifdo: if(rExp(1) < 2.d0/T) then
                 flag(x2, y2) = 1
                 top = top+1
                 stack(top,:) = (/ x2, y2 /)
              end if ifdo
           end if ifcan
        end do neighbor
     end do clst
     ! flip all spins in the cluster
     forall (x = 1:L, y = 1:L, flag(x,y) == 1) 
        map(x,y) = -map(x,y)
     end forall
     ! save current M state
     j = i-heating
     if (j > 0) then
        ms(j) = 0.d0
        do x=1,L
           do y = 1,L
              ms(j) = ms(j) + dble(map(x,y))
           end do
        end do
        if (ms(j) < 0) then
           ms(j) = -ms(j)
        end if
     end if
  end do
  ms = ms/dble(N)
end subroutine wolff

!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test VSL Error Codes !!
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vsl_test(errcode)
  integer(kind=4), intent(in) :: errcode
  if (errcode /= VSL_ERROR_OK) then
     print *, "error ", errcode, VSL_ERROR_OK
     stop
  end if
end subroutine vsl_test

!!!!!!!!!!!!!!!!!!
!! Mean and Var !!
!!!!!!!!!!!!!!!!!!
subroutine meanAndVar(ms, n, mean, var)
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: ms
  real(kind=8), intent(out) :: mean, var
  integer :: i
  mean = 0.d0
  var = 0.d0
  do i=1,n
     mean = mean + ms(i)
     var = var + ms(i) * ms(i)
  end do
  mean = mean / dble(n)
  var = (var/dble(n) - mean * mean)
end subroutine meanAndVar

!!!!!!!!!!!!!!!!!!!!!
!! Autocorrelation !!
!!!!!!!!!!!!!!!!!!!!!
subroutine autocor(ms, n, cor)
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: ms
  real(kind=8), dimension(n), intent(out) :: cor
  integer :: i, j
  real(kind=8) :: mean, var
  call meanAndVar(ms, n, mean, var)
  do i=1, n
     cor(i) = 0.d0
     do j=0, n-i
        cor(i) = cor(i) + (ms(j)-mean)*(ms(i+j) - mean)
     end do
     cor(i) = cor(i)/dble(n-i+1)/var
  end do
end subroutine autocor
