include 'mkl_vsl.f90'

program MKL_VSL_TEST

use MKL_VSL_TYPE
use MKL_VSL

REAL(kind = 8) r(100000)
real(kind = 8) s
real(kind = 8) a, b

TYPE (VSL_STREAM_STATE) :: stream

integer*4 errcode
integer*4 i,j
integer brng, method, seed, n

n = 100000
s = 0.0
a = 0.0
b = 1.0
method = VSL_BRNG_METHOD_UNIFORM_STD
brng = VSL_BRNG_MT2203
seed = 777

errcode = vslnewstream( stream, brng, seed )

open(1,file="test.dat")
do i=1, 10
   errcode = vdrnguniform( method, stream, n, r, a, b )
   do j = 2, n
      s = s+r(j)
      write (1,*) r(j), r(j-1)
   end do
end do
close(1)

errcode = vsldeletestream( stream)

print *, "mean of normal distribution", s

end program MKL_VSL_TEST
