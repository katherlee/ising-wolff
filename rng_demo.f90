include 'mkl_vsl.f90'

program MKL_VSL_TEST

use MKL_VSL_TYPE
use MKL_VSL

REAL(kind = 8), dimension(10000):: r
real(kind = 8) a, b

TYPE (VSL_STREAM_STATE) :: stream

integer*4 errcode
integer*4 i,j
integer brng, method, seed, n

n = 100000
s = 0.0
a = 0.0
b = 1.0
method = VSL_RNG_METHOD_EXPONENTIAL_ICDF
brng = VSL_BRNG_MT2203
seed = 777

errcode = vslnewstream( stream, brng, seed )

errcode = vdrngexponential( method, stream, 10000, r, a, b )

open(1,file='test.txt')
do i=1,10000
write (1,*) r(i)
end do
close(1)
errcode = vsldeletestream( stream)

end program MKL_VSL_TEST
