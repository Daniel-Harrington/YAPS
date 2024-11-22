
!_____________________________________________________________________________
! modules
!
! pgf95 -Mcuda -ta=nvidia,cc50,time -fast -O2 -Minfo=par -c cufft_module.f95

module precision
! Precision control
integer, parameter, public :: Single = kind(0.0)   ! Single precision
integer, parameter, public :: Double = kind(0.0d0) ! Double precision
!integer, parameter, public :: fp_kind = Double
integer, parameter, public :: fp_kind = Single
end module precision


!
! Define the interface to the NVIDIA CUFFT routines
!

module cufft_interface

integer, public :: CUFFT_FORWARD = -1
integer, public :: CUFFT_INVERSE =  1
integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

! 
!cufftResult
! cufftPlan1d (cufftHandle *plan, int nx,cufftType type,int batch)
!
interface cufftPlan1d
subroutine cufftPlan1d(plan, nx, type, batch) bind(C,name='cufftPlan1d')
use iso_c_binding
integer(c_int):: plan
integer(c_int),value:: nx, batch,type
end subroutine cufftPlan1d
end interface cufftPlan1d

! 
!cufftResult
! cufftPlan1d (cufftHandle *plan, int rank, int *n,
! int *inembed, int istride, int idist,
! int *onembed, int ostride, int odist, cufftType type, int batch)
! 
interface cufftPlanMany
subroutine cufftPlanMany(plan, rank, n, type, batch) bind(C,name='cufftPlanMany')
use iso_c_binding
integer(c_int):: plan, n(*)
integer(c_int),value:: rank, batch, type
end subroutine cufftPlanMany
end interface cufftPlanMany

interface cufftPlan2d
subroutine cufftPlan2d(plan, nx, ny, type) bind(C,name='cufftPlan2d')
use iso_c_binding
integer(c_int):: plan
integer(c_int),value:: nx, ny, type
end subroutine cufftPlan2d
end interface cufftPlan2d

interface cufftPlan3d
subroutine cufftPlan3d(plan, nx, ny, nz, type) bind(C,name='cufftPlan3d')
use iso_c_binding
integer(c_int):: plan
integer(c_int),value:: nx, ny, nz, type
end subroutine cufftPlan3d
end interface cufftPlan3d
 
!
! cufftDestroy(cufftHandle plan)
! 
interface cufftDestroy
subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
use iso_c_binding
integer(c_int),value:: plan
end subroutine cufftDestroy
end interface cufftDestroy
 
!
! cufftExecR2C(cufftHandle plan,
! cufftComplex *idata,
! cufftComplex *odata,
! int direction)
!
interface cufftExecR2C
subroutine cufftExecR2C(plan, idata, odata, direction) &
& bind(C,name='cufftExecR2C')
use iso_c_binding
use precision
integer(c_int),value:: direction
integer(c_int),value:: plan
complex(fp_kind),device:: idata(*),odata(*)
end subroutine cufftExecR2C
end interface cufftExecR2C
 
!
! cufftExecC2C(cufftHandle plan,
! cufftComplex *idata,
! cufftComplex *odata,
! int direction)
! 
interface cufftExecC2C
subroutine cufftExecC2C(plan, idata, odata, direction) &
& bind(C,name='cufftExecC2C')
use iso_c_binding
use precision
integer(c_int),value:: direction
integer(c_int),value:: plan
complex(fp_kind),device:: idata(*),odata(*)
end subroutine cufftExecC2C
end interface cufftExecC2C
 
!
! cufftExecZ2Z(cufftHandle plan,
! cufftDoubleComplex *idata,
! cufftDoubleComplex *odata,
! int direction);
! 
interface cufftExecZ2Z
subroutine cufftExecZ2Z(plan, idata, odata, direction) &
& bind(C,name='cufftExecZ2Z')
use iso_c_binding
use precision
integer(c_int),value:: direction
integer(c_int),value:: plan
complex(fp_kind),device:: idata(*),odata(*)
end subroutine cufftExecZ2Z
end interface cufftExecZ2Z

end module cufft_interface



!__________________________________________________________________
!
!   arithmetic load; sample compute kernel
!__________________________________________________________________
module simpleOps
contains
attributes(global) subroutine calc(a, b, n)
    implicit none
!    real(kind=4) :: a(:), b(:)
    complex(kind=4) :: a(:), b(:)
    integer, value :: n
    integer :: i, m
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x	
    if (i <= n) then 
		do m = 1, 80
		 a(i) = a(i)+ 1e-7* (3.1415+m/1500.) ! ~13 ops
		 b(i) = b(i)+ 2e-7* (1.1415+m*0.124)
 		 a(i) = a(i)+ 1e-7* (b(i)+m/2500.) + .123
		end do ! m
	end if
	call syncthreads
  end subroutine calc
end module simpleOps


!______________________________________________________________________________
!compile with: 
! pgf95 -mp -Mcuda=fastmath,cc35,cc50,cc60,fma,unroll,flushz,lineinfo -ta=nvidia \
! -tp=haswell -fast -O2 -Minfo=par -mcmodel=medium prog_name.f95 -o prog_name.x \
! -L/usr/local/cuda-9.2/lib64 -lcufft -lcupti
!
! dp test requires recompilation of modules with fp_kind = Double
!______________________________________________________________________________

program fft_test
	use precision
	use cufft_interface
	use simpleOps
	implicit none
	integer, parameter :: thr_w = 256
	complex(fp_kind),pinned,allocatable::   a(:),  b(:)
	complex(fp_kind),device,allocatable:: a_d(:),b_d(:)
	integer:: n0, nq, it, i, j
	integer,value :: n
	integer:: plan, maxit, istat, pinnedFlag
	integer:: omp_get_thread_num, omp_get_max_threads
	double precision::  t0,t1, omp_get_wtime, omp_get_wtick
	real :: tim1, tim2, t_p, rate, rate_x, gflops
 	integer:: cudaDeviceSynchronize, device

  do device = 0,1	
!	print*,' which device? [0,1]'
!	read(*,*) device
	call cudaSetDevice(device)

	print*,omp_get_thread_num(), '/', omp_get_max_threads()
	t0 = omp_get_wtime()
!	n0 = 5**5 * 3**5 *2**7   
!	n0 = 5**3 * 3**5 *2**12 
!	n0 = 5**2 * 3**7 *2**11  
!	n0 = 5**0 * 3**8 *2**14  
!	n0 = 5**0 * 3**8 *2**14  
!	n0 = 5**0 * 3**1 *2**25  
	n0 = 5**0 * 3**0 *2**27 
! allocate arrays on the host
	 allocate (a(n0),b(n0),STAT=istat, PINNED=pinnedFlag)
! allocate arrays on the device
	 allocate (a_d(n0),b_d(n0))
!initialize arrays on host
	 a = 1e-6; b = 1.
!copy arrays to device
	 a_d = a
! time
	tim1 = omp_get_wtime() - t0
	t0 = omp_get_wtime() 
	print*,' alloc,Xfer to dev in',tim1,' +-',real(omp_get_wtick())
	print*,' alloc+copy MB/s :',n0*4e-6/tim1
! Print max size of array
	print*,' max size = n0 = ', n0/1024/1024,' M'

!-----loop over sizes of transform---------------------------------	 
  print *,' N, trans[GB/s],  tim1[ms],  Gpt/s,    ns/pt,   GFLOPS'
  do nq = 0,20
	n= n0/2**nq 	 
 
! Initialize the plan
	call cufftPlan1D(plan,n,CUFFT_C2C,1)
	a_d(1:10) = a(1:10); 		b_d(1:10) = b(1:10)
	t0 = omp_get_wtime() 

! Compute kernel 	
	!call calc<<<(1+n/thr_w),thr_w>>>(a_d,b_d,n)
! FFT forward transform
	call cufftExecC2C(plan,a_d,b_d,CUFFT_FORWARD)
	istat = cudaDeviceSynchronize()	
	!call calc<<<(1+n/thr_w),thr_w>>>(a_d,b_d,n)
! FFT backward transform
	call cufftExecC2C(plan,b_d,a_d,CUFFT_INVERSE) 
	istat = cudaDeviceSynchronize()
! timing of one calc 
	t1 = omp_get_wtime()
	tim1 = (t1-t0)/2
 
! Copy results back to host
    b(1:n) = b_d(1:n)
	tim2 = omp_get_wtime() - t1  	! time of xfer b_d -> cpu
! Print initial array & timing (104 sp OPS assumed)
	 t_p = tim1*1e9/(n)  			! time[ns]/pt 
	 rate 	= 2*n/tim1/1e9    		! calc est. Gpt/s (fft) 
 	 rate_x = 2*n*4e-9/tim2  		! GB/s (array b to cpu)
	 gflops = 1e-9*(5*n*alog(n*1.0)/alog(2.))/tim1
	 print ('(i9, f7.2, f10.3, f10.3, f9.3, f10.2)'), &
			n,rate_x,tim1*1e3,rate,t_p,gflops
! Destroy the plan
 	 call cufftDestroy(plan)
  end do ! nq
!----------------------------------------------------------------
  call system('nvidia-smi')
  print*,' istat, pinnedflag', istat, pinnedflag
!release memory on the host
  deallocate (a, b)
!release memory on the device
  deallocate (a_d, b_d) 
 end do ! device



end program fft_test







