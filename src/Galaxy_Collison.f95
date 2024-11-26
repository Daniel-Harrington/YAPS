!#########################################
!Compile with
!    pgf95 -mp -Mcuda=fastmath,cc35,cc50,cc60,fma,unroll,flushz,lineinfo -ta=nvidia -tp=haswell -fast -O2 -Minfo=all -mcmodel=medium Galaxy_Collison.f95 -L/usr/local/cuda-9.2/lib64 -lcufft -lcupti
!
!#########################################
!_____________________________________________________________________________
! modules
!
! Define the interface to the NVIDIA CUFFT routines

!

module precision
    ! Precision control
    integer, parameter, public :: Single = kind(0.0)   ! Single precision
    integer, parameter, public :: Double = kind(0.0d0) ! Double precision
    !integer, parameter, public :: fp_kind = Double
    integer, parameter, public :: fp_kind = Single
end module precision

module device_ops
contains
attributes(global) subroutine particle_to_grid_cuda_dummy(density_grid_r_d, nx, ny, nz)
    implicit none

    ! Arguments
    integer, value, intent(in) :: nx, ny, nz
    real, dimension(:,:,:), device :: density_grid_r_d

    ! Local variables
    integer :: i, j, k

    ! Compute thread indices (adjusted for Fortran's 1-based indexing)
    i = (blockIdx%x - 1) * blockDim%x + threadIdx%x + 1
    j = (blockIdx%y - 1) * blockDim%y + threadIdx%y + 1
    k = (blockIdx%z - 1) * blockDim%z + threadIdx%z + 1

    ! Bounds check to avoid out-of-range memory access
    if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        density_grid_r_d(i, j, k) = 1.0
    end if
end subroutine particle_to_grid_cuda_dummy

attributes(global) subroutine compute_gravities(gravity_grid_c_d,density_grid_c_d,nx,ny,nz)
    implicit none
    complex, Dimension(:,:,:):: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:):: gravity_grid_c_d
    real :: G = 1 ! Natural Units
    real,parameter:: pi = atan(1.0)*4 
    real:: constants
    complex:: p_term
    integer,value::  nx,ny,nz

    ! Iteration
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag

    ! Wave number stuff

    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition

    real::del_kx,del_ky,del_kz

    k_x = (blockIdx%x-1)*blockDim%x + threadIdx%x
    k_y = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k_z = (blockIdx%z-1)*blockDim%z + threadIdx%z

    
    !######################################################
    !Compute Gravitational Accelerations in Fourier Space
    !#################################################

    ! The transformed density grid in the reduced dimensions
    ! is then mapped back onto itself as the acceleration grid
    ! replacing the density grid it once was
    ! this saves memory since we need to recompute density each
    ! step anyways
    
    
    constants =  -4*pi*G
    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition
    del_kx =  2*pi/nx
    del_ky = 2*pi/ny
    del_kz = 2*pi/nz

    if (k_x <= nx .and. k_y<=ny .and. k_z <= (nz/2 +1)) then
            ! Splitting positive and negative frequencies for x-y equivalents
            if (k_x < nx/2) then
                k1 = del_kx*k_x
            else
                k1 = del_kx*k_x - nx
            endif
            if (k_y < ny/2) then
                k2 = del_ky*k_y
            else
                k2 = del_ky*k_y - ny
            endif
            
            ! z freq always > 0

            k3 = del_kz*k_z


            !k * p(k) *(-4)*pi*G/|K|^2

            K = (abs(k1)**2 + abs(k2)**2 + abs(k3)**2)**(-1)



            !compute once
            p_term = density_grid_c_d(k_x, k_y, k_z) * K * constants

            ! Sets x,y,z accelerations in fourier space on device grid
            gravity_grid_c_d(1, k_x, k_y, k_z) = k1 * p_term
            gravity_grid_c_d(2, k_x, k_y, k_z) = k2 * p_term
            gravity_grid_c_d(3, k_x, k_y, k_z) = k3 * p_term
        endif
    call syncthreads
end subroutine compute_gravities

attributes(global) subroutine normalize3d(arr,nx,ny,nz,factor)
    implicit none
    integer,value::nx,ny,nz
    real, dimension(nx,ny,nz)::arr
    integer::i,j,K
    
    real,value::factor
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z

    if ( k <= nx .and. j<=ny .and. k <= nz) then
        arr(i,j,k) = arr(i,j,k) *factor
    endif
end subroutine
attributes(global) subroutine calculate_U(density_grid_c_d,nx,ny,nz,U)
    implicit none
    complex, Dimension(:,:,:):: density_grid_c_d
    real,parameter:: pi = atan(1.0)*4 
    integer,value::  nx,ny,nz

    ! Iteration
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag
    real,value:: U


    ! Wave number stuff

    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition

    real::del_kx,del_ky,del_kz

    k_x = (blockIdx%x-1)*blockDim%x + threadIdx%x
    k_y = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k_z = (blockIdx%z-1)*blockDim%z + threadIdx%z

    
    !######################################################
    !Compute Potential in Fourier Space
    !#################################################


    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition
    del_kx =  2*pi/nx
    del_ky = 2*pi/ny
    del_kz = 2*pi/nz

    if (k_x <= nx .and. k_y<=ny .and. k_z < (nz/2 +1)) then
            ! Splitting positive and negative frequencies for x-y equivalents
            if (k_x < nx/2) then
                k1 = del_kx*k_x
            else
                k1 = del_kx*k_x - nx
            endif
            if (k_y < ny/2) then
                k2 = del_ky*k_y
            else
                k2 = del_ky*k_y - ny
            endif
        
            ! z freq always > 0

            k3 = del_kz*k_z


            !k * p(k) *(-4)*pi*G/|K|^2

            K = 1/(abs(k1)**2 + abs(k2)**2 + abs(k3)**2)

            if (k_x > 0) then       

                ! may have errors, debugging
                p_mag = density_grid_c_d(k_x,k_y,k_z)
                U = U + p_mag*K
            end if

        endif
    call syncthreads
end subroutine calculate_U
attributes(global) subroutine integration_step(particles_d, N,dt)
    implicit none
    integer :: i
    integer, value::N
    real, dimension(:,:) :: particles_d
    real, value:: dt


    i = (blockIdx%x-1)*blockDim%x + threadIdx%x

    !******************************
    ! 2nd order Leapfrog Integration
    !******************************

    ! keep in mind, but avoiding the copy for such a huge
    ! set of 100M particles on device
    ! x = particles_d(1:3,:)
    ! v = particles_d(4:6,:)
    ! a = particles_d(7:9,:)

    ! note for other guys, i think remember <= here since fortran arrays are inclusive of N
    if (i<= N) then
        ! kick
        particles_d(4:6,i) = particles_d(4:6,i) + particles_d(7:9,i)*dt*0.5

        !drift
        particles_d(1:3,i) = particles_d(1:3,i)+particles_d(4:6,i)*dt

        !kick
        particles_d(4:6,i) =particles_d(4:6,i)+particles_d(7:9,i)*dt*0.5

    endif

end subroutine integration_step
attributes(global) subroutine calculate_KE(particles_d, N,m,smbh1_m,smbh2_m,KE)

    real, dimension(:,:) :: particles_d
    integer,value::N
    real,value:: m,smbh1_m,smbh2_m

    real,value:: KE
    real,dimension(3),device::v_i
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    

    if (i==1) then
        
        ! Add KE for Supermassive seperately 
        ! assuming it is particle 1 (index 1)

        ! get velocities
        v_i = particles_d(4:6,1)
        KE = KE + smbh1_m* 0.5* sum(v_i**2)
    else if(i==2)then
        ! Add KE for Supermassive seperately 
        ! assuming it is particle 2 (index 2)
        ! get velocities
        v_i = particles_d(4:6,2)
        KE = KE + smbh2_m* 0.5* sum(v_i**2)
    else if (i<=N) then
        ! get velocities

        v_i = particles_d(4:6,i)

        KE = KE + m* 0.5* sum(v_i**2)
    
    endif
end subroutine calculate_KE
! A in works Subroutine that does particle to gid on GPU
attributes(global) subroutine particle_to_grid_cuda(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    use cudafor
    implicit none
    integer, value :: N, nx, ny, nz, nx2, ny2, nz2
    real(kind(0.0)), value :: dx, dy, dz,smbh1_m, smbh2_m
    real(kind(0.0)),dimension(:,:):: particles_d
    real,dimension(:,:,:)::density_grid_r_d

    ! Thread and block indices
    integer :: idx, ix, iy, iz, thread_id,istat
    real(kind(0.0)) :: x, y, z, m
    real(kind(0.0)) :: x_rel, y_rel, z_rel
    real(kind(0.0)) :: wx0, wx1, wy0, wy1, wz0, wz1
    real(kind(0.0)) :: x_min, y_min, z_min, x_max, y_max, z_max, delta, delta_z, x_i, y_j, z_k
    

    ! predefined
    x_min = -1.5
    x_max = 1.5
    y_min = -1.5
    y_max = 1.5
    z_min = -1.5
    z_max = 1.5
    delta = (x_max - x_min) / ((nx/2)-1)
    delta_z = (z_max - z_min) / ((nz/2)-1)

    nx2 = nx/4 !offsets for zero padding
    ny2 = ny/4
    nz2 = nz/4

    density_grid_r_d = 0.0

    !compute global thread ID
    thread_id = (blockIdx%x -1) * blockDim%x + threadIdx%x
    if (thread_id >= N) return

    !read particle positions 
    x = particles_d(1,thread_id)
    y = particles_d(2,thread_id)
    z = particles_d(3,thread_id)

    ! Assign mass based on particle ID
    if (thread_id==1) then
        m = smbh1_m 
    end if 
    if (thread_id==2) then
        m = smbh2_m 
    else
        m = 1.0 /N 
    end if 

    
    ! determine grid cell indicies 
    ix = int(floor((x-x_min)/delta)) + 1
    iy = int(floor((y-y_min)/delta)) + 1
    iz = int(floor((z-z_min)/delta_z)) + 1

    !clamp indecies within bounds 
    if (ix < 1) ix = 1
    if (ix >= nx/2) ix = nx/2 -1
    if (iy < 1) iy = 1
    if (iy >= ny) iy = ny/2 -1
    if (iz < 1) iz = 1
    if (iz >= nz) iz = nz/2 -1

    x_i = x_min + (ix - 1) * delta
    y_j = y_min + (iy - 1) * delta
    z_k = z_min + (iz - 1) * delta_z

    x_rel = (x-x_i)/delta
    y_rel = (y-y_j)/delta
    z_rel = (z-z_k)/delta_z

    ! Claculate weights
    wx0 = 1.0 - x_rel 
    wx1 = x_rel 
    wy0 = 1.0 - y_rel 
    wy1 = y_rel 
    wz0 = 1.0 - z_rel 
    wz1 = z_rel

    ! Update density feiled (atomic operations to prevent race condition)

    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2, iz+nz2), m * wx0 * wy0 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2, iz+nz2), m * wx1 * wy0 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2+1, iz+nz2), m * wx0 * wy1 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2+1, iz+nx2), m * wx1 * wy1 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2, iz+nz2+1), m * wx0 * wy0 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2, iz+nz2+1), m * wx1 * wy0 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2+1, iz+nz2+1), m * wx0 * wy1 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2+1, iz+nz2+1), m * wx1 * wy1 * wz1)

end subroutine particle_to_grid_cuda

    attributes(global) subroutine grid_to_particle_cuda(acceleration_grid, particles, N, nx, ny, nz, dx, dy, dz,smbh1_m, smbh2_m)
    implicit none

    integer, value :: N,nx,ny,nz
    real(kind(0.0)), value :: dx,dy,dz
    real(kind(0.0)), dimension(:,:),device :: particles
    real,dimension(:,:,:,:),device:: acceleration_grid

    !thread and block indecies 
    integer :: i,ix,iy,iz,ix_shifted,iy_shifted,iz_shifted,thread_id
    real, value :: smbh1_m,smbh2_m
    real(kind(0.0)) :: x,y,z,m
    real(kind(0.0)) :: x_rel,y_rel,z_rel
    real(kind(0.0)) :: wx0, wx1, wy0, wy1, wz0, wz1
    real(kind(0.0)) :: x_min, y_min, z_min, x_max, y_max, z_max, delta, delta_z, x_i, y_j, z_k
    real(kind(0.0)) :: acc_x , acc_y, acc_z

    ! Predefined 
    x_min = -1.5
    x_max = 1.5
    y_min = -1.5
    y_max = 1.5
    z_min = -1.5
    z_max = 1.5
    delta = (x_max - x_min) / real(nx/2 - 1)
    delta_z = (z_max - z_min) / real(nz/2 - 1)

    ! Compte global thread ID
    thread_id = (blockIdx%x - 1) * blockDim%x + threadIdx%x
    if (thread_id > N) return 

    ! Read particle positons 
    x = particles(1, thread_id)
    y = particles(2, thread_id)
    z = particles(3, thread_id)

    ! Assign mass based on particle ID
    if (thread_id == 1) then 
        m = smbh1_m
    end if 
    if  (thread_id == 2) then
        m = smbh2_m
    else 
        m = 1/N 
    end if 

    ! Ignore particles outside the range [-1.5, 1.5]
    acc_x = 0.0
    acc_y = 0.0
    acc_z = 0.0

    ! Determine grid cell indecies 
    ix = int(floor((x - x_min) / delta)) + 1
    iy = int(floor((x - y_min) / delta)) + 1
    iz = int(floor((z - z_min) / delta_z)) + 1

    ! Clamp indices within bounds
    if (ix < 1) ix = 1
    if (ix >= nx / 2) ix = nx / 2 - 1
    if (iy < 1) iy = 1
    if (iy >= ny / 2) iy = ny / 2 - 1
    if (iz < 1) iz = 1
    if (iz >= nz / 2) iz = nz / 2 - 1

    x_i = x_min + (ix - 1) * delta
    y_j = y_min + (iy - 1) * delta
    z_k = z_min + (iz - 1) * delta_z
    
    ! Calculate relative distances
    x_rel = (x - x_i) / delta
    y_rel = (y - y_j) / delta
    z_rel = (z - z_k) / delta_z

    ! Calculate weights 
    wx0 = 1.0 - x_rel
    wx1 = x_rel
    wy0 = 1.0 - y_rel
    wy1 = y_rel
    wz0 = 1.0 - z_rel
    wz1 = z_rel

    !initialize accelerations 
    acc_x = 0.0
    acc_y = 0.0
    acc_z = 0.0

    ix_shifted = ix+nx/4
    iy_shifted = iy+ny/4
    iz_shifted = iz+nz/4

    ! Interpolate acceleration from the grid to the particle position
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted, iz_shifted)/m * wx0 * wy0 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted, iz_shifted)/m * wx1 * wy0 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted + 1, iz_shifted)/m * wx0 * wy1 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted + 1, iz_shifted)/m * wx1 * wy1 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted, iz_shifted + 1)/m * wx0 * wy0 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted, iz_shifted + 1)/m * wx1 * wy0 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted + 1, iz_shifted + 1)/m * wx0 * wy1 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1)/m * wx1 * wy1 * wz1

    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted, iz_shifted)/m * wx0 * wy0 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted, iz_shifted)/m * wx1 * wy0 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted + 1, iz_shifted)/m * wx0 * wy1 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted + 1, iz_shifted)/m * wx1 * wy1 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted, iz_shifted + 1)/m * wx0 * wy0 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted, iz_shifted + 1)/m * wx1 * wy0 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted + 1, iz_shifted + 1)/m * wx0 * wy1 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1)/m * wx1 * wy1 * wz1

    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted, iz_shifted)/m * wx0 * wy0 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted, iz_shifted)/m * wx1 * wy0 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted + 1, iz_shifted)/m * wx0 * wy1 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted + 1, iz_shifted)/m * wx1 * wy1 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted, iz_shifted + 1)/m * wx0 * wy0 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted, iz_shifted + 1)/m * wx1 * wy0 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted + 1, iz_shifted + 1)/m * wx0 * wy1 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1)/m * wx1 * wy1 * wz1

    ! update particle 
    particles(7, thread_id) = acc_x
    particles(8, thread_id) = acc_y
    particles(9, thread_id) = acc_z
end subroutine grid_to_particle_cuda 
end module device_ops
    

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
! cufftResult
! cufftPlan1d (cufftHandle *plan, int rank, int *n,
! int *inembed, int istride, int idist,
! int *onembed, int ostride, int odist, cufftType type, int batch)
! 

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
! cufftReal *idata,
! cufftComplex *odata,
!
interface cufftExecR2C
subroutine cufftExecR2C(plan, idata, odata) &
& bind(C,name='cufftExecR2C')
use iso_c_binding
use precision
integer(c_int),value:: plan
real(fp_kind),device:: idata(*)
complex(fp_kind),device::odata(*)
end subroutine cufftExecR2C
end interface cufftExecR2C

!
! cufftExecC2R(cufftHandle plan,
! cufftComplex *idata,
! cufftReal *odata,
!
interface cufftExecC2R
subroutine cufftExecC2R(plan, idata, odata) &
& bind(C,name='cufftExecR2C')
use iso_c_binding
use precision
integer(c_int),value:: plan
complex(fp_kind),device:: idata(*)
real(fp_kind),device::odata(*)
end subroutine cufftExecC2R
end interface cufftExecC2R

end module cufft_interface
    

module your_mom
    implicit none
    
contains
subroutine check_energy(density_grid_r_d,density_grid_c_d,nx,ny,nz,particles_d,N,m,smbh1_m,smbh2_m,E)
    !###########################################################
    ! Instructions:
    !      This function is NON Desctructive, ie not operating
    !      in-place on the density grid, call it after the recompute of the
    !      particle densities and before the integration step
    !      
    !      nx,nz,ny are grid dimensions
    !       
    !      pass in the supermassive black holes mass
    !      into smbh_m
    !
    !       density_grid : real(nx,ny,nz)
    !            nx,ny,nz      : int,int, N
    !           particles      : real(9,N)
    !           smbh_m         : real
    !               N          : int
    !               E          : real
    !
     !##########################################################
    use precision
    use cufft_interface
    use device_ops
    implicit none
    
    !########################
    !   Host Initializations
    !#########################
    
    integer::  nx,ny,nz,N
    
    integer :: blockDimX, blockDimY, blockDimZ
    integer :: gridDimX, gridDimY, gridDimZ


    ! Energy compute_acceleration
    real::U,E,KE,m,smbh1_m,smbh2_m
    
    integer::V

     ! Constants

    real,parameter:: pi = atan(1.0)*4 

    real :: G = 1 ! Natural Units

    ! Iteration
    integer::i

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan

    !###################################
    !   Device Initialization
    !###################################    
    !################################### 

    real,Dimension(:,:),device::particles_d
    real,Dimension(:,:,:),device::density_grid_r_d
    complex,Dimension(:,:,:),device::density_grid_c_d

    
    !#######################################
    !   Forward FFT
    !#######################################
    !print*,"beginning fft"
   

    
    ! 3D R2C Fourier Transform plan setup
    call cufftPlan3d(plan,nx,ny,nz,CUFFT_R2C)
    !print*,"FINISHED     planning"

    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)
    !print*,"FINISHED     fft"

    !######################################################
    !Compute Gravitational Potential in Fourier Space
    !#####################################  ############
    
    ! Define block dimensions
    blockDimX = 8
    blockDimY = 8
    blockDimZ = 8

    gridDimX = (nx + blockDimX - 1) / blockDimX
    gridDimY = (ny + blockDimY - 1) / blockDimY
    gridDimZ = ((nz/2 +1) + blockDimZ - 1) / blockDimZ


 
    ! get Volume of cube
    V = nx*ny*nz

    ! Reset U,KE

    U=0.0
    KE = 0.0
    !print*, "got to just before potential launch"
    call calculate_U<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(density_grid_c_d,nx,ny,nz,U)
    call cudaDeviceSynchronize()

    !print*, "Calculated potential"
    call calculate_KE<<<256,256>>>(particles_d,N,m,smbh1_m,smbh2_m,KE)
    call cudaDeviceSynchronize()

    !print*, "Calculated KE"

    !Destroy Plan
    call cufftDestroy(plan)

    U = (2*pi*G/V)*U

    
    ! combine energies
    E = U + KE

    
end subroutine check_energy

subroutine fft_step(density_grid_r_d,density_grid_c_d,gravity_grid_r_d,gravity_grid_c_d,nx,ny,nz,N)
    !###########################################################
    ! Instructions:
    !      Pass in a density grid to fill density_grid
    !      that will be modified inplace to become accelerations
    !      in a cubic cell region then mapped onto the particles 
    !      corresponding to that region
    !      
    !      nx,nz,ny are grid dimensions
    !
    !       density_grid : real(nx,ny,nz)
    !            nx,ny,nz      : int,int, N
    !
     !##########################################################

    use precision
    use cufft_interface
    use device_ops 
    implicit none
    
    !########################
    !   Host Initializations
    !
    !   Only small administrative things
    !   for composing kernel calls
    !   and cufft functions
    !###########################
    
    integer,intent(in),value::  nx,ny,nz, N

   
    ! Normalization Factor

    real :: factor
  

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan
    
    !########################
    !   Device Initializations
    !#########################


   
    ! Real and complex density on gpu

    real, Dimension(:,:,:), allocatable, device :: density_grid_r_d
    complex, Dimension(:,:,:), allocatable,device:: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:), allocatable,device:: gravity_grid_c_d
    real, Dimension(:,:,:,:), allocatable, device :: gravity_grid_r_d

    integer :: blockDimX, blockDimY, blockDimZ
    integer :: gridDimX, gridDimY, gridDimZ

    ! Define block dimensions
    blockDimX = 8
    blockDimY = 8
    blockDimZ = 8

    gridDimX = (nx + blockDimX - 1) / blockDimX
    gridDimY = (ny + blockDimY - 1) / blockDimY
    gridDimZ = ((nz/2 +1) + blockDimZ - 1) / blockDimZ


    
    !#######################################
    !   Forward FFT
    !#######################################

   
    

    ! 3D R2C Fourier Transform plan setup
    call cufftPlan3d(plan, nx,ny,nz,CUFFT_R2C)
    

    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)


    call compute_gravities<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(gravity_grid_c_d,density_grid_c_d,nx,ny,nz)
    call cudaDeviceSynchronize()

    !print *, "compute_gravities finished"
    !#######################################
    !   Inverse FFT
    !#######################################

    !Inverse 3D C2R Fourier Transform execution on the Gravity Cube
    call cufftExecC2R(plan,gravity_grid_c_d,gravity_grid_r_d)
    call cudaDeviceSynchronize()

    ! !print*, "Density Grid"
    ! !print*, density_grid
    ! !print*, "Gravity Grid"
    ! !print*, gravity_grid
    ! Normalize Gravity Cube in real space(divide by N/)
    

    ! TODO: Check im not crazy and i should be using N
    factor = 1/N
    ! Apply factor ONLY to the acceleration dimensions not the index ones
    
    
    call normalize3d<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(gravity_grid_r_d,nx,ny,nz,factor)
  
    !print *, "normalized"
    !Destroy Plan
    call cufftDestroy(plan)

    
end subroutine fft_step

    
end module your_mom



subroutine initialize_particles2(particle_arr,N,Ra)
    !
    !   Iniatiates the particle positions to form a single galaxy after 2 have been merged 
    !   The larger SMBH in the center and a smaller one now orbating near 
    ! 
    implicit none
    real :: r, theta, pitch_angle, arm_separation, random_offset, rotation_velocity
    Integer, intent(in) :: N, Ra
    Real,Dimension(9,N),intent(out) ::  particle_arr
    real, parameter :: pi = atan(1.0)*4 
    real, parameter :: angle = pi/4 ! Anlge of galaxy 2 relative to galaxy 1
    real, parameter :: offset = 5
    Real::x,y,z,v_x,v_y,v_z,a_x,a_y,a_z
    Integer::i, particles_in_galaxy
    real :: spiral_factor, cos_angle, sin_angle
       
    ! data keyword cleanly sets all these to 0.0 as we need for 
    ! ever initial velocity and acceleration in the cloud
    data v_x,v_y,v_z,a_x,a_y,a_z /6*0.0/
 
    pitch_angle = 15.0 * pi / 180.0      ! Pitch angle in radians
    spiral_factor = 1.0 / tan(pitch_angle)  ! Controls spiral tightness
    arm_separation = pi / 2               ! Separation between arms (4 arms)

    particles_in_galaxy = (N-2)/2 ! equal size galaxies for now 

    cos_angle = cos(angle)
    sin_angle = sin(angle)

    ! Initialize first SMBH
    particle_arr(1:3, 1) = (/ 0.0, 0.0, 0.0 /)  ! Position
    particle_arr(4:6, 1) = (/ 0.0, 0.0, 0.0 /)  ! Velocity
    particle_arr(7:9, 1) = (/ 0.0, 0.0, 0.0 /)  ! Acceleration

    ! Initialize second SMBH
    call random_number(r)
    r = Ra/4 + r * (Ra - Ra/4)
    call random_number(random_offset)
    random_offset = (random_offset - 0.5) * Ra / 10
    ! compute_gravitiesulate theta for a logarithmic spiral
    theta = spiral_factor * log(r) + mod(i, 4) * arm_separation + random_offset / r
    
    particle_arr(1:3, 2) = (/ (r + random_offset) * cos(theta), (r + random_offset) * sin(theta), (2.0 * random_offset - 1.0) * 0.01 * Ra /)  ! Position
    rotation_velocity = sqrt(1000/(1*r + abs(random_offset)))
    particle_arr(4:6, 2) = (/ rotation_velocity * sin(theta), -rotation_velocity * cos(theta), 0.0 /)  ! Velocity
    particle_arr(7:9, 2) = (/ 0.0, 0.0, 0.0 /)  ! Acceleration

    ! Particles in first galaxy
    do i = 3, N
        ! Set radial distance r within the range [Ra/4, Ra] with random variation
        call random_number(r)
        r = Ra/4 + r * (Ra - Ra/4)

        ! Generate random offset for more natural spread around the arms
        call random_number(random_offset)
        random_offset = (random_offset - 0.5) * Ra / 10

        ! compute_gravitiesulate theta for a logarithmic spiral
        theta = spiral_factor * log(r) + mod(i, 4) * arm_separation + random_offset / r

        ! Set x, y, z for a spiral galaxy
        x = (r + random_offset) * cos(theta)
        y = (r + random_offset) * sin(theta)
        z = (2.0 * random_offset - 1.0) * 0.01 * Ra  ! slight z offset for thickness

        ! Assign rotation velocities, decreasing with distance from center
        rotation_velocity = sqrt(1000/(1*r + abs(random_offset)))  ! Example galaxy-like rotation curve
        v_x = rotation_velocity * sin(theta)
        v_y = -rotation_velocity * cos(theta)
        v_z = 0.0

        particle_arr(:,i) = (/x,y,z,v_x,v_y,v_z,a_x,a_y,a_z/)

    end do 

    ! Open a file with a unique unit number
    
    ! open(unit=10, file='particledata.csv', status="replace", action="write")

    ! ! Write header
    ! write(10, '(A)') "x,y,z,v_x,v_y,v_z,a_x,a_y,a_z"

    ! ! Write data
    ! do i = 1, N
    !     write(10, '(9(F12.6, ","))') particle_arr(:, i)
    ! end do

    ! ! Close the file
    ! close(10)
end subroutine initialize_particles2


subroutine initialize_particles(particle_arr,N,Ra)
    !
    !   Iniatiates the particle positions to form 2 spiral galaxies
    !   of uniform density and with inital velocities suitable for stable orbit 
    ! 
    implicit none
    real :: r, theta, pitch_angle, arm_separation, random_offset, rotation_velocity
    Integer, intent(in) :: N, Ra
    Real,Dimension(9,N),intent(out) ::  particle_arr
    real, parameter :: pi = atan(1.0)*4 
    real, parameter :: angle = pi/4 ! Anlge of galaxy 2 relative to galaxy 1
    real, parameter :: offset = 5
    Real::x,y,z,v_x,v_y,v_z,a_x,a_y,a_z
    real :: x_rot, y_rot
    Integer::i, particles_in_galaxy
    real :: spiral_factor, cos_angle, sin_angle
       
    ! data keyword cleanly sets all these to 0.0 as we need for 
    ! ever initial velocity and acceleration in the cloud
    data v_x,v_y,v_z,a_x,a_y,a_z /6*0.0/
 
    pitch_angle = 15.0 * pi / 180.0      ! Pitch angle in radians
    spiral_factor = 1.0 / tan(pitch_angle)  ! Controls spiral tightness
    arm_separation = pi / 2               ! Separation between arms (4 arms)

    particles_in_galaxy = (N-2)/2 ! equal size galaxies for now 

    cos_angle = cos(angle)
    sin_angle = sin(angle)

    ! Initialize first SMBH
    particle_arr(1:3, 1) = (/ 0.0, 0.0, 0.0 /)  ! Position
    particle_arr(4:6, 1) = (/ 0.0, 0.0, 0.0 /)  ! Velocity
    particle_arr(7:9, 1) = (/ 0.0, 0.0, 0.0 /)  ! Acceleration

    ! Initialize second SMBH
    particle_arr(1:3, 2) = (/ offset * cos_angle, offset * sin_angle, 0.0 /)  ! Position
    particle_arr(4:6, 2) = (/ 0.0, 0.0, 0.0 /)  ! Velocity
    particle_arr(7:9, 2) = (/ 0.0, 0.0, 0.0 /)  ! Acceleration

    ! Particles in first galaxy
    do i = 3, particles_in_galaxy + 2
        ! Set radial distance r within the range [Ra/4, Ra] with random variation
        call random_number(r)
        r = Ra/4 + r * (Ra - Ra/4)

        ! Generate random offset for more natural spread around the arms
        call random_number(random_offset)
        random_offset = (random_offset - 0.5) * Ra / 10

        ! compute_gravitiesulate theta for a logarithmic spiral
        theta = spiral_factor * log(r) + mod(i, 4) * arm_separation + random_offset / r

        ! Set x, y, z for a spiral galaxy
        x = (r + random_offset) * cos(theta)
        y = (r + random_offset) * sin(theta)
        z = (2.0 * random_offset - 1.0) * 0.01 * Ra  ! slight z offset for thickness

        ! Assign rotation velocities, decreasing with distance from center
        rotation_velocity = sqrt(1000/(1*r + abs(random_offset)))  ! Example galaxy-like rotation curve
        v_x = rotation_velocity * sin(theta)
        v_y = -rotation_velocity * cos(theta)
        v_z = 0.0

        particle_arr(:,i) = (/x,y,z,v_x,v_y,v_z,a_x,a_y,a_z/)

    end do 

    ! Generate particles for the second galaxy
    do i = particles_in_galaxy + 3, N
        call random_number(r)
        r = Ra/2 + r * (Ra - Ra/2)

        call random_number(random_offset)
        random_offset = (random_offset - 0.5) * Ra / 10

        theta = spiral_factor * log(r) + mod(i, 4) * arm_separation + random_offset / r

        x = (r + random_offset) * cos(theta)
        y = (r + random_offset) * sin(theta)
        z = (2.0 * random_offset - 1.0) * 0.01 * Ra  

        ! Apply rotation to second galaxy
        x_rot = cos_angle * x - sin_angle * y + offset * cos_angle
        y_rot = sin_angle * x + cos_angle * y + offset * sin_angle

        rotation_velocity = sqrt(1000/(1*r + abs(random_offset)))
        v_x = rotation_velocity * sin(theta)
        v_y = -rotation_velocity * cos(theta)
        v_z = 0.0

        particle_arr(:, i) = (/ x_rot, y_rot, z, v_x, v_y, v_z, 0.0, 0.0, 0.0 /)
    end do

    ! Open a file with a unique unit number
    
    open(unit=10, file='particledata.csv', status="replace", action="write")

    ! Write header
    write(10, '(A)') "x,y,z,v_x,v_y,v_z,a_x,a_y,a_z"

    ! Write data
    do i = 1, N
        write(10, '(9(F12.6, ","))') particle_arr(:, i)
    end do

    ! Close the file
    close(10)
end subroutine initialize_particles

subroutine particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
    !
    ! Returns density grip of all particles
    ! TODO : dx dy dz allready contained in particles likely unecessary inputs
    ! can likely be made a bit more concise w vector ops
    ! Notes:
    ! Neat use of cycle
    !
    implicit none
    integer, intent(in) :: N, nx, ny, nz
    real, intent(in) :: particles(9, N), dx, dy, dz
    real, intent(out) :: density(nx, ny, nz)

    integer :: i, j, k, ix, iy, iz !particle index and grid indices
    real :: x, y, z, m
    real :: x_rel, y_rel, z_rel !relative distance of particle in cell
    real :: wx0, wx1, wy0, wy1, wz0, wz1 !interpolation weights
    real :: x_min, y_min, z_min, x_lower, y_lower, z_lower

    density = 0.0 
    x_min = -1.0
    y_min = -1.0
    z_min = -1.0

    m = 1.0 !assign mass to each particle

    !read particle position from initial array
    do i = 1, N
        x = particles(1, i)
        y = particles(2, i)
        z = particles(3, i)
    
        !ignore particles outside the range [-1.0, 1.0]
        if (x < -1.0 .or. x > 1.0 .or. y < -1.0 .or. y > 1.0 .or. z < -1.0 .or. z > 1.0) cycle
        
        !determine the grid index from particle position
        ix = floor((x - x_min) / (2.0*dx)) + 1
        iy = floor((y - y_min) / (2.0*dy)) + 1
        iz = floor((z - z_min) / (2.0*dz)) + 1
    
        if (ix < 1) ix = 1
        if (ix >= nx) ix = nx-1
        if (iy < 1) iy = 1
        if (iy >= ny) iy = ny-1
        if (iz < 1) iz = 1
        if (iz >= nz) iz = nz-1
    
        !compute_gravitiesulate lower bound of the particle
        x_lower = x_min + (ix-1) * (2*dx)
        y_lower = y_min + (iy - 1) * (dy*2.0)
        z_lower = z_min + (iz - 1) * (dz*2.0)

        !compute relative positions
        x_rel = (x - x_lower) / (2*dx)
        y_rel = (y - y_lower) / (2*dy)
        z_rel = (z - z_lower) / (2*dz)
    
        ! !print *, "Particle position:", x, y, z 
        ! !print *, "Grid index", ix, iy, iz
        ! !print *, "Relative position:", x_rel, y_rel, z_rel 
        ! !print *, "__________________________"
        
        !compute_gravitiesulate weights


        wx0 = 1.0 - x_rel 
        wx1 = x_rel 
        wy0 = 1.0 - y_rel 
        wy1 = y_rel 
        wz0 = 1.0 - z_rel 
        wz1 = z_rel 
        
        !update density field
        density(ix, iy, iz) = density(ix, iy, iz) + m * wx0 * wy0 * wz0 

        density(ix+1, iy, iz) = density(ix+1, iy, iz) + m * wx1 * wy0 * wz0 
        density(ix, iy+1, iz) = density(ix, iy+1, iz) + m * wx0 * wy1 * wz0 
        density(ix+1, iy+1, iz) = density(ix+1, iy+1, iz) + m * wx1 * wy1 * wz0

        density(ix, iy, iz+1) = density(ix, iy, iz+1) + m * wx0 * wy0 * wz1
        density(ix+1, iy, iz+1) = density(ix+1, iy, iz+1) + m * wx1 * wy0 * wz1
        density(ix, iy+1, iz+1) = density(ix, iy+1, iz+1) + m * wx0 * wy1 * wz1 
        density(ix+1, iy+1, iz+1) = density(ix+1, iy+1, iz+1) + m * wx1 * wy1 * wz1
        
    end do


end subroutine particle_to_grid
 

subroutine grid_to_particle(acceleration_grid,particles, N, nx, ny, nz, dx, dy, dz, wx0, wy0, wz0, wx1, wy1, wz1 )
    !
    !Returns accelerations of all particles
    !

    implicit none 
    integer, intent(in),value :: N, nx, ny, nz
    real, intent(in),value:: dx, dy, dz
    real , device:: particles(9, N)

    integer :: i, j, k, ix, iy, iz !particle index and grid indices
    real :: x, y, z, m
    real :: x_rel, y_rel, z_rel !relative distance of particle in cell
    real :: wx0, wx1, wy0, wy1, wz0, wz1 !interpolation weights
    real :: x_min, y_min, z_min
    real,dimension(:,:,:,:),device:: acceleration_grid
    real:: acc_x,acc_y,acc_z
    x_min = -1.0
    y_min = -1.0
    z_min = -1.0

    do i = 1, N
        x = particles(1, i)
        y = particles(2, i)
        z = particles(3, i)
    
        !ignore particles outside the range [-1.0, 1.0]
        if (x < -1.0 .or. x > 1.0 .or. y < -1.0 .or. y > 1.0 .or. z < -1.0 .or. z > 1.0) cycle
    
        !determine the grid index from particle position
        ix = floor((x - x_min) / (2.0*dx)) + 1
        iy = floor((y - y_min) / (2.0*dy)) + 1
        iz = floor((z - z_min) / (2.0*dz)) + 1
    
        if (ix < 1) ix = 1
        if (ix >= nx) ix = nx-1
        if (iy < 1) iy = 1
        if (iy >= ny) iy = ny-1
        if (iz < 1) iz = 1
        if (iz >= nz) iz = nz-1

        !initialize acceleration components
        acc_x = 0.0
        acc_y = 0.0
        acc_z = 0.0

        !interpolate acceleration from the grid to the particle position
        acc_x = acc_x + acceleration_grid(1,ix, iy, iz) * wx0 * wy0 * wz0
        acc_x = acc_x + acceleration_grid(1,ix + 1, iy, iz) * wx1 * wy0 * wz0
        acc_x = acc_x + acceleration_grid(1,ix, iy + 1, iz) * wx0 * wy1 * wz0
        acc_x = acc_x + acceleration_grid(1,ix + 1, iy + 1, iz) * wx1 * wy1 * wz0
        acc_x = acc_x + acceleration_grid(1,ix, iy, iz + 1) * wx0 * wy0 * wz1
        acc_x = acc_x + acceleration_grid(1,ix + 1, iy, iz + 1) * wx1 * wy0 * wz1
        acc_x = acc_x + acceleration_grid(1,ix, iy + 1, iz + 1) * wx0 * wy1 * wz1
        acc_x = acc_x + acceleration_grid(1,ix + 1, iy + 1, iz + 1) * wx1 * wy1 * wz1

        acc_y = acc_y + acceleration_grid(2,ix, iy, iz) * wx0 * wy0 * wz0
        acc_y = acc_y + acceleration_grid(2,ix + 1, iy, iz) * wx1 * wy0 * wz0
        acc_y = acc_y + acceleration_grid(2,ix, iy + 1, iz) * wx0 * wy1 * wz0
        acc_y = acc_y + acceleration_grid(2,ix + 1, iy + 1, iz) * wx1 * wy1 * wz0
        acc_y = acc_y + acceleration_grid(2,ix, iy, iz + 1) * wx0 * wy0 * wz1
        acc_y = acc_y + acceleration_grid(2,ix + 1, iy, iz + 1) * wx1 * wy0 * wz1
        acc_y = acc_y + acceleration_grid(2,ix, iy + 1, iz + 1) * wx0 * wy1 * wz1
        acc_y = acc_y + acceleration_grid(2,ix + 1, iy + 1, iz + 1) * wx1 * wy1 * wz1

        acc_z = acc_z + acceleration_grid(3,ix, iy, iz) * wx0 * wy0 * wz0
        acc_z = acc_z + acceleration_grid(3,ix + 1, iy, iz) * wx1 * wy0 * wz0
        acc_z = acc_z + acceleration_grid(3,ix, iy + 1, iz) * wx0 * wy1 * wz0
        acc_z = acc_z + acceleration_grid(3,ix + 1, iy + 1, iz) * wx1 * wy1 * wz0
        acc_z = acc_z + acceleration_grid(3,ix, iy, iz + 1) * wx0 * wy0 * wz1
        acc_z = acc_z + acceleration_grid(3,ix + 1, iy, iz + 1) * wx1 * wy0 * wz1
        acc_z = acc_z + acceleration_grid(3,ix, iy + 1, iz + 1) * wx0 * wy1 * wz1
        acc_z = acc_z + acceleration_grid(3,ix + 1, iy + 1, iz + 1) * wx1 * wy1 * wz1

        ! Update particle acceleration components
        particles(7, i) = acc_x
        particles(8, i) = acc_y
        particles(9, i) = acc_z
    end do

    end subroutine grid_to_particle

! subroutine track_a_particle(particle_id, timestep,filename)

!     integer :: particle_id 
!     real, dimension(9,:), intent(in) :: particles 
!     integer, intent(in) :: timestep
!     real :: x,y,z
!     real :: vx,vy,vz
    
!     !Extract particles's position and velocity 
!     x = particles(1,particle_id)
!     y = particles(2,particle_id)
!     Z = particles(3,particle_id)
!     vx = particles(4,particle_id)
!     vy = particles(5,particle_id)
!     vz = particles(6,particle_id)

!     !open the file 
!     open(20, file="track_particle.csv", status="unknown", action="write")

!     !Write to file
!     if (timestep==0) then 
!         write(20, '(A)') "timestep,x,y,z,vx,vy,vz"  ! Write header
!     end if 
!     write write(20, '(I8, 6(F12.6, ","))') timestep, x, y, z, vx, vy, vz

!     !close file 
!     close(20)
! end subroutine Track_a_particle 

program nbody_sim
    use precision
    use device_ops
    use cufft_interface
    use your_mom
    implicit none
    integer, parameter::N = 257
    integer, parameter:: nx =4 , ny = 4, nz = 4
    real, Dimension(nx,ny,nz):: density_grid_test
    real, Dimension(3,nx,ny,nz):: gravity_grid_test

    integer:: checkpoint,steps,k,i,ierr
    real:: m,smbh1_m,smbh2_m,E_0,E,dx,dy,dz
    real, dimension(9,N)::particles
    real, parameter::dt = 10e-5 ! Needed to keep Energy change way below 10^-5
    real,dimension(3)::p
    logical::animate
    integer:: particle_to_track = 50
    integer :: blockDimX, blockDimY, blockDimZ
    integer :: gridDimX, gridDimY, gridDimZ


    ! ############################################
    ! DEVICE MEMORY SETUP
    ! #######################################
    ! Particles on GPU
    real, Dimension(:,:),allocatable, device::particles_d
    real,device::dt_d,nx_d,ny_d,nz_d,m_d,smbh1_m_d,smbh2_m_d
    integer,device:: N_d

    ! Real and complex density on gpu
    real, Dimension(:,:,:), allocatable, device :: density_grid_r_d
    complex, Dimension(:,:,:), allocatable,device:: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:), allocatable,device:: gravity_grid_c_d
    real, Dimension(:,:,:,:), allocatable, device :: gravity_grid_r_d

    
    ! Allocate input and output arrays of cufft on the device memory (gpu)
    allocate(density_grid_r_d(nx, ny, nz), stat=ierr)
    if (ierr /= 0) then
        !print *, "Allocation of density_grid_r_d failed with error code:", ierr
        stop
    end if
    
    
    ! https://docs.nvidia.com/cuda/cufft/index.html#multidimensional-transforms
    allocate(density_grid_c_d(nx,ny,(nz/2 +1)))


    ! allocate the gravity complex and real grids
    allocate(gravity_grid_c_d(3,nx,ny,(nz/2 +1)))
    allocate(gravity_grid_r_d(3,nx,ny,nz))

    ! allocation particles on device
    allocate(particles_d(9,N))

    
    !##############################################
    !
    !   HOST(CPU) INITIALIZATIONS
    !
    !##############################################

    
    call initialize_particles2(particles,N,1)

    smbh1_m = 1.0  ! just set to whatever it is
    smbh2_m = smbh1_m/10
    m = 1/N
    E = 0
    dx = 1.0/(nx-1) 
    dy = 1.0/(ny-1)
    dz = 1.0/(nz-1)
    !##############################################
    !
    !   Device(GPU) INITIALIZATIONS
    !
    !##############################################
    
    call cudaSetDevice(0)
    nx_d = nx
    ny_d = ny
    nz_d = nz
    smbh1_m_d = smbh1_m
    smbh2_m_d = smbh2_m
    m_d = m
    N_d = N
  
    dt_d = dt


    ! copy particles from host to device memory
    particles_d = particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! beyond this every major step stays in device memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initial E_0
    print*, 'Got past initialization'



    blockDimX = 16
    blockDimY = 16
    blockDimZ = 16
    gridDimX = ceiling(real(nx) / blockDimX)
    gridDimY = ceiling(real(ny) / blockDimY)
    gridDimZ = ceiling(real(nz) / blockDimZ)
        


    

    print*, 'Got past particle to grid'
    ! call particle_to_grid_cuda<<<256,256>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    call particle_to_grid_cuda<<<(N-1)/256,256>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    call cudaDeviceSynchronize()
    !call check_energy(density,nx,ny,nz,particles,N,smbh_m,E_0)
    !print*, 'Got past check energy - lol no'

    do i=1, 2
        ! These 2 will go inside a do loop until end condition
        ! call particle_to_grid_cuda<<<256,256>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
        call particle_to_grid_cuda<<<(N-1)/256,256>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
        call cudaDeviceSynchronize()

        density_grid_test = density_grid_r_d
        print*, "Density Cube:"
        print*, density_grid_test
        ! add an if for however many steps 
        ! again, like fft stays on gpu but composes with a fft call
        !call check_energy(density_grid_r_d, density_grid_c_d, nx, ny, nz, particles_d, N, m, smbh1_m, smbh2_m, E)

        !print*, 'Got past second energy check'


        ! fills the real gravity grid
        ! still stays on gpu
        ! host/cpu - style function is just composing cuda kernel functions

        call fft_step(density_grid_r_d,density_grid_c_d,gravity_grid_r_d,gravity_grid_c_d, nx,ny,nz,N)
        
        print*, "got past fft_step"
        gravity_grid_test = gravity_grid_r_d
        print*, "All HAIL Gravity Cube:"
        print*, gravity_grid_test
        !! here zac call your grid to particles kernel
        !! heres and example you can change dimensions and stuff
        call grid_to_particle_cuda<<<(N-1)/256,256>>>gravity_grid_r_d,particles_d,N_d,nx, ny, nz,dx, dy, dz,smbh1_m,smbh2_m)
        call cudaDeviceSynchronize()
        
        print*, "Got past grid to particle"
        ! integration step pushes all positions
        ! ill need to revisit thread count block size just going quick
        ! to get structure
        call integration_step<<<(N-1)/256,256>>>(particles_d,N,dt)
        call cudaDeviceSynchronize()

        particles = particles_d
        print*, "Particles"
        do i = 1, N
           print*, particles(:,i) 
        end do
        !print*, "Done step: ", i

        ! need a step to [pass back & write out

    end do
    ! Deallocations


    !deallocate density grids
    deallocate(density_grid_r_d,density_grid_c_d)
    

    ! deallocate gravity grids
    deallocate(gravity_grid_r_d,gravity_grid_c_d)

    ! deallocate particles
    deallocate(particles_d)
 
end program nbody_sim