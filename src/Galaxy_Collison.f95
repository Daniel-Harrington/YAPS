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
    


subroutine check_energy(density_grid,nx,ny,nz,particles,N,smbh_m,E)
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
    implicit none
    
    !########################
    !   Host Initializations
    !#########################
    
    integer,intent(in)::  nx,ny,nz,N
    ! I can operate in place on the density grid since it needs to be
    ! recomputed anyways, this way avoids creating another huge 3d array of grid_dim^3 cells
    ! atleast on the cpu
    real, Dimension( nx,ny,nz) :: density_grid
    
    real, Dimension(9,N)::particles
    real, Dimension(3):: v_i

    ! Energy compute_accelerationsulation 
    real::U,E,KE
    integer::V

     ! Constants

    real :: G = 1 ! Natural Units
    real,parameter:: pi = atan(1.0)*4 
    real:: factor,constants
    real::m,smbh_m
   

    ! Iteration
    integer::i
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan

    ! Wave number stuff

    real::del_kx,del_ky,del_kz


    !###################################
    !   Device Initialization
    !###################################
    real, Dimension(:,:,:), allocatable, device :: density_grid_r_d
    complex(fp_kind), Dimension(:,:,:), allocatable,device:: density_grid_c_d
    call cudaSetDevice(0)



    !#######################################
    !   Forward FFT
    !#######################################
    
    


    ! Allocate input and output arrays on the device memory (gpu)

    allocate(density_grid_r_d(nx,ny,nz))

    density_grid_r_d = density_grid ! transfer from host to device
    

    ! https://docs.nvidia.com/cuda/cufft/index.html#multidimensional-transforms
    allocate(density_grid_c_d(nx,ny,nz/2 +1))

    ! 3D R2C Fourier Transform plan setup
    call cufftPlan3d(plan, nx,ny,nz,CUFFT_R2C)

    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)

    !######################################################
    !Compute Gravitational Potential in Fourier Space
    !#################################################

 
    ! get Volume of cube
    V = nx*ny*nz
    
    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition
    del_kx =  2*pi/nx
    del_ky = 2*pi/ny
    del_kz = 2*pi/nz

    ! Reset U,KE

    U=0.0
    KE = 0.0
    do k_x=1,nx
        do k_y=1,ny
            do k_z=1,nx/2 + 1
                

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
    
            end do
        end do
    end do

    !Destroy Plan
    call cufftDestroy(plan)


    
    ! Release memory
    deallocate (density_grid_r_d,density_grid_c_d)

    U = (2*pi*G/V)*U

    m = 1/N


    ! Add KE for Supermassive seperately 
    ! assuming it is particle 1 (index 0)

    ! get velocities
    v_i = particles(4:6,0)
    KE = KE + smbh_m* 0.5* sum(v_i**2)


    ! get KE 
    do i=1,N
            ! get velocities

        v_i = particles(4:6,0)

        KE = KE + m* 0.5* sum(v_i**2)
    end do
    
    ! combine energies
    E = U + KE

    
end subroutine check_energy

module device_ops
    contains
    attributes(global) subroutine compute_accelerations(gravity_grid_c_d,density_grid_c_d,nx,ny,nz)
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

        ! Cuda Variables for plan process
        !cufftHandle plan identifier and error
        integer::status,plan

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

                K = (abs(k1)**2 + abs(k2)**2 + abs(k3)**2)**(-1)



                !compute once
                p_term = density_grid_c_d(k_x, k_y, k_z) * K * constants

                ! Sets x,y,z accelerations in fourier space on device grid
                gravity_grid_c_d(1, k_x, k_y, k_z) = k1 * p_term
                gravity_grid_c_d(2, k_x, k_y, k_z) = k2 * p_term
                gravity_grid_c_d(3, k_x, k_y, k_z) = k3 * p_term
            endif
        call syncthreads
      end subroutine compute_accelerations
end module device_ops
    


subroutine fft_step(density_grid,nx,ny,nz,particles,N)
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
    !           particles      : real(9,N)
    !               N          : int
    !
     !##########################################################

    use precision
    use cufft_interface
    use device_ops
    implicit none
    
    !########################
    !   Host Initializations
    !###########################
    
    integer,intent(in)::  nx,ny,nz,N
    ! I can operate in place on the density grid since it needs to be
    ! recomputed anyways, this way avoids creating another huge 3d array of grid_dim^3 cells
    ! atleast on the cpu
    real, Dimension( nx,ny,nz) :: density_grid
    real, Dimension( 3,nx,ny,nz) :: gravity_grid
    real, Dimension(9,N)::particles

    ! Constants

  
    real:: factor


    ! Iteration
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan

    ! Wave number stuff

    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition

    real::del_kx,del_ky,del_kz
    
    
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
    blockDimX = 32
    blockDimY = 32
    blockDimZ = 1

    gridDimX = (nx + blockDimX - 1) / blockDimX
    gridDimY = (ny + blockDimY - 1) / blockDimY
    gridDimZ = ((nz/2 +1) + blockDimZ - 1) / blockDimZ


    call cudaSetDevice(0)
    
    !#######################################
    !   Forward FFT
    !#######################################


    ! Allocate input and output arrays on the device memory (gpu)
    allocate(density_grid_r_d(nx,ny,nz))
    density_grid_r_d = density_grid ! transfer from host to device
    
    ! https://docs.nvidia.com/cuda/cufft/index.html#multidimensional-transforms
    allocate(density_grid_c_d(nx,ny,(nz/2 +1)))


    ! allocate the gravity complex and real grids
    allocate(gravity_grid_c_d(3,nx,ny,(nz/2 +1)))
    allocate(gravity_grid_r_d(3,nx,ny,nz))

    

    ! 3D R2C Fourier Transform plan setup
    call cufftPlan3d(plan, nx,ny,nz,CUFFT_R2C)

    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)

    call compute_accelerations<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(gravity_grid_c_d,density_grid_c_d,nx,ny,nz)


    !#######################################
    !   Inverse FFT
    !#######################################

    !Inverse 3D C2R Fourier Transform execution on the Gravity Cube
    call cufftExecC2R(plan,gravity_grid_c_d,gravity_grid_r_d)

    ! Move from device to host
    gravity_grid = gravity_grid_r_d

    ! print*, "Density Grid"
    ! print*, density_grid
    ! print*, "Gravity Grid"
    ! print*, gravity_grid
    ! Normalize Gravity Cube in real space(divide by N/)
    
    factor = 1/N ! precompute to do multiplication instead of division on array ops

    ! Apply factor ONLY to the acceleration dimensions not the index ones
    gravity_grid(1:3, :, :, :) = gravity_grid(1:3, :, :, :) * factor

    
    !Destroy Plan
    call cufftDestroy(plan)

    !release memory on the device
    deallocate(density_grid_r_d,density_grid_c_d,gravity_grid_r_d, gravity_grid_c_d)

    ! ################################
    ! Update particles accelerations
    ! ###############################

    ! Choose a way of mapping cube of force back to particles
    ! inside cube, since density array has less
    ! in my code with same dimension this was easy
    ! just particles(7:9,:) = accel_array
    ! but that only works on same size arr

    particles(7:9,:) = 0!   
    

    
end subroutine fft_step


subroutine integration_step(density_grid,nx,ny,nz, particles, N, dt)
    use precision
    use cufft_interface
    implicit none

    integer,intent(in)::  nx,ny,nz,N
    real, Dimension( nx,ny,nz) :: density_grid

    integer :: i, j
    real, dimension(9, N) :: particles


    real :: dt

   

    !******************************
    ! 2nd order Leapfrog Integration
    !******************************

    ! keep in mind, but avoiding the copy for such a huge
    ! set of 100M particles
    ! x = particles(1:3,:)
    ! v = particles(4:6,:)
    ! a = particles(7:9,:)

    call fft_step(density_grid,nx,ny,nz, particles, N)

   

    ! kick
    particles(4:6,:) = particles(4:6,:) + particles(7:9,:)*dt*0.5

    !drift
    particles(1:3,:) = particles(1:3,:)+particles(4:6,:)*dt

    !kick
    particles(4:6,:) =particles(4:6,:)+particles(7:9,:)*dt*0.5

    

    
end subroutine integration_step

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
    ! compute_accelerationsulate theta for a logarithmic spiral
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

        ! compute_accelerationsulate theta for a logarithmic spiral
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

        ! compute_accelerationsulate theta for a logarithmic spiral
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

subroutine particle_to_grid(density_full, particles, N, nx, ny, nz, dx, dy, dz)
    implicit none
    integer, intent(in) :: N, nx, ny, nz
    real, intent(in) :: particles(9, N), dx, dy, dz
    real, intent(out) :: density_full(nx, ny, nz)
    real :: density(nx/2, ny/2, nz/2)

    integer :: i, j, k, ix, iy, iz !particle index and grid indices
    real :: x, y, z, m, ux, uy, uz
    real :: x_rel, y_rel, z_rel !relative distance of particle in cell
    real :: wx0, wx1, wy0, wy1, wz0, wz1 !interpolation weights
    real :: x_min, y_min, z_min, x_max, y_max, z_max, delta, x_i, y_j, z_k
    integer :: nx_inner, ny_inner, nz_inner, kstart, kend, jstart, jend, istart, iend
    integer :: id, jd, kd

    density = 0.0 
    x_min = -1.5
    x_max = 1.5
    y_min = -1.5
    y_max = 1.5
    z_min = -1.5
    z_max = 1.5
    delta = (x_max - x_min) / ((nx/2)-1) 


    !read particle position from initial array
    do i = 1, N
        if (i == 1 .or. i ==2) then 
            m = 1.0  !assign mass to each particle
        else
            m = 1.0
        end if 
        x = particles(1, i)
        y = particles(2, i)
        z = particles(3, i)
    
        !ignore particles outside the range [-1.0, 1.0]
        if (x < -1.5 .or. x > 1.5 .or. y < -1.5 .or. y > 1.5 .or. z < -1.5 .or. z > 1.5) cycle

        ix = int(floor((x - x_min) / delta)) + 1
        iy = int(floor((y - y_min) / delta)) + 1
        iz = int(floor((z - z_min) / delta)) + 1
    
        if (ix < 1) ix = 1
        if (ix >= nx/2) ix = nx/2 -1
        if (iy < 1) iy = 1
        if (iy >= ny) iy = ny/2 -1
        if (iz < 1) iz = 1
        if (iz >= nz) iz = nz/2 -1

        x_i = x_min + (ix - 1) * delta
        y_j = x_min + (iy - 1) * delta
        z_k = x_min + (iz - 1) * delta
    
        !calcualte relative distance of particle in the cell
        x_rel = (x - x_i) / delta
        y_rel = (y - y_j) / delta 
        z_rel = (z - z_k) / delta
    
        !calculate weights
        wx0 = 1.0 - x_rel 
        wx1 = x_rel 
        wy0 = 1.0 - y_rel 
        wy1 = y_rel 
        wz0 = 1.0 - z_rel 
        wz1 = z_rel 

        !print *, "Particle position:", x, y, z 
        !print *, "Grid index", ix, iy, iz
        !print *, 'Position', x_i, y_j, z_k
        !print *, "Relative position:", x_rel, y_rel, z_rel 
        !print *, "weights", wx0, wx1, wy0, wy1, wz0, wz1
        !print *, "__________________________"
    
    
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

    nx_inner = nx / 2
    ny_inner = ny / 2
    nz_inner = nz / 2
    
    ! Calculate the start and end indices for the loops
    istart = nx / 4 + 1
    iend = istart + nx_inner - 1
    jstart = ny / 4 + 1
    jend = jstart + ny_inner - 1
    kstart = nz / 4 + 1
    kend = kstart + nz_inner - 1
    
    ! Initialize density_full to zero (if not already done)
    density_full = 0.0
    
    ! Loop over the indices to copy density into density_full
    do k = kstart, kend
        kd = k - kstart + 1  ! Corresponding index in density
        do j = jstart, jend
            jd = j - jstart + 1
            do i = istart, iend
                id = i - istart + 1
                density_full(i, j, k) = density(id, jd, kd)
            end do
        end do
    end do

end subroutine particle_to_grid
 

subroutine grid_to_particle(acceleration_grid,particles, N, nx, ny, nz, dx, dy, dz, wx0, wy0, wz0, wx1, wy1, wz1 )
    !
    !Returns accelerations of all particles
    !

    implicit none 
    integer, intent(in) :: N, nx, ny, nz
    real :: particles(9, N), dx, dy, dz

    integer :: i, j, k, ix, iy, iz !particle index and grid indices
    real :: x, y, z, m
    real :: x_rel, y_rel, z_rel !relative distance of particle in cell
    real :: wx0, wx1, wy0, wy1, wz0, wz1 !interpolation weights
    real :: x_min, y_min, z_min
    real,dimension(3,nx,ny,nz):: acceleration_grid
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

subroutine track_a_particle(particle_id, timestep,filename)

    integer :: particle_id 
    real, dimension(9,:), intent(in) :: particles 
    integer, intent(in) :: timestep
    real :: x,y,Z
    real :: x,y,z
    
    !Extract particles's position and velocity 
    x = particles(1,particle_id)
    y = particles(2,particle_id)
    Z = particles(3,particle_id)
    vx = particles(4,particle_id)
    vy = particles(5,particle_id)
    vz = particles(6,particle_id)

    !open the file 
    open(20, file="track_particle.csv", status="unknown", action="write")

    !Write to file
    if (timestep==0) then 
        write(20, '(A)') "timestep,x,y,z,vx,vy,vz"  ! Write header
    end if 
    write write(20, '(I8, 6(F12.6, ","))') timestep, x, y, z, vx, vy, vz

    !close file 
    close(20)
end subroutine Track_a_particle 

program nbody_sim
    use precision
    use cufft_interface
    implicit none
    integer, parameter::N = 100000000
    integer, parameter:: nx =512 , ny = 512, nz = 256
    integer:: checkpoint,steps,k,i
    real:: smbh_m
    real, dimension(9,N)::particles, particle_arr
    real, parameter::dt = 10e-5 ! Needed to keep Energy change way below 10^-5
    real:: E_0,E,Rm,Vm,t_c,curr_time,Rm_0,anim_time, dx, dy, dz, density(nx, ny, nz)
    real,dimension(3)::p
    logical::animate
    integer:: particle_to_track = 50

    dx = 1.0/(nx-1)
    dy = 1.0/(ny-1)
    dz = 1.0/(nz-1)

    smbh_m = 1 ! just set to whatever it is
    animate = .true.
    checkpoint = 10 !s
    
    call initialize_particles2(particles,N,1)

    ! initial E_0
    print*, 'Got past initialization'
    !call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
    print*, 'Got past particle to grid'

    !call check_energy(density,nx,ny,nz,particles,N,smbh_m,E_0)
    print*, 'Got past check energy - lol no'

    do i=1, 1000
        ! These 2 will go inside a do loop until end condition
        !call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)

    ! These 2 will go inside a do loop until end condition
    call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
    print*, 'Got past particle to grid'
    ! add an if for however many steps 
    call check_energy(density,nx,ny,nz,particles,N,smbh_m,E)
    print*, 'Got past second energy check'

        call integration_step(density, nx, ny, nz, particles, N, dt)
        print*, "Done step: ", i

    
    end do

    ! need a step to write out
 
end program nbody_sim 