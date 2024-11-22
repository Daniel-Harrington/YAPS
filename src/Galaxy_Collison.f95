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

    ! Energy Calculation 
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
    do k_x=0,nx
        do k_y=0,ny
            do k_z=0,nx/2 + 1
                

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

subroutine compute_accelerations(density_grid,nx,ny,nz,particles,N)
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

    real :: G = 1 ! Natural Units
    real,parameter:: pi = atan(1.0)*4 
    real:: factor,constants

    complex:: p_term

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

    print*,"Got here"
    call cudaSetDevice(0)
    print*,"Got here"

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

    do k_x=0,nx
        do k_y=0,ny
            do k_z=0,nx/2 + 1
                

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



                !compute once
                p_term = density_grid_c_d(k_x, k_y, k_z) * K * constants

                ! Sets x,y,z accelerations in fourier space on device grid
                gravity_grid_c_d(1, k_x, k_y, k_z) = k1 * p_term
                gravity_grid_c_d(2, k_x, k_y, k_z) = k2 * p_term
                gravity_grid_c_d(3, k_x, k_y, k_z) = k3 * p_term


            end do
        end do
    end do
    

    !#######################################
    !   Inverse FFT
    !#######################################

    !Inverse 3D C2R Fourier Transform execution on the Gravity Cube
    call cufftExecC2R(plan,gravity_grid_c_d,gravity_grid_r_d)

    ! Move from device to host
    gravity_grid = gravity_grid_r_d


    ! Normalize Gravity Cube in real space(divide by N/)
    
    factor = 1/N ! precompute to do multiplication instead of division on array ops

    ! Apply factor ONLY to the acceleration dimensions not the index ones
    gravity_grid(1:3, :, :, :) = gravity_grid(1:3, :, :, :) * factor


    !Destroy Plan
    call cufftDestroy(plan)

    !release memory on the device
    deallocate (density_grid_r_d,density_grid_c_d,gravity_grid_r_d, gravity_grid_c_d)

    ! ################################
    ! Update particles accelerations
    ! ###############################

    ! Choose a way of mapping cube of force back to particles
    ! inside cube, since density array has less
    ! in my code with same dimension this was easy
    ! just particles(7:9,:) = accel_array
    ! but that only works on same size arr

    particles(7:9,:) = 0!   
    

    
end subroutine compute_accelerations


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

    call compute_accelerations(density_grid,nx,ny,nz, particles, N)

   

    ! kick
    particles(4:6,:) = particles(4:6,:) + particles(7:9,:)*dt*0.5

    !drift
    particles(1:3,:) = particles(1:3,:)+particles(4:6,:)*dt

    !kick
    particles(4:6,:) =particles(4:6,:)+particles(7:9,:)*dt*0.5

    

    
end subroutine integration_step

subroutine initiate_particles(particle_arr,N,Ra)
    !
    !   Iniatiates the particle positions to form a spherical cloud
    !   of uniform density and with inital velocities of 0
    ! 
    implicit none
    real :: r, theta, pitch_angle, arm_separation, random_offset, rotation_velocity
    Integer, intent(in) :: N, Ra
    Real,Dimension(9,N),intent(out) ::  particle_arr
    real, parameter :: pi = atan(1.0)*4 
    Real::x,y,z,v_x,v_y,v_z,a_x,a_y,a_z
    Integer::i
    real :: spiral_factor
    
    
    ! data keyword cleanly sets all these to 0.0 as we need for 
    ! ever initial velocity and acceleration in the cloud
    data v_x,v_y,v_z,a_x,a_y,a_z /6*0.0/
 

    pitch_angle = 15.0 * pi / 180.0      ! Pitch angle in radians
    spiral_factor = 1.0 / tan(pitch_angle)  ! Controls spiral tightness
    arm_separation = pi / 2               ! Separation between arms (4 arms)


    do i = 1,N
        ! Set radial distance r within the range [Ra/4, Ra] with random variation
        call random_number(r)
        r = Ra/4 + r * (Ra - Ra/4)

        ! Generate random offset for more natural spread around the arms
        call random_number(random_offset)
        random_offset = (random_offset - 0.5) * Ra / 10

        ! Calculate theta for a logarithmic spiral
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
    
    open(unit=10, file='particledata.csv', status="replace", action="write")

    ! Write header
    write(10, '(A)') "x,y,z,v_x,v_y,v_z,a_x,a_y,a_z"

    ! Write data
    do i = 1, N
        write(10, '(9(F12.6, ","))') particle_arr(:, i)
    end do

    ! Close the file
    close(10)
end subroutine initiate_particles

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
    
        !calculate lower bound of the particle
        x_lower = x_min + (ix-1) * (2*dx)
        y_lower = y_min + (iy - 1) * (dy*2.0)
        z_lower = z_min + (iz - 1) * (dz*2.0)

        !compute relative positions
        x_rel = (x - x_lower) / (2*dx)
        y_rel = (y - y_lower) / (2*dy)
        z_rel = (z - z_lower) / (2*dz)
    
        print *, "Particle position:", x, y, z 
        print *, "Grid index", ix, iy, iz
        print *, "Relative position:", x_rel, y_rel, z_rel 
        print *, "__________________________"
    
        !calculate weights


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

program nbody_sim
    use precision
    use cufft_interface
    implicit none
    integer, parameter::N = 1000
    integer, parameter:: nx =20 , ny = 20, nz = 20
    integer:: checkpoint,steps,k
    real:: smbh_m
    real, dimension(9,N)::particles, particle_arr
    real, parameter::dt = 10e-7 ! Needed to keep Energy change way below 10^-5
    real:: E_0,E,Rm,Vm,t_c,curr_time,Rm_0,anim_time, dx, dy, dz, density(nx, ny, nz)
    real,dimension(3)::p
    logical::animate

    dx = 1.0/(nx-1)
    dy = 1.0/(ny-1)
    dz = 1.0/(nz-1)

    smbh_m = 1 ! just set to whatever it is
    animate = .true.
    checkpoint = 10 !s
    
    call initiate_particles(particles,N,1)

    ! initial E_0
    print*, 'Got past initialization'
    call check_energy(density,nx,ny,nz,particles,N,smbh_m,E_0)
    print*, 'Got past check energy'


    ! These 2 will go inside a do loop until end condition
    call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
    print*, 'Got past particle to grid'
    ! add an if for however many steps 
    call check_energy(density,nx,ny,nz,particles,N,smbh_m,E)
    print*, 'Got past second energy check'

    call integration_step(density, nx, ny, nz, particles, N, dt)
    print*, 'Got past integration step energy check'

    ! need a step to write out
 
end program nbody_sim 
