


subroutine compute_accelerations(density_accel_grid,nx,ny,nz,particles,N)
    !###########################################################
    ! Instructions:
    !      Pass in a densitry grid to fill density_accel_grid
    !      that will be modified inplace to become accelerations
    !      in a cubic cell region then mapped onto the particles 
    !      corresponding to that region
    !      
    !      nx,nz,ny are grid dimensions
    !
     !##########################################################

    ! cuda fft lib requirements
    use cudafor
    use cfft
    implicit none
    
    !########################
    !   Host Initializations
    !###########################
    
    integer,intent(in)::  nx,ny,nz,N
    
    ! I can operate in place on the density grid since it needs to be
    ! recomputed anyways, this way avoids creating another huge 3d array of grid_dim^3 cells
    ! atleast on the cpu
    real, Dimension( nx,ny,nz) :: density_accel_grid
    
    real, Dimension(9,N)::particles

    ! Constants

    real :: G = 1 ! Natural Units
    real,parameter:: pi = atan(1.0)*4 
    real:: factor,constants

    ! Iteration
    integer::d1,d2,d3
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K
    ! Wave number stuff

    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition

    real::del_kx,del_ky,del_kz
    
    del_kx =  2*pi/nx
    del_ky = 2*pi/ny
    del_kz = 2*pi/nz

    !#######################################
    !   Forward FFT
    !#######################################

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan


    ! Allocate input and output arrays on the device memory (gpu)
    real, Dimension(:,:,:), allocatable, device :: gpu_accel_grid_r
    allocate(gpu_accel_grid_r(nx,ny,nz))
    
    complex, Dimension(:,:,:), allocatable,device:: gpu_accel_grid_c

    ! https://docs.nvidia.com/cuda/cufft/index.html#multidimensional-transforms
    allocate(gpu_accel_grid_c(nx,ny,nz/2 +1))

    ! 3D R2C Fourier Transform plan setup
    status = cufftPlan3d(plan, nx,ny,nz,CUFFT_R2C)

    if (status .ne. CUFFT_SUCCESS)then
        print*, "Error creating R2C GPU Plan :", status
        stop
    endif
    ! 3D R2C Fourier Transform execution
    status = cufftExecR2C(plan,gpu_accel_grid_r,gpu_accel_grid_c)

    if (status .ne. CUFFT_SUCCESS)then
        print*, "Error executing R2C GPU Plan :", status
        stop
    endif

    !######################################################
    !Compute Gravitational Accelerations in Fourier Space
    !#################################################

    ! The transformed density grid in the reduced dimensions
    ! is then mapped back onto itself as the acceleration grid
    ! replacing the density grid it once was
    ! this saves memory since we need to recompute density each
    ! step anyways
    

    constants =  -4*pi*G
    do k_x=0,nx
        do k_y=0,ny
            do k_z=0,nx/2
                

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

            
                gpu_accel_grid_c(k_x) = k1*gpu_accel_grid_c(k_x) *K*constants
                gpu_accel_grid_c(k_y) = k2*gpu_accel_grid_c(k_y) *K*constants
                gpu_accel_grid_c(k_z) = k3*gpu_accel_grid_c(k_z) *K*constants
    
            end do
        end do
    end do



    !#######################################
    !   Inverse FFT
    !#######################################

    !Inverse 3D C2R Fourier Transform execution on the Gravity Cube
    status = cufftExecC2R(plan,gpu_accel_grid_c,gpu_accel_grid_r)

    if (status .ne. CUFFT_SUCCESS)then
        print*, "Error executing C2R :", status
        stop
    endif
    

    ! Replace density_accel_grid with accels on host memory

    density_accel_grid = gpu_accel_grid_r


    ! Normalize Gravity Cube in real space(divide by N/)
    
    factor = 1/N ! precompute to do multiplication instead of division on array ops

    density_accel_grid = density_accel_grid*factor

    !Destroy Plan
    status = cufftDestroy(plan)
    if (status /= CUFFT_SUCCESS) then
        print *, 'Error destroying cuFFT plan'
        stop
    end if

    ! ################################
    ! Update particles accelerations
    ! ###############################

    ! Choose a way of mapping cube of force back to particles
    ! inside cube, since density array has less
    ! in my code with same dimension this was easy
    ! just particles(7:9,:) = accel_array
    ! but that only works on same size arr

    particles(7:9,:) = ! mapped array (or do in-place)
    

    
end subroutine compute_accelerations


subroutine integration_step(density_grid,nx,ny,nz, particles, N, dt)
    implicit none

    integer,intent(in)::  nx,ny,nz,N
    real, Dimension( nx,ny,nz) :: density_grid

    integer :: i, j
    real, dimension(9, N) :: particles
    ! Sub-arrays for clarity
    real, dimension(3, N) :: x, v, a

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

    open(unit = 10, file='densityfield.csv', status="replace", action="write")

    write(10, '(A)') "x,y,z,density"

    do k = 1, nz
        z = (k - 1.0) * dz
        do j = 1, ny
            y = (j - 1.0) * dy
            do i = 1, nx
                x = (i - 1.0) * dx
                write(10, '(F8.3, ",", F8.3, ",", F8.3, ",", F10.5)') x, y, z, density(i, j, k)
            end do
        end do
    end do

    close(10)

end subroutine particle_to_grid
 
program nbody_sim
    implicit none
    integer, parameter::N = 1000
    integer, parameter:: nx =20 , ny = 20, nz = 20
    integer:: checkpoint,steps,k
    real, dimension(9,N)::particles, particle_arr
    real, parameter::dt = 10e-7 ! Needed to keep Energy change way below 10^-5
    real:: E_0,E,Rm,Vm,t_c,curr_time,Rm_0,anim_time, dx, dy, dz, density(nx, ny, nz)
    real,dimension(3)::p
    logical::animate

    dx = 1.0/(nx-1)
    dy = 1.0/(ny-1)
    dz = 1.0/(nz-1)

    animate = .true.
    checkpoint = 10 !s

    
    call initiate_particles(particles,N,1)


    ! These 2 will go inside a do loop until end condition
    call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
    
    call integration_step(density, nx, ny, nz, particles, N, dt)
    ! need a step to monitor energy
    ! need a step to write out
 
end program nbody_sim 