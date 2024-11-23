
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
end subroutine initiate_particles

subroutine particle_to_grid(density_full, particles, N, nx, ny, nz, dx, dy, dz)
    use cudafor
    implicit none
    integer, intent(in) :: N, nx, ny, nz
    real, intent(in) :: particles(9, N), dx, dy, dz
    real, intent(out) :: density_full(nx, ny, nz)
  

    real, device :: d_particles(9, N)
    real, device :: d_density(nx/2, ny/2, nz/2)
  

    real :: density(nx/2, ny/2, nz/2)
  

    integer :: threadsPerBlock, blocksPerGrid
    integer :: i, j, k, istart, iend, jstart, jend, kstart, kend
    integer :: nx_inner, ny_inner, nz_inner
    integer :: id, jd, kd
  

    density = 0.0

    d_particles = particles
  

    d_density = 0.0
  

    threadsPerBlock = 256
    blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock
  

    call particle_to_grid_kernel<<<blocksPerGrid, threadsPerBlock>>>(d_density, d_particles, N, nx, ny, nz)
  

    density = d_density
  

    nx_inner = nx / 2
    ny_inner = ny / 2
    nz_inner = nz / 2
  
    istart = nx / 4 + 1
    iend = istart + nx_inner - 1
    jstart = ny / 4 + 1
    jend = jstart + ny_inner - 1
    kstart = nz / 4 + 1
    kend = kstart + nz_inner - 1
  
    density_full = 0.0
  
    do k = kstart, kend
      kd = k - kstart + 1
      do j = jstart, jend
        jd = j - jstart + 1
        do i = istart, iend
          id = i - istart + 1
          density_full(i, j, k) = density(id, jd, kd)
        end do
      end do
    end do

    !output to csv
  
    open(unit = 10, file='densityfield.csv', status="replace", action="write")

    write(10, '(A)') "x,y,z,density"

    do k = 1, nz
        z = k
        do j = 1, ny
            y = j
            do i = 1, nx
                x = i
                write(10, '(F8.3, ",", F8.3, ",", F8.3, ",", F20.8)') x, y, z, density_full(i, j, k)
            end do
        end do
    end do

    close(10)

    open(unit = 10, file='densityfield2.csv', status="replace", action="write")

    write(10, '(A)') "x,y,z,density"

    do k = 1, nz/2
        z = k
        do j = 1, ny/2
            y = j
            do i = 1, nx/2
                x = i
                write(10, '(F8.3, ",", F8.3, ",", F8.3, ",", F20.8)') x, y, z, density(i, j, k)
            end do
        end do
    end do

    close(10)
  
  end subroutine particle_to_grid
  
  attributes(global) subroutine particle_to_grid_kernel(d_density, d_particles, N, nx, ny, nz)
    use cudafor
    implicit none
    integer, value :: N, nx, ny, nz
    real, device :: d_density(nx/2, ny/2, nz/2)
    real, device :: d_particles(9, N)
  
    integer :: idx
    integer :: ix, iy, iz
    real :: x, y, z, m
    real :: x_rel, y_rel, z_rel
    real :: wx0, wx1, wy0, wy1, wz0, wz1
    real :: x_min, y_min, z_min, x_max, y_max, z_max, delta, delta_z
    real :: x_i, y_j, z_k
  
    idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
  
    if (idx <= N) then
      m = 1.0  ! Assign mass to each particle
  
      x = d_particles(1, idx)
      y = d_particles(2, idx)
      z = d_particles(3, idx)
  
      x_min = -1.5
      x_max = 1.5
      y_min = -1.5
      y_max = 1.5
      z_min = -1.5
      z_max = 1.5
      delta = (x_max - x_min) / ((nx/2)-1)
      delta_z = (z_max - z_min) / ((nz/2)-1)
  

      if (x >= x_min .and. x <= x_max .and. y >= y_min .and. y <= y_max .and. z >= z_min .and. z <= z_max) then
        ix = int(floor((x - x_min) / delta)) + 1
        iy = int(floor((y - y_min) / delta)) + 1
        iz = int(floor((z - z_min) / delta_z)) + 1
  

        ix = max(1, min(ix, nx/2 - 1))
        iy = max(1, min(iy, ny/2 - 1))
        iz = max(1, min(iz, nz/2 - 1))
  
        x_i = x_min + (ix - 1) * delta
        y_j = y_min + (iy - 1) * delta
        z_k = z_min + (iz - 1) * delta_z
  

        x_rel = (x - x_i) / delta
        y_rel = (y - y_j) / delta
        z_rel = (z - z_k) / delta_z
  

        wx0 = 1.0 - x_rel
        wx1 = x_rel
        wy0 = 1.0 - y_rel
        wy1 = y_rel
        wz0 = 1.0 - z_rel
        wz1 = z_rel
  

        call atomicAdd(d_density(ix, iy, iz), m * wx0 * wy0 * wz0)
        call atomicAdd(d_density(ix+1, iy, iz), m * wx1 * wy0 * wz0)
        call atomicAdd(d_density(ix, iy+1, iz), m * wx0 * wy1 * wz0)
        call atomicAdd(d_density(ix+1, iy+1, iz), m * wx1 * wy1 * wz0)
        call atomicAdd(d_density(ix, iy, iz+1), m * wx0 * wy0 * wz1)
        call atomicAdd(d_density(ix+1, iy, iz+1), m * wx1 * wy0 * wz1)
        call atomicAdd(d_density(ix, iy+1, iz+1), m * wx0 * wy1 * wz1)
        call atomicAdd(d_density(ix+1, iy+1, iz+1), m * wx1 * wy1 * wz1)
      end if
    end if
  end subroutine particle_to_grid_kernel
 
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
    call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)

 
end program nbody_sim 