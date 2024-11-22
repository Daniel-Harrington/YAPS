
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

subroutine particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)
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
    call particle_to_grid(density, particles, N, nx, ny, nz, dx, dy, dz)

 
end program nbody_sim 