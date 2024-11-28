!#########################################
!Compile with
!    pgf95 -g -mp -Mcuda=fastmath,cc35,cc50,cc60,fma,unroll,flushz,lineinfo -ta=nvidia -tp=haswell -O2 -Minfo=all -mcmodel=medium Galaxy_Collison.f95 -L/usr/local/cuda-9.2/lib64 -lcufft -lcupti
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
    complex, Dimension(:,:,:),device:: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:),device:: gravity_grid_c_d
    real :: G = 1 ! Natural Units
    real,parameter:: pi = atan(1.0)*4 
    real:: constants
    complex:: p_term
    integer,value::  nx,ny,nz
    complex, parameter :: i_c = (0.0, 1.0)  ! complex i
    ! Iteration
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag

    ! Wave number stuff

    ! From dividing a full 2pi wave over the length of each cell
    ! we get these deltas https://en.wikipedia.org/wiki/Wave_vector#Definition

    real::del_kx,del_ky,del_kz
    real::epsilon
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
    
    
    epsilon = 10e-12
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
                k1 = del_kx*(k_x - nx)
            endif
            if (k_y < ny/2) then
                k2 = del_ky*(k_y)
            else
                k2 = del_ky*(k_y - ny)
            endif
            
            ! z freq always > 0

            k3 = del_kz*k_z


            !k * p(k) *(-4)*pi*G/|K|^2

            K = 1/(abs(k1)**2 + abs(k2)**2 + abs(k3)**2 + epsilon)



            !compute once
            p_term = density_grid_c_d(k_x, k_y, k_z) * K * constants

            ! Sets x,y,z accelerations in fourier space on device grid
            gravity_grid_c_d(1, k_x, k_y, k_z) = i_c*k1 * p_term
            gravity_grid_c_d(2, k_x, k_y, k_z) = i_c*k2 * p_term
            gravity_grid_c_d(3, k_x, k_y, k_z) = i_c*k3 * p_term
        endif
end subroutine compute_gravities

attributes(global) subroutine normalize3d_and_shift(gravity_grid_r_d,gravity_grid_r_d_shifted,nx,ny,nz,factor)
    implicit none
    integer,value::nx,ny,nz
    real, dimension(:,:,:,:),device::gravity_grid_r_d
    real, dimension(:,:,:,:),device::gravity_grid_r_d_shifted
    integer::i,j,K,k_shifted
    real,value::factor
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z


    if ( i<= nx .and. j<=ny .and. k <= nz) then
        
        k_shifted = k + nz/4

        gravity_grid_r_d_shifted(1,i,j,k_shifted) = gravity_grid_r_d(1,i,j,k)*factor

        gravity_grid_r_d_shifted(2,i,j,k_shifted) = gravity_grid_r_d(2,i,j,k)*factor

   

        gravity_grid_r_d_shifted(3,i,j,k_shifted) = gravity_grid_r_d(3,i,j,k)*factor
    endif
    call syncthreads

end subroutine
attributes(global) subroutine calculate_U(density_grid_c_d,nx,ny,nz,U)
    use cudafor
    implicit none
    complex, Dimension(:,:,:),device:: density_grid_c_d
    real,parameter:: pi = atan(1.0)*4 
    integer,value::  nx,ny,nz

    ! Iteration
    integer::k_x,k_y,k_z
    complex::k1,k2,k3
    real::K, p_mag,epsilon
    real,device::U
    real, shared :: U_shared

    integer::istat


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
    U_shared = 0.0
    call syncthreads


    epsilon = 10e-12    

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
                k1 = del_kx*(k_x - nx)
            endif
            if (k_y < ny/2) then
                k2 = del_ky*k_y
            else
                k2 = del_ky*(k_y - ny)
            endif
        
            ! z freq always > 0

            k3 = del_kz*k_z


            !k * p(k) *(-4)*pi*G/|K|^2

            K = 1/(abs(k1)**2 + abs(k2)**2 + abs(k3)**2 + epsilon)

            if (k_x > 0) then       
                ! may have errors, debugging
                p_mag = density_grid_c_d(k_x,k_y,k_z)
                U_shared  = U_shared  - p_mag*K
            end if
        endif
    call syncthreads
    if (threadIdx%x == 1) then
        istat = atomicAdd(U, U_shared)
    end if
end subroutine calculate_U
attributes(global) subroutine integration_step(particles_d, N,dt)
    implicit none
    integer :: i
    integer, value::N
    real, dimension(:,:):: particles_d
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
    call syncthreads

end subroutine integration_step
attributes(global) subroutine calculate_KE(particles_d, N,m,smbh1_m,smbh2_m,KE)

    real, dimension(:,:),device:: particles_d
    integer,value::N
    real,value:: m,smbh1_m,smbh2_m

    real,device::KE
    real, shared :: KE_shared

    integer::istat

    real,dimension(3),device::v_i
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    
    KE_shared = 0.0
    call syncthreads
    if (i==1) then
        
        ! Add KE for Supermassive seperately 
        ! assuming it is particle 1 (index 1)

        ! get velocities
        v_i = particles_d(4:6,1)
        KE_shared = KE_shared + smbh1_m* 0.5* sum(v_i**2)
    else if(i==2)then
        ! Add KE for Supermassive seperately 
        ! assuming it is particle 2 (index 2)
        ! get velocities
        v_i = particles_d(4:6,2)
        KE_shared = KE_shared + smbh2_m* 0.5* sum(v_i**2)
    else if (i<=N) then
        ! get velocities

        v_i = particles_d(4:6,i)

        KE_shared = KE_shared + m* 0.5* sum(v_i**2)
    
    endif
    call syncthreads
    if (threadIdx%x == 1) then
        istat = atomicAdd(KE, KE_shared)
    end if

end subroutine calculate_KE
! A in works Subroutine that does particle to gid on GPU
attributes(global) subroutine particle_to_grid_cuda(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    use cudafor
    implicit none
    integer, value :: N, nx, ny, nz, nx2, ny2, nz2
    real(kind(0.0)), value :: dx, dy, dz,smbh1_m, smbh2_m
    real(kind(0.0)),dimension(:,:),device:: particles_d
    real,dimension(:,:,:),device::density_grid_r_d

    ! Thread and block indices
    integer :: idx, ix, iy, iz, thread_id,istat
    real(kind(0.0)) :: x, y, z, m
    real(kind(0.0)) :: x_rel, y_rel, z_rel
    real(kind(0.0)) :: wx0, wx1, wy0, wy1, wz0, wz1
    real(kind(0.0)) :: x_min, y_min, z_min, x_max, y_max, z_max, delta, delta_z, x_i, y_j, z_k
    

    ! predefined
    x_min = -15.0
    x_max = 15.0
    y_min = -15.0
    y_max = 15.0
    z_min = -10.0
    z_max = 10.0
    delta = (x_max - x_min) / ((nx/2)-1)
    delta_z = (z_max - z_min) / ((nz/2)-1)

    nx2 = nx/4 !offsets for zero padding
    ny2 = ny/4
    nz2 = nz/4


    !compute global thread ID
    thread_id = (blockIdx%x -1) * blockDim%x + threadIdx%x
    if (thread_id >= N) return

    !read particle positions 
    x = particles_d(1,thread_id)
    y = particles_d(2,thread_id)
    z = particles_d(3,thread_id)


    ! Assign mass based on particle ID
    if (thread_id == 1) then 
        m = smbh1_m
    else if  (thread_id == 2) then
        m = smbh2_m
    else 
        m = 0.5/N 
    end if

    
    ! determine grid cell indicies 
    ix = int(floor((x-x_min)/delta)) + 1
    iy = int(floor((y-y_min)/delta)) + 1
    iz = int(floor((z-z_min)/delta_z)) + 1

    ! if (ix < 1 .or. iy < 1 .or. iz < 1) then 
    !     print *, 'particle out of range'
    ! end if 

    !clamp indecies within bounds 
    !if (ix < 1) ix = 1
    !if (ix >= nx/2) ix = nx/2 -1
    !if (iy < 1) iy = 1
    !if (iy >= ny) iy = ny/2 -1
    !if (iz < 1) iz = 1
    !if (iz >= nz) iz = nz/2 -1

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

    if (wx0 < 0.0 .or. wx1 < 0.0) then
        ix = ix -1 
        x_i = x_min + (ix - 1) * delta
       x_rel = (x-x_i)/delta
       wx0 = 1.0 - x_rel 
        wx1 = x_rel
    end if
    if (wy0 < 0.0 .or. wy1 < 0.0) then
        iy = iy -1 
       y_j = y_min + (iy - 1) * delta
        y_rel = (y-y_j)/delta
        wy0 = 1.0 - y_rel 
        wy1 = y_rel 
    end if
    if (wz0 < 0.0 .or. wz1 < 0.0) then
        iz = iz -1
        z_k = z_min + (iz - 1) * delta_z
        z_rel = (z-z_k)/delta_z
        wz0 = 1.0 - z_rel 
        wz1 = z_rel
    end if


    !if (wx0 < 0.0 .or. wx1 < 0.0) then
    !    ix = ix +2 
    !    x_i = x_min + (ix - 1) * delta
    !    x_rel = (x-x_i)/delta
    !    wx0 = 1.0 - x_rel 
    !    wx1 = x_rel
    !end if
    !if (wy0 < 0.0 .or. wy1 < 0.0) then
    !    iy = iy +2 
    !    y_j = y_min + (iy - 1) * delta
    !   y_rel = (y-y_j)/delta
    !    wy0 = 1.0 - y_rel 
    !    wy1 = y_rel 
    !end if
    !if (wz0 < 0.0 .or. wz1 < 0.0) then
    !    iz = iz +2
    !    z_k = z_min + (iz - 1) * delta_z
    !    z_rel = (z-z_k)/delta_z
    !    wz0 = 1.0 - z_rel 
    !    wz1 = z_rel
    !end if


    ! if (wx0 < 0.0 .or. wx1 < 0.0) then
    !     print*, (x-x_min)/delta, x, x_i, ix
    ! end if
    ! if (wy0 < 0.0 .or. wy1 < 0.0) then
    !     print*, (y-y_min)/delta, y,y_j, iy
    ! end if
    ! if (wz0 < 0.0 .or. wz1 < 0.0) then
    !     print*, (z-z_min)/delta, z,z_k, iz
    ! end if

    ! Update density feiled (atomic operations to prevent race condition)

    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2, iz+nz2), m * wx0 * wy0 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2, iz+nz2), m * wx1 * wy0 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2+1, iz+nz2), m * wx0 * wy1 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2+1, iz+nz2), m * wx1 * wy1 * wz0)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2, iz+nz2+1), m * wx0 * wy0 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2, iz+nz2+1), m * wx1 * wy0 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2, iy+ny2+1, iz+nz2+1), m * wx0 * wy1 * wz1)
    istat = atomicadd(density_grid_r_d(ix+nx2+1, iy+ny2+1, iz+nz2+1), m * wx1 * wy1 * wz1)
    call syncthreads

end subroutine particle_to_grid_cuda

    attributes(global) subroutine grid_to_particle_cuda(acceleration_grid, particles, N, nx, ny, nz, dx, dy, dz,smbh1_m, smbh2_m)
    implicit none

    integer, value :: N,nx,ny,nz
    real(kind(0.0)), value :: dx,dy,dz,smbh1_m,smbh2_m
    real(kind(0.0)), dimension(:,:),device:: particles
    real,dimension(:,:,:,:),device:: acceleration_grid

    !thread and block indecies 
    integer :: i,ix,iy,iz,ix_shifted,iy_shifted,iz_shifted,thread_id
    real(kind(0.0)) :: x,y,z,m
    real(kind(0.0)) :: x_rel,y_rel,z_rel
    real(kind(0.0)) :: wx0, wx1, wy0, wy1, wz0, wz1
    real(kind(0.0)) :: x_min, y_min, z_min, x_max, y_max, z_max, delta, delta_z, x_i, y_j, z_k
    real(kind(0.0)) :: acc_x , acc_y, acc_z

    ! Predefined 
    x_min = -15.0
    x_max = 15.0
    y_min = -15.0
    y_max = 15.0
    z_min = -10.0
    z_max = 10.0
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
    else if  (thread_id == 2) then
        m = smbh2_m
    else 
        m = 0.5/N 
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
    !if (ix < 1) ix = 1
    !if (ix >= nx / 2) ix = nx / 2 - 1
    !if (iy < 1) iy = 1
    !if (iy >= ny / 2) iy = ny / 2 - 1
    !if (iz < 1) iz = 1
    !if (iz >= nz / 2) iz = nz / 2 - 1

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

    if (wx0 < 0.0 .or. wx1 < 0.0) then
        ix = ix -1 
        x_i = x_min + (ix - 1) * delta
        x_rel = (x-x_i)/delta
        wx0 = 1.0 - x_rel 
        wx1 = x_rel
    end if
    if (wy0 < 0.0 .or. wy1 < 0.0) then
        iy = iy -1 
        y_j = y_min + (iy - 1) * delta
        y_rel = (y-y_j)/delta
        wy0 = 1.0 - y_rel 
        wy1 = y_rel 
    end if
    if (wz0 < 0.0 .or. wz1 < 0.0) then
        iz = iz -1
        z_k = z_min + (iz - 1) * delta_z
        z_rel = (z-z_k)/delta_z
        wz0 = 1.0 - z_rel 
        wz1 = z_rel
    end if
    

    !initialize accelerations 
    acc_x = 0.0
    acc_y = 0.0
    acc_z = 0.0

    ix_shifted = ix+nx/4
    iy_shifted = iy+ny/4
    iz_shifted = iz+nz/4
    if (ix_shifted < 1) ix = 1
    if (ix_shifted >= nx / 2) ix = nx /2 -1
    if (iy_shifted < 1) iy = 1
    if (iy_shifted >= ny / 2) iy = ny / 2 - 1
    if (iz_shifted < 1) iz = 1
    if (iz_shifted >= nz / 2) iz = nz / 2 - 1

    ! Interpolate acceleration from the grid to the particle position
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted, iz_shifted) * wx0 * wy0 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted, iz_shifted) * wx1 * wy0 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted + 1, iz_shifted) * wx0 * wy1 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted + 1, iz_shifted) * wx1 * wy1 * wz0
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted, iz_shifted + 1) * wx0 * wy0 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted, iz_shifted + 1) * wx1 * wy0 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted, iy_shifted + 1, iz_shifted + 1)* wx0 * wy1 * wz1
    acc_x = acc_x + acceleration_grid(1, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1) * wx1 * wy1 * wz1

    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted, iz_shifted) * wx0 * wy0 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted, iz_shifted) * wx1 * wy0 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted + 1, iz_shifted) * wx0 * wy1 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted + 1, iz_shifted)* wx1 * wy1 * wz0
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted, iz_shifted + 1) * wx0 * wy0 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted, iz_shifted + 1) * wx1 * wy0 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted, iy_shifted + 1, iz_shifted + 1) * wx0 * wy1 * wz1
    acc_y = acc_y + acceleration_grid(2, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1) * wx1 * wy1 * wz1

    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted, iz_shifted) * wx0 * wy0 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted, iz_shifted) * wx1 * wy0 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted + 1, iz_shifted) * wx0 * wy1 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted + 1, iz_shifted) * wx1 * wy1 * wz0
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted, iz_shifted + 1) * wx0 * wy0 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted, iz_shifted + 1) * wx1 * wy0 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted, iy_shifted + 1, iz_shifted + 1) * wx0 * wy1 * wz1
    acc_z = acc_z + acceleration_grid(3, ix_shifted + 1, iy_shifted + 1, iz_shifted + 1) * wx1 * wy1 * wz1

    ! update particle 
    particles(7, thread_id) = acc_x
    particles(8, thread_id) = acc_y
    particles(9, thread_id) = acc_z
    
    call syncthreads

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
    

module particle_kernels
    implicit none
    
contains
subroutine check_energy(density_grid_r_d,density_grid_c_d,nx,ny,nz,particles_d,N,m,smbh1_m,smbh2_m,E,plan)
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
    use cudafor
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

    real :: G = 1.0 ! Natural Units

    ! Iteration
    integer::i

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan,blockDim,istat

    !###################################
    !   Device Initialization
    !###################################    
    !################################### 

    real,Device::U_d,KE_d
    real,Dimension(:,:),allocatable,device::particles_d
    real,Dimension(:,:,:),allocatable,device::density_grid_r_d
    complex,Dimension(:,:,:),allocatable,device::density_grid_c_d

    
    !#######################################
    !   Forward FFT
    !#######################################
    !print*,"beginning fft"
   


    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)
    !print*,"FINISHED     fft"

    !######################################################
    !Compute Gravitational Potential in Fourier Space
    !#####################################  ############
    
    ! Define block dimensions
    blockDim = 1024
    blockDimX = 8
    blockDimY = 8
    blockDimZ = 8

    gridDimX = (nx + blockDimX - 1) / blockDimX
    gridDimY = (ny + blockDimY - 1) / blockDimY
    gridDimZ = ((nz/2 +1) + blockDimZ - 1) / blockDimZ


 
    ! get Volume of cube
    V = nx*ny*nz

    ! Reset U,KE

    U_d = 0.0
    KE_d = 0.0
    U=0.0
    KE = 0.0
    !print*, "got to just before potential launch"
    call calculate_U<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(density_grid_c_d,nx,ny,nz,U_d)
    istat = cudaDeviceSynchronize()

    !print*, "Calculated potential"
    call calculate_KE<<<(N-1)/blockDim,blockDim>>>(particles_d,N,m,smbh1_m,smbh2_m,KE_d)
    istat = cudaDeviceSynchronize()

    !print*, "Calculated KE"
    U=U_d
    KE = KE_d

    U = (2.0*pi*G/V)*U

    
    ! combine energies
            E = U + KE

    
end subroutine check_energy

subroutine fft_step(density_grid_r_d,density_grid_c_d,gravity_grid_r_d,gravity_grid_r_d_shifted,gravity_grid_c_d,nx,ny,nz,N,plan)
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
    use cudafor
    implicit none
    
    !########################
    !   Host Initializations
    !
    !   Only small administrative things
    !   for composing kernel calls
    !   and cufft functions
    !###########################
    
    integer::  nx,ny,nz, N

   
    ! Normalization Factor

    real :: factor
  

    ! Cuda Variables for plan process
    !cufftHandle plan identifier and error
    integer::status,plan,istat
    
    !########################
    !   Device Initializations
    !#########################


   
    ! Real and complex density on gpu
    real, Dimension(:,:,:), allocatable,device :: density_grid_r_d

    complex, Dimension(:,:,:), allocatable,device:: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:), allocatable,device:: gravity_grid_c_d
    real, Dimension(:,:,:,:), allocatable,device :: gravity_grid_r_d
    real, Dimension(:,:,:,:), allocatable, device,intent(inout):: gravity_grid_r_d_shifted

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

   
    


    ! 3D R2C Fourier Transform execution
    call cufftExecR2C(plan,density_grid_r_d,density_grid_c_d)
    istat = cudaDeviceSynchronize()	


    call compute_gravities<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(gravity_grid_c_d,density_grid_c_d,nx,ny,nz)
	istat = cudaDeviceSynchronize()	

    !print *, "compute_gravities finished"
    !#######################################
    !   Inverse FFT
    !#######################################

    !Inverse 3D C2R Fourier Transform execution on the Gravity Cube
    call cufftExecC2R(plan,gravity_grid_c_d,gravity_grid_r_d)
	istat = cudaDeviceSynchronize()	

    ! !print*, "Density Grid"
    ! !print*, density_grid
    ! !print*, "Gravity Grid"
    ! !print*, gravity_grid
    ! Normalize Gravity Cube in real space(divide by N/)
    

    ! TODO: Check im not crazy and i should be using N
    factor = 1.0/real((nx*ny*nz))
    ! Apply factor ONLY to the acceleration dimensions not the index ones

    gridDimX = (nx + blockDimX - 1) / blockDimX
    gridDimY = (ny + blockDimY - 1) / blockDimY
    gridDimZ = (nz+ blockDimZ - 1) / blockDimZ

    gravity_grid_r_d_shifted = 0.0
    call normalize3d_and_shift<<<[gridDimX, gridDimY, gridDimZ], [blockDimX, blockDimY, blockDimZ]>>>(gravity_grid_r_d,gravity_grid_r_d_shifted,nx,ny,nz,factor)

    istat = cudaDeviceSynchronize()	

    !print*, "inside fft"

    
end subroutine fft_step

    
end module particle_kernels



subroutine initialize_particles2(particles, N, Ra, disk_mass, smbh1_mass, smbh2_mass, R_disk, R_cl, G, rho_c)
    implicit none

    ! Inputs
    integer, intent(in) :: N
    real, intent(in) :: Ra, disk_mass, smbh1_mass, smbh2_mass, R_disk, R_cl, G, rho_c
    real, dimension(9, N), intent(out) :: particles

    ! Local variables
    integer :: i
    real :: r, theta, z, v_c, v_x, v_y, v_z, random_factor

    call random_seed()

    ! Initialize particles
    do i = 1, N
        if (i == 1) then
            ! Central SMBH (M1)
            particles(1:3, i) = (/ 0.0, 0.0, 0.0 /)
            particles(4:6, i) = (/ 0.0, 0.0, 0.0 /)
            particles(7:9, i) = (/ 0.0, 0.0, 0.0 /)

        elseif (i == 2) then
            ! Secondary SMBH (M2)
            particles(1:3, i) = (/ Ra * 0.5, 0.0, 0.0 /)
            particles(4:6, i) = (/ 0.0, sqrt(G * smbh1_mass / max(Ra * 0.5, 1e-6)), 0.0 /)
            particles(7:9, i) = (/ 0.0, 0.0, 0.0 /)

        else
            ! Disk particles
            call random_number(r)
            r = 1+(R_disk-1)*r
            if (r < 1e-6) r = 1e-6  ! Avoid very small values

            call random_number(theta)
            theta = theta * 2.0 * atan(1.0) * 4.0

            call random_number(z)
            z = (z - 0.5) * 2.0 * 0.08 * R_disk

            v_c = sqrt(G * (smbh1_mass + smbh2_mass + compute_M_cl(r, rho_c, R_cl) + disk_mass * r / R_disk) / r)
            if (v_c /= v_c) stop "NaN detected in v_c"

            call random_number(random_factor)
            v_x = -v_c * sin(theta) * (1.0 + (random_factor - 0.5) * 0.16)
            call random_number(random_factor)
            v_y = v_c * cos(theta) * (1.0 + (random_factor - 0.5) * 0.16)
            v_z = 0.0
            
            particles(:, i) = (/ r * cos(theta), r * sin(theta), z, v_x, v_y, v_z, 0.0, 0.0, 0.0 /)
            
            
        end if
    end do


contains
    real function compute_M_cl(r, rho_c, R_cl)
        implicit none
        real, intent(in) :: r, rho_c, R_cl
        real :: term1, term2

        term1 = 0.5 * log(1.0 + (r / max(R_cl, 1e-6))**2)
        term2 = (r / max(R_cl, 1e-6))**2 / (1.0 + (r / max(R_cl, 1e-6))**2)
        compute_M_cl = 4.0 * atan(1.0) * 4.0 * rho_c * R_cl**3 * (term1 + term2)
    end function compute_M_cl
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
    use cudafor
    use device_ops
    use cufft_interface
    use particle_kernels
    implicit none
    integer, parameter::N = 10000
    integer, parameter:: nx =36 , ny = 36, nz = 16  
    real, Dimension(nx,ny,nz):: density_grid_test
    real, Dimension(3,nx,ny,nz):: gravity_grid_test
    real, parameter :: R_disk = 10.0, R_cl = 1.0, Ra = 10.0
    real, parameter :: disk_mass = 0.5, smbh1_m = 0.5, smbh2_m = 0.01
    real, parameter :: G = 1.0  ! Gravitational constant in normalized units
    real :: rho_c

    integer:: checkpoint,steps,k,i,ierr,t,u,v,w
    real:: m,E_0,E,dx,dy,dz
    real, dimension(9,N)::particles
    real, parameter::dt = 0.1 !Needed to keep Energy change way below 10^-5
    real,dimension(3)::p
    logical::animate
    integer:: particle_to_track = 50
    integer :: blockDimX, blockDimY, blockDimZ,blockDim
    integer :: gridDimX, gridDimY, gridDimZ

    integer:: plan,istat


    ! ############################################
    ! DEVICE MEMORY SETUP
    ! #######################################
    ! Particles on GPU
    real, Dimension(:,:),allocatable, device::particles_d
    

    ! Real and complex density on gpu
    real, Dimension(:,:,:), allocatable, device :: density_grid_r_d
    complex, Dimension(:,:,:), allocatable,device:: density_grid_c_d

    ! Real and complex gravities on gpu
    complex, Dimension(:,:,:,:), allocatable,device:: gravity_grid_c_d
    real, Dimension(:,:,:,:), allocatable, device:: gravity_grid_r_d_shifted
    real, Dimension(:,:,:,:), allocatable, device :: gravity_grid_r_d
    character(len=50) :: t_string
    
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
    allocate(gravity_grid_r_d_shifted(3,nx,ny,nz))


    ! allocation particles on device
    allocate(particles_d(9,N))

    !##############################################
    !
    !   HOST(CPU) INITIALIZATIONS
    !
    !##############################################

    
    rho_c = (disk_mass / (4.0 * atan(1.0) * 4.0 * R_cl**3)) * (3.0 / (2.0 * log(2.0)))

    ! Initialize particles
    call initialize_particles2(particles, N, Ra, disk_mass, smbh1_m, smbh2_m, R_disk, R_cl, G, rho_c)


    m = 0.5/N
    E_0 = 0.0
    dx = 1.0/(nx-1) 
    dy = 1.0/(ny-1)
    dz = 1.0/(nz-1)
    !##############################################
    !
    !   Device(GPU) INITIALIZATIONS
    !
    !##############################################
    
    istat = cudaSetDevice(0)



    ! copy particles from host to device memory
    particles_d = particles
    print*, "Particles (only 10)"
        do k = 1, 10
            print*, particles(:,k) 
        enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! beyond this every major step stays in device memory
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initial E_0
    print*, 'Got past initialization'


    blockDim = 1024
    blockDimX = 8
    blockDimY = 8
    blockDimZ = 8
    gridDimX = (nx-1+blockDimX)/blockDimX
    gridDimY =(nx-1+blockDimY)/blockDimY
    gridDimZ = (nx-1+blockDimZ)/blockDimZ
        

    
    ! 3D R2C Fourier Transform plan setup
    call cufftPlan3d(plan, nx,ny,nz,CUFFT_R2C)
    istat = cudaDeviceSynchronize()	

    

    print*, 'Got past particle to grid'
    ! call particle_to_grid_cuda<<<blockDim,blockDim>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    call particle_to_grid_cuda<<<(N+blockDim-1)/blockDim,blockDim>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
    istat = cudaDeviceSynchronize()	

    
    !print*, 'Got past check energy - lol no'
    

    do i=1, 200

        if (i==1) then
            call check_energy(density_grid_r_d,density_grid_c_d,nx,ny,nz,particles_d,N,m,smbh1_m,smbh2_m,E_0,plan)
            print*, E_0
        endif
        ! These 2 will go inside a do loop until end condition
        ! call particle_to_grid_cuda<<<blockDim,blockDim>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
        density_grid_r_d = 0.0
        istat = cudaDeviceSynchronize()	

        call particle_to_grid_cuda<<<(N+blockDim-1)/blockDim,blockDim>>>(density_grid_r_d, particles_d, N, nx, ny, nz, dx, dy, dz,smbh1_m,smbh2_m)
        istat = cudaDeviceSynchronize()	


        !density_grid_test = density_grid_r_d
        
        !print*, "Density Cube:"
        !print*, density_grid_test
        ! add an if for however many steps 
        ! again, like fft stays on gpu but composes with a fft call
        call check_energy(density_grid_r_d,density_grid_c_d,nx,ny,nz,particles_d,N,m,smbh1_m,smbh2_m,E,plan)


        !print*, 'Got past second energy check'


        ! fills the real gravity grid
        ! still stays on gpu
        ! host/cpu - style function is just composing cuda kernel functions

        
        call fft_step(density_grid_r_d,density_grid_c_d,gravity_grid_r_d,gravity_grid_r_d_shifted,gravity_grid_c_d, nx,ny,nz,N,plan)
        istat = cudaDeviceSynchronize()	

   
        !print*, "Got past grid to particle"
       

        !! here zac call your grid to particles kernel
        !! heres and example you can change dimensions and stuff
        
        call grid_to_particle_cuda<<<(N+blockDim-1)/blockDim,blockDim>>>(gravity_grid_r_d_shifted,particles_d,N,nx, ny, nz,dx, dy, dz,smbh1_m,smbh2_m)
        istat = cudaDeviceSynchronize()	

        ! Check for errors
        
        ! print*, "after grid to p Particles (only 10)"
        ! particles = particles_d

        ! do k = 1, N
        !     print*, particles(:,k) 
        ! end do

      
      
        ! integration step pushes all positions
        ! ill need to revisit thread count block size just going quick
        ! to get structure
        call integration_step<<<(N+blockDim-1)/blockDim,blockDim>>>(particles_d,N,dt)
        istat = cudaDeviceSynchronize()	

     
        print*, "Done step: ", i
        print*, "Relative Energy Change: ",(E-E_0)/E_0
        ! need a step to [pass back & write out

    end do

    ! Deallocations

    particles = particles_d
    print*, "Particles (only 10)"
        do k = 1, 10
            print*, particles(:,k) 
        enddo

    !deallocate density grids
    deallocate(density_grid_r_d,density_grid_c_d)
    

    ! deallocate gravity grids
    deallocate(gravity_grid_r_d,gravity_grid_c_d,gravity_grid_r_d_shifted)

    ! deallocate particles
    deallocate(particles_d)
        !Destroy Plan
    call cufftDestroy(plan)

 
end program nbody_sim