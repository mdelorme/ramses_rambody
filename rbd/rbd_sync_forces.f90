subroutine rbd_sync_forces
  use amr_commons
  use rbd_commons
  use pm_commons
  use rbd_parameters
  implicit none
  
  include 'mpif.h'
  
  integer :: ierr, i, j, k
  integer :: nbufpart
  integer :: n_particles
  real(dp), allocatable, dimension(:,:) :: dbl_buf
  integer,  allocatable, dimension(:)   :: int_buf
  real(dp), allocatable, dimension(:)   :: dbg_tbl
  real(dp), dimension(1:3) :: force, dpt
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2, scale_M
  real(dp), parameter :: cm2pc  = 3.24078D-19 
  real(dp), parameter :: s2Myr  = 3.16881D-14
  real(dp) :: conv_factor
  real(dp) :: min_ext, max_ext, nrm

  real(dp) :: tot, dist, c, tmp, to_cgs, to_pc, nb6_mesh_scale
  real(dp), dimension(3) :: ri, ai
  integer :: total_part, pid, ilevel, jpart, icpu, igrid, jgrid, ipart, check_tot
  integer, dimension(1:3) :: mesh_desc

  real(dp) :: fr0, min_frc, max_frc, fr, min_dbg, max_dbg, fa0
  real(dp) :: rij
  real(dp), dimension(1:3) :: dbg_rij
  
  logical, allocatable, dimension(:) :: dbg_ids
  character(len=1024) :: dbg_fn

  ! Retrieving units
  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
  scale_M = scale_d * scale_l**3.0

  ! /\ Force in Ramses is an acceleration, not a real force !
  ! Force is in Ramses units, so we convert them to CGS then to Nbody physical units
  to_cgs = scale_l / scale_T**2.0
  to_pc  = cm2pc   / s2Myr**2.0
  
  rbd_send_buf = 0
  check_tot = 0

  allocate(dbl_buf(3, rbd_mesh_np))
  allocate(int_buf(rbd_mesh_np))

  
  rbd_fc = 0
  rbd_dt = 0
  c = boxlen * 0.5

  ! If the current process is the root (of Rambody), then we gather all the information from the other processes
  ! into one array, that we send to NBody 6
  ! Otherwise, the info on current process is sent to the root
  !
  ! All data is converted to Nbody6's physical units (pc, Msun, Myr) before being sent.
  if (myid .eq. 1) then
     n_particles = 0 

     ! Filling Rambody comm buffer with CPU 1's own particles
     do ilevel=1, nlevelmax
        do icpu=1, ncpu
           igrid=headl(icpu, ilevel)
           do jgrid=1,numbl(icpu, ilevel) 
              ipart = rbd_headp(igrid)
              check_tot = check_tot + rbd_id(ipart)

              do jpart=1, rbd_numbp(igrid)
                 n_particles = n_particles + 1

                 ! Storing position of the particles
                 rbd_send_buf(1:3, rbd_id(ipart)) = rbd_mesh_pos(:, rbd_id(ipart)) - rbd_mesh_pos(:,1)
                 rbd_send_buf(4:6, rbd_id(ipart)) = rbd_fp(:, ipart)

                 ! Duplicating the info for output dumping
                 rbd_mesh_force(1:3, rbd_id(ipart)) = rbd_fp(:, ipart)

                 ipart = rbd_nextp(ipart)
              end do

              igrid = next(igrid)
           end do
        end do
     end do

     ! TODO : Group this
     call MPI_AllReduce(MPI_IN_PLACE, rbd_fc, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_RAMSES, ierr)
     call MPI_AllReduce(MPI_IN_PLACE, rbd_dt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_RAMSES, ierr)

     ! Getting the rest of the mesh from other processes
     do i=2, ncpu
        ! Number of particles, ids and forces
        call MPI_Recv(n_particles, 1, MPI_INT, i-1, 0, MPI_COMM_RAMSES, MPI_STATUS_IGNORE, ierr)
        call MPI_Recv(int_buf, n_particles, MPI_INT, i-1, 1, MPI_COMM_RAMSES, MPI_STATUS_IGNORE, ierr)
        call MPI_Recv(dbl_buf, 3*n_particles, MPI_DOUBLE_PRECISION, i-1, 2, MPI_COMM_RAMSES, MPI_STATUS_IGNORE, ierr)

        ! Then we group everything in the send buffer
        do j=1, n_particles
           pid = int_buf(j)

           rbd_send_buf(1:3,pid) = rbd_mesh_pos(:,pid) - rbd_mesh_pos(:,1)
           rbd_send_buf(4:6,pid) = dbl_buf(:,j)
           ! Duplicating the info for output dumping
           rbd_mesh_force(1:3, pid) = dbl_buf(:, j)
        end do
     end do

     ! Converting to Nbody6 physical units
     do i=1, rbd_mesh_np
        rbd_send_buf(1:3,i) = rbd_send_buf(1:3,i) * scale_l * cm2pc
        rbd_send_buf(4:6,i) = rbd_send_buf(4:6,i) * to_cgs  * to_pc
     end do

     call timer('rbd_wait_force', 'start')
     ! And sending everything to Nbody6
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 16 (-> Mesh scale)'
        call flush(6)
     end if

     nb6_mesh_scale = rbd_mesh_scale * scale_l * cm2pc
     call MPI_Send(nb6_mesh_scale, 1, MPI_DOUBLE_PRECISION, 0, 16, MPI_COMM_RAMBODY, ierr)
     
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 9 (-> N mesh particles)'
        call flush(6)
     end if
     
     call MPI_Send(rbd_mesh_np, 1, MPI_INTEGER, 0, 9, MPI_COMM_RAMBODY, ierr)
     
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 10 (-> Mesh particles)'
        call flush(6)
     end if

     call MPI_Send(rbd_send_buf, rbd_mesh_np*6, MPI_DOUBLE_PRECISION, 0, 10, MPI_COMM_RAMBODY, ierr)

     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 12 (-> Mesh description)'
        call flush(6)
     end if

     ! TODO : Redo this with nbody6
     mesh_desc(1) = rbd_mesh_np
     mesh_desc(2) = rbd_mesh_Nx
     mesh_desc(3) = 0
     
     call MPI_Send(mesh_desc, 3, MPI_INTEGER, 0, 12, MPI_COMM_RAMBODY, ierr)

     ! Debug comm
     rbd_dbg_buf(:,1) = (rbd_mesh_pos(:,1) - c) * scale_l * cm2pc 

     call MPI_Send(rbd_dbg_buf, 3*n_dbg_buf, MPI_DOUBLE_PRECISION, 0, 99, MPI_COMM_RAMBODY, ierr)

     ! Receiving now the force from Nbody6
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 17 (<- Nbody force)'
        call flush(6)
     end if
     call MPI_Recv(rbd_mesh_force, rbd_mesh_np*3, MPI_DOUBLE_PRECISION, 0, 17, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)

     ! Converting to Ramses units
     rbd_mesh_force(:,1:rbd_mesh_np) = rbd_mesh_force(:,1:rbd_mesh_np) / to_cgs / to_pc
  else
     ! All the other processes only gather the force and send it to the root process
     n_particles = 0

     ! Filling in the communicator
     do ilevel=1, nlevelmax
        do icpu=1, ncpu
           igrid=headl(icpu, ilevel)
           do jgrid=1,numbl(icpu, ilevel) 
              ipart = rbd_headp(igrid)
              do jpart=1, rbd_numbp(igrid)
                 n_particles = n_particles + 1
                 int_buf(n_particles) = rbd_id(ipart)
                 dbl_buf(:, n_particles) = rbd_fp(:, ipart)
                 ipart = rbd_nextp(ipart)
              end do

              igrid = next(igrid)
           end do
        end do
     end do

     ! And communicating everything 
     ! TODO : Group this
     call MPI_AllReduce(MPI_IN_PLACE, rbd_fc, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_RAMSES, ierr)
     call MPI_AllReduce(MPI_IN_PLACE, rbd_dt, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_RAMSES, ierr)
     
     call MPI_Send(n_particles, 1, MPI_INT, 0, 0, MPI_COMM_RAMSES, ierr)
     call MPI_Send(int_buf, n_particles, MPI_INT, 0, 1, MPI_COMM_RAMSES, ierr)
     call MPI_Send(dbl_buf, n_particles*3, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_RAMSES, ierr)
     call timer('rbd_wait_force', 'start')
  end if

  ! And broadcasting the final force mesh
  write(6,*) 'RBD : Broadcasting mesh force to everyone'
  call MPI_Bcast(rbd_mesh_force, rbd_mesh_np*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
  call timer('rbd_sync', 'start')
  
  deallocate(dbl_buf)
  deallocate(int_buf)

  ! This is just for testing/publication purposes
  !call rbd_build_force_profile
end subroutine rbd_sync_forces

subroutine rbd_sync_forces_restart
  use amr_commons
  use rbd_commons
  use rbd_parameters
  implicit none
  
  include 'mpif.h'
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2, scale_M
  real(dp), parameter :: cm2pc  = 3.24078D-19 
  real(dp), parameter :: s2Myr  = 3.16881D-14
  real(dp) :: to_cgs, to_pc
  integer, dimension(1:3) :: mesh_desc
  integer :: i, ierr

  ! This should be called only by the root process in Ramses !

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
  scale_M = scale_d * scale_l**3.0

  to_cgs = scale_l / scale_T**2.0
  to_pc  = cm2pc   / s2Myr**2.0

  if (rbd_dbg_comm) then
     write(6,*) 'RBD : Comm tag 13 (-> Restart mesh)'
  end if

  call MPI_Send(rbd_mesh_np, 1, MPI_INTEGER, 0, 13, MPI_COMM_RAMBODY, ierr)


  do i=1, rbd_mesh_np
     rbd_send_buf(1:3, i) = (rbd_mesh_pos(:, i) - rbd_mesh_pos(:, 1))  * scale_l * cm2pc
     rbd_send_buf(4:6, i) = rbd_mesh_force(:, i) * to_cgs * to_pc
  end do


  if (rbd_dbg_comm) then
     write(6,*) 'RBD : Comm tag 14 (-> Mesh position and forces)'
  end if

  call MPI_Send(rbd_send_buf, rbd_mesh_np*6, MPI_DOUBLE_PRECISION, 0, 14, MPI_COMM_RAMBODY, ierr)

  if (rbd_dbg_comm) then
     write(6,*) 'RBD : Comm tag 15 (-> Mesh description)'
  end if

  mesh_desc(1) = rbd_mesh_np
  mesh_desc(2) = rbd_mesh_Nx
  mesh_desc(3) = 0

  call MPI_Send(mesh_desc, 3, MPI_INTEGER, 0, 15, MPI_COMM_RAMBODY, ierr)

end subroutine
