subroutine rbd_sync_cluster
  use amr_commons
  use rbd_commons
  use rbd_parameters
  use mpi_mod
  implicit none
  
  
  integer :: ierr, i, k, ipart, nx_loc, ilevel
  integer :: nbufpart
  integer, dimension(1:nvector) :: cc
  real(dp) :: mtot
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp), parameter :: pc2cm   = 3.086e18
  real(dp), parameter :: kms2cms = 1.0e5
  real(dp), parameter :: M2g     = 1.988e33
  real(dp), dimension(1:3) :: com ! center of mass
  real(dp) :: min_ext, max_ext, nrm, dx, dx_loc, scale
  real(dp), dimension(1:3, 1:nvector) :: xx_dp

  integer :: n_escapers

  logical :: debug_force = .false.
  integer, parameter :: np = 100
  real, parameter :: pi = acos(-1.0)
  real, parameter :: da = 2.0 * pi / np
  real :: angle = 0.0
  
  ! First we synchronise escapers
  call rbd_sync_escapers

  ! Retrieving units
  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  call timer('rbd_wait_cluster', 'start')

  ! Only getting the cluster on cpu 1
  if (myid == 1) then     
     write(6,*) 'RBD : Getting particle count'
     call flush(6)
     ! How many particles should we receive ?
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 6 (<- N_particles)'
        call flush(6)
     end if

     call MPI_Recv(nbufpart, 1, MPI_INTEGER, 0, 6, MPI_COMM_RAMBODY,&
          & MPI_STATUS_IGNORE, ierr)

     ! Big comm to get all the particles in the cluster
     rbd_recv_buf = 0
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 7 (<- NB6 particles)'
        call flush(6)
     end if

     call MPI_Recv(rbd_recv_buf, nbufpart*7, MPI_DOUBLE_PRECISION, 0,&
          & 7, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)

     ! Mass sum
     mtot = 0.0
     do i=1, nbufpart
        mtot = mtot + rbd_recv_buf(7, i)

	do k=1,3
          nb6_xp(k,i) = rbd_recv_buf(k,i)   * pc2cm   / scale_l !+ rbd_xc(k)
          nb6_vp(k,i) = rbd_recv_buf(k+3,i) * kms2cms / scale_v + rbd_vc(k)
       end do

        nb6_mp(i) = rbd_recv_buf(7,i) * M2g / (scale_d * scale_l**3)
     end do

     com = 0

     nb6_npart = nbufpart
     nb6_tot_mass = sum(nb6_mp(1:nb6_npart))
     write(6,*) 'Total mass of the cluster :', nb6_tot_mass
  end if

  ! Communicating the total mass of the cluster to the other processes
  call MPI_Bcast(nb6_tot_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)

  call timer('rbd_sync', 'start')
end subroutine rbd_sync_cluster

subroutine rbd_sync_gc(add_virtual)
  use rbd_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  
  implicit none

  integer :: ierr, j, i
  logical, intent(in) :: add_virtual

  ! Finding GC owner
  call rbd_compute_gc_owner(add_virtual)

  ! Broadcasting position and velocity
  !write(6,*) 'Broadcasting Guiding center info'
  !call flush(6)

  call MPI_Bcast(rbd_xc,      3, MPI_DOUBLE_PRECISION, rbd_gc_owner-1, MPI_COMM_RAMSES, ierr)
  call MPI_Bcast(rbd_last_xc, 3, MPI_DOUBLE_PRECISION, rbd_gc_owner-1, MPI_COMM_RAMSES, ierr)
  call MPI_Bcast(rbd_vc,      3, MPI_DOUBLE_PRECISION, rbd_gc_owner-1, MPI_COMM_RAMSES, ierr)

end subroutine rbd_sync_gc
