subroutine rbd_init
  use amr_commons
  use rbd_commons
  use pm_commons
  implicit none
  
  include 'mpif.h'
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2

  ! Used for restarting
  integer :: ilun, i, ierr
  character(LEN=5):: nchar, ncharcpu
  character(LEN=80):: fileloc
  logical::ok, is_restart
  integer :: rbd_mesh_Nx2
  real(dp), parameter :: kpc2cm=3.086e21

  allocate(rbd_recv_buf (7,rbdpartmax)) ! 7 variables : 3 pos, 3 vel, 1 mass
  allocate(rbd_send_buf (6,rbdpartmax)) ! 3 variables : fx, fy, fz
  allocate(rbd_xp (ndim, rbdpartmax))
  allocate(rbd_vp (ndim, rbdpartmax))
  allocate(rbd_fp (ndim, rbdpartmax))
  allocate(rbd_mp (rbdpartmax))
  allocate(rbd_id (rbdpartmax))
  allocate(rbd_escapers (7, rbd_escapers_buf))

  ! Nbody6 particles, non mesh, only stars
  allocate(nb6_xp(ndim, rbdpartmax))
  allocate(nb6_vp(ndim, rbdpartmax))
  allocate(nb6_mp(rbdpartmax))
  
  allocate(rbd_nextp  (rbdpartmax))
  allocate(rbd_prevp  (rbdpartmax))
  allocate(rbd_levelp (rbdpartmax))
  
  allocate(rbd_headp(1:ngridmax))
  allocate(rbd_tailp(1:ngridmax))
  allocate(rbd_numbp(1:ngridmax))

  rbd_headp=0
  rbd_tailp=0
  rbd_numbp=0

  rbd_next_sync = -1.0

  ! Setting initial conditions
  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  rbd_xc = rbd_xc0 / scale_l + boxlen * 0.5
  rbd_dbg_buf(:,1) = rbd_xc
  rbd_vc = rbd_vc0 / scale_v
  rbd_fp = 0 ! We set initial force to zero, so we can change the velocity of the cluster guide at t0

  ! Generating Rambody mesh on scale [-0.5, 0,5]
  rbd_mesh_np = rbd_mesh_Nx**3+1
  rbd_mesh_last_scale = 1.0
  rbd_mesh_scale      = 1.0
  rbd_last_xc = rbd_xc
  
  allocate(rbd_mesh_pos(3, rbd_mesh_np))
  ! Mesh force is only required for the original process
  allocate(rbd_mesh_force(3, rbd_mesh_np))
  
  ! Debug
  bar_sync_count = 0

  is_restart = (nrestart > 0)

  if (is_restart .and. rbd_restart) then
     if (myid == 1) then
        call title(nrestart,nchar)

        if(IOGROUPSIZEREP>0)then 
           call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
           fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/rbd_mesh_'//TRIM(nchar)//'.out'
        else
           fileloc='output_'//TRIM(nchar)//'/rbd_mesh_'//TRIM(nchar)//'.out'
        endif
        
        call title(myid,nchar)
        fileloc=TRIM(fileloc)
        inquire(file=fileloc, exist=ok)
        if(.not. ok)then
           write(*,*)'Restart failed:'
           write(*,*)'File '//TRIM(fileloc)//' not found'
           call clean_stop
        end if
        
        ilun=myid+10
        
        open(unit=ilun,file=fileloc,form='unformatted')
        read(ilun) rbd_xc
        read(ilun) rbd_vc
        read(ilun) rbd_gc_id

        write(6,*) 'RBD : Restart :'
        write(6,*) ' - GC ID :', rbd_gc_id
        write(6,*) ' - GC X :', rbd_xc
        write(6,*) ' - GC V :', rbd_vc

        read(ilun) rbd_mesh_np, rbd_mesh_Nx2
        
        if (rbd_mesh_Nx /= rbd_mesh_Nx2) then
           write(*,*) 'Restart failed !'
           write(*,*) 'Mesh incompatible between restarts :'
           write(*,*) 'Current - Restart'
           write(*,*) 'NP : ', rbd_mesh_Nx, rbd_mesh_Nx2
           call clean_stop()
        end if
        
        do i=1, 3
           read(ilun) rbd_mesh_pos(i, 1:rbd_mesh_np)
        end do
        
        do i=1, 3
           read(ilun) rbd_mesh_force(i, 1:rbd_mesh_np)
        end do
        
        call rbd_sync_forces_restart
     end if

     call MPI_Bcast(rbd_xc, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
     call MPI_Bcast(rbd_vc, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
     call MPI_Bcast(rbd_gc_id, 1, MPI_INTEGER, 0, MPI_COMM_RAMSES, ierr)
  else
     call rbd_init_gc(is_restart)
  end if

  ! Sending epsilon
  if (myid == 1) then
     if (rbd_dbg_comm) then
        write(6,*) 'Comm tag 18 (-> smoothing length)'
        call flush(6)
     end if
     call MPI_Send(rbd_epsilon, 1, MPI_DOUBLE_PRECISION, 0, 18, MPI_COMM_RAMBODY, ierr)
  end if
  next_pid = rbd_gc_id + 1

end subroutine rbd_init

subroutine rbd_init_gc(is_restart)
  use rbd_commons
  use amr_commons
  use pm_commons
  implicit none

  include 'mpif.h'

  logical, intent(in) :: is_restart
  integer :: i, j, k, pid, ierr, nx_loc, nxny, ix, iy, iz
  real(dp), dimension(1:1, 1:3) :: xx_dp
  integer, dimension(1:1) :: ind_grid, ind_part
  logical, dimension(1:1) :: ok=.true.
  integer, dimension(1:1) :: cc
  real(dp) :: scale
  real(dp),dimension(1:3) :: skip_loc
  integer :: ip, ilevel, icpu, gc_id_all
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all

  
  nxny = nx*ny
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  xx_dp(1, 1) = rbd_xc(1) / boxlen
  xx_dp(1, 2) = rbd_xc(2) / boxlen
  xx_dp(1, 3) = rbd_xc(3) / boxlen

  !call cmp_cpumap(xx_dp, cc, 1)
  
  ! We take the maximum particle id
  ! WARNING : This one does not work if there are some unloaded particles in the IC files !
  ! then the maximum idp does not correspond to npartmax
  !call MPI_Allreduce(npart, ip, 1, MPI_INT, MPI_SUM, MPI_COMM_RAMSES, ierr)
  
  ip = 0
  do i=1, npart
     ip = MAX(ip, idp(i))
  end do
  call MPI_Reduce(ip, rbd_gc_id, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_RAMSES, ierr)

  ! Doesn't seem to work for restarts as the grid is already refined/load balanced
  !if (cc(1) == myid) then
  ! Instead we use whichever process has the level 1 grid
  if (numbl(myid, 1) == 1) then
     write(6,*) 'Current process', myid, ' has ', numbl(myid, 1), 'grids on level 1'
     rbd_gc_id = rbd_gc_id+1
     npart = npart + 1
     
     idp(npart) = rbd_gc_id
     gc_id_all = rbd_gc_id
     
     xp(npart, :) = rbd_xc
     vp(npart, :) = rbd_vc
     mp(npart) = 0.0
     typep(npart)%family = FAM_TRACER_STAR
     typep(npart)%tag = 0
     levelp(npart) = levelmin
     write(6,*) 'RBDDBG ! Here on process ', myid, 'GCID = ', rbd_gc_id
  else
     gc_id_all = 0
  endif

  npart_cpu=0; npart_all=0
  npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_RAMSES,ierr)
#else
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_RAMSES,ierr)
#endif
  npart_cpu(1)=npart_all(1)
#endif
  do icpu=2,ncpu
     npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
  end do
  
  call MPI_AllReduce(gc_id_all, rbd_gc_id, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_RAMSES, ierr)


end subroutine rbd_init_gc

