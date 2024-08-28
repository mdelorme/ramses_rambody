subroutine rbd_sync_escapers
  use amr_commons
  use pm_commons
  use rbd_commons
  use mpi_mod
  
  implicit none

  integer  :: i, j, ierr
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp), parameter :: pc2cm   = 3.086e18
  real(dp), parameter :: kms2cms = 1.0e5
  real(dp), parameter :: M2g     = 1.988e33
  
  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  call timer('rbd_wait_escapers', 'start')
  if (myid == 1) then
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 4 (<- N_escapers)'
        call flush(6)
     end if

     call MPI_Recv(rbd_n_escapers, 1, MPI_INTEGER, 0, 4, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)
  end if

  call MPI_Bcast(rbd_n_escapers, 1, MPI_INTEGER, 0, MPI_COMM_RAMSES, ierr)
  write(6,*) 'RBD : Broadcasted/Received escaper count : ', rbd_n_escapers

  rbd_tot_escapers = rbd_tot_escapers + rbd_n_escapers

  ! TODO : Check memory is enough to hold the new particles ?

  if (rbd_n_escapers > 0 .and. rbd_store_escapers) then
     ! Receiving particles and dispatching
     if (myid == 1) then
        if (rbd_dbg_comm) then
           write(6,*) 'RBD : Comm tag 5 (<- Escapers)'
           call flush(6)
        end if
        
        call MPI_Recv(rbd_recv_buf, rbd_n_escapers*7, MPI_DOUBLE_PRECISION, 0, 5, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)
      
        do i=1, rbd_n_escapers
           do j=1,3
              rbd_escapers(j,i) = rbd_recv_buf(j,i) * pc2cm / scale_l
              rbd_escapers(j+3,i) = rbd_recv_buf(j+3,i) * kms2cms / scale_v
           end do
           ! Do we really care about this one ?
           rbd_escapers(7,i) = rbd_recv_buf(7,i) * M2g / (scale_d * scale_l**3.0)
        end do
     end if

     call MPI_Bcast(rbd_escapers, rbd_n_escapers*7, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
  end if
  call timer('rbd_sync', 'start')
end subroutine rbd_sync_escapers

subroutine rbd_push_escapers
  use amr_commons
  use pm_commons
  use rbd_commons
  use mpi_mod

  implicit none

  logical, allocatable, dimension(:)  :: pushed
  integer, dimension(1:nvector), save :: ind_grid, ind_cell, ind_part
  logical, dimension(1:nvector), save :: ok, ok_new=.true.
  integer, allocatable, dimension(:)  :: ind_grid_new, ind_part_new, ind_level_part
  real(dp),dimension(1:3) :: skip_loc, cellx, escx
  integer  :: nx_loc, ix, iy, iz, ip, op, idp_new, counter
  real(dp) :: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v, scale
  real(dp) :: dx_loc, dx_min, dx
  real(dp),dimension(1:twotondim,1:3)::xc
  integer  :: igrid, jgrid, icpu, ncache, ngrid, iskip, ind, i, j, ilevel, ipart, nnew, npart_loc, nnew_tot, ierr
  logical  :: valid

  nnew_tot = 0
  
  if (rbd_n_escapers > 0) then
     ! Conversion factor from user units to cgs units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     
     allocate(pushed(rbd_n_escapers))
     allocate(ind_grid_new(rbd_n_escapers))
     allocate(ind_part_new(rbd_n_escapers))
     allocate(ind_level_part(rbd_n_escapers))
     pushed = .false.

     ! Finding what id we should start at
     counter = 0
     ! Debug : Counting particle
     do ilevel=levelmin, nlevelmax
        ncache = active(ilevel)%ngrid
        
        do jgrid=1, ncache
           igrid = active(ilevel)%igrid(jgrid)
           if (numbp(igrid) > 0) then
              counter = counter + numbp(igrid)
           end if
        end do
     end do

     ! Finding where the particles are on the grid
     nnew = 0
     ilevel = levelmin
     dx=0.5D0**(ilevel-1)
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     skip_loc(1)=dble(icoarse_min)
     skip_loc(2)=dble(jcoarse_min)
     skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     dx_min=(0.5D0**nlevelmax)*scale
     
     ncache = active(ilevel)%ngrid
     
     ! Vectorization
     do igrid=1, ncache, nvector
        ngrid=MIN(nvector, ncache-igrid+1)
        
        ! Looping over cells
        do i=1, ngrid
           ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
           ! Calculate the position of the cell
           do j=1,3
              cellx(j) = xg(ind_grid(i),j) * scale
           end do
           
           ! We loop through all particles and flag the ones in the cell
           do ipart=1, rbd_n_escapers
              ! We skip particles already on the grid
              if (pushed(ipart)) cycle
              
              valid = .true.
              do j=1,3
                 escx(j) = rbd_escapers(j,ipart) + rbd_xc(j)
                 if (escx(j) < cellx(j) - 0.5*dx_loc .or. escx(j) > cellx(j) + 0.5*dx_loc) then
                    valid = .false.
                    exit
                 end if
              end do
              
              if (valid) then
                 nnew = nnew + 1
                 ind_grid_new(nnew) = ind_grid(i)
                 ind_part_new(nnew) = ipart
                 ind_level_part(nnew) = ilevel
                 !write(6,*) dbg_counter, ' PUSHING ! ', ind_level_part(ipart), ilevel, nnew
                 pushed(ipart) = .true.
              end if
           end do
        end do
     end do

     ! Now we push everything
     do ipart=1, nnew, nvector
        npart_loc = min(nvector, nnew-ipart)
        call remove_free(ind_part, npart_loc)
        call add_list(ind_part, ind_grid_new(ipart:ipart+npart_loc), ok_new, npart_loc)

        do i=1, npart_loc
           ip = ind_part(i)              ! Index on this process
           op = ind_part_new(ipart+i-1)  ! Index in the escapers table

           if (ip > npartmax) then
              write(6,*) 'Not enough memory, increase npartmax'
              call clean_stop
           end if

           ! Setting up position
           do j=1, 3
              xp(ip, j) = rbd_escapers(j, op)   + rbd_xc(j)
              vp(ip, j) = rbd_escapers(j+3, op) + rbd_vc(j)
           end do
           mp(ip) = rbd_escapers(7, op)

           ! TODO : Set this accordingly to NBody6
           if (star) then
              tp(ip) = 0.0
              if (metal) then
                 zp(ip) = 0.0
              end if
           end if

           idp(ip) = next_pid + ipart + i - 1
           typep(ip)%family = FAM_TRACER_STAR ! TODO : Change that to STAR
           typep(ip)%tag = 1 ! Setting a specific tag to the escapers so we can find them easily

           !rbd_dbg_esc(:,ipart) = xp(ip,:)
           
           levelp(ip) = ind_level_part(ipart+i-1)
        end do
     end do

     next_pid = next_pid + rbd_n_escapers

     deallocate(pushed)
     deallocate(ind_grid_new)
     deallocate(ind_part_new)
     deallocate(ind_level_part)

     write(6,*) 'RBD : Pushed ', nnew, ' escapers locally'
     call MPI_Reduce(nnew, nnew_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_RAMSES, ierr)
     if (myid == 1) then
        write(6,*) 'RBD : Pushed ', nnew_tot, ' escapers in total (', rbd_n_escapers, ')'
        if (nnew_tot /= rbd_n_escapers) then
           write(6,*) 'ERROR : Number of pushed particles is different than number of escapers !'
           call flush(6)
           call MPI_ABORT(MPI_COMM_WORLD, 165, ierr)
        end if
     end if

     ! Debug : Counting particle
     rbd_n_escapers = 0
  end if  
end subroutine rbd_push_escapers
