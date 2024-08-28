!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_init_tree
  use rbd_commons
  use rbd_parameters
  use amr_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------
  ! This subroutine build the particle linked list at the 
  ! coarse level for ALL the particles in the box.
  ! This routine should be used only as initial set up for
  ! the particle tree.
  !------------------------------------------------------
  integer::ipart,idim,i,nxny,ilevel,j
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale
  integer::igrid, jgrid, part_tot, ierr, N_part

  if(verbose)write(*,*)' RBD : Initializing tree'
  rbd_total_added = 0
  rbd_total_removed = 0

  ! Local constants
  nxny=nx*ny
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  !----------------------------------
  ! Initialize particle linked list
  !----------------------------------
  ! Hard reset of everything
  rbd_prevp = 0
  rbd_nextp = 0
  
  rbd_numbp = 0
  rbd_headp = 0
  rbd_tailp = 0
  
  rbd_prevp(1)=0; rbd_nextp(1)=2
  do ipart=2,rbdpartmax-1
     rbd_prevp(ipart)=ipart-1
     rbd_nextp(ipart)=ipart+1
  end do
  
  rbd_prevp(rbdpartmax)=rbdpartmax-1;
  rbd_nextp(rbdpartmax)=0

  ! Free memory linked list
  rbd_headp_free=rbd_npart+1
  rbd_tailp_free=rbdpartmax
  rbd_numbp_free=rbd_tailp_free-rbd_headp_free+1

  if(rbd_numbp_free>0)then
     rbd_prevp(rbd_headp_free)=0
  end if
  rbd_nextp(rbd_tailp_free)=0
  call MPI_ALLREDUCE(rbd_numbp_free,rbd_numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_RAMSES,info)

  !--------------
  ! Periodic box
  !--------------

  ! Maybe not relevant here since particles are reinitialized
  do idim=1,ndim
     do j=1,rbd_npart
        ipart = j
        if(rbd_xp(idim, ipart)/scale+skip_loc(idim)<0.0d0) &
             & rbd_xp(idim, ipart)=rbd_xp(idim, ipart)+(xbound(idim)-skip_loc(idim))*scale
        if(rbd_xp(idim, ipart)/scale+skip_loc(idim)>=xbound(idim)) &
             & rbd_xp(idim, ipart)=rbd_xp(idim, ipart)-(xbound(idim)-skip_loc(idim))*scale
     end do
  end do
 
  !----------------------------------
  ! Reset all linked lists at level 1
  !----------------------------------
  do i=1,active(1)%ngrid
     rbd_headp(active(1)%igrid(i))=0
     rbd_tailp(active(1)%igrid(i))=0
     rbd_numbp(active(1)%igrid(i))=0
  end do
  do icpu=1,ncpu
     do i=1,reception(icpu,1)%ngrid
        rbd_headp(reception(icpu,1)%igrid(i))=0
        rbd_tailp(reception(icpu,1)%igrid(i))=0
        rbd_numbp(reception(icpu,1)%igrid(i))=0
     end do
  end do

  !------------------------------------------------
  ! Build linked list at level 1 by vector sweeps
  !------------------------------------------------
  N_part = rbd_npart
  rbd_npart = 0
  do ipart=1,N_part,nvector
     npart1=min(nvector,N_part-ipart+1)
     ! Gather particles
     do i=1,npart1
        ! Previous ind_part(i) = ipart+i-1
        ind_part(i) = ipart+i-1
        !ind_part(i)=rbd_id(ipart+i-1)
     end do
     ! Compute coarse cell
#if NDIM>0
     do i=1,npart1
        ix(i)=int(rbd_xp(1, ind_part(i))/scale+skip_loc(1))
     end do
#endif
#if NDIM>1
     do i=1,npart1
        iy(i)=int(rbd_xp(2, ind_part(i))/scale+skip_loc(2))
     end do
#endif
#if NDIM>2
     do i=1,npart1
        iz(i)=int(rbd_xp(3, ind_part(i))/scale+skip_loc(3))
     end do
#endif
     ! Compute level 1 grid index
     error=.false.
     do i=1,npart1
        ind_grid(i)=son(1+ix(i)+nx*iy(i)+nxny*iz(i))
        if(ind_grid(i)==0)error=.true.
     end do
     if(error)then 
        write(*,*)'Error in rbd_init_tree'
        write(*,*)'Particles appear in unrefined regions'
        call clean_stop
     end if
     ! Add particle to level 1 linked list
     call rbd_add_list(ind_part,ind_grid,ok,npart1)
  end do

  ! Sort particles down to levelmin
  do ilevel=1, nlevelmax
     call rbd_make_tree_fine(ilevel)
     call rbd_virtual_tree_fine(ilevel)
     call rbd_kill_tree_fine(ilevel)
  end do

  call rbd_compute_level_min

end subroutine rbd_init_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_make_tree_fine(ilevel)
  use rbd_commons
  use amr_commons
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine checks if particles have moved from their parent grid
  ! to one of the 3**ndim neighboring sister grids. The particle is then 
  ! disconnected from the parent grid linked list, and connected to the
  ! corresponding sister grid linked list. If the sister grid does
  ! not exist, the particle is left to its original parent grid.
  ! Particles must not move to a distance greater than direct neighbors
  ! boundaries. Otherwise an error message is issued and the code stops.
  !-----------------------------------------------------------------------
  integer::idim,nx_loc
  real(dp)::dx,scale
  real(dp),dimension(1:3)::xbound
  real(dp),dimension(1:3)::skip_loc
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=rbd_numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=rbd_headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <--- Very important !!!
              next_part=rbd_nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              ! Gather nvector particles
              if(ip==nvector)then
                 call rbd_check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call rbd_check_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

111 format('   Entering rbd_make_tree_fine for level ',I2)

end subroutine rbd_make_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_check_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use rbd_commons
  use pm_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by make_tree_fine.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  real(dp),dimension(1:3)::xbound
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_father
  ! Particle-based arrays
  integer,dimension(1:nvector),save::ind_son,igrid_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel    
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_father(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_father,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Compute particle position in 3-cube
  error=.false.
  ind_son(1:np)=1
  ok(1:np)=.false.
  do idim=1,ndim
     do j=1,np
        i=int((rbd_xp(idim, ind_part(j))/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
        if(i<0.or.i>2)error=.true.
        ind_son(j)=ind_son(j)+i*3**(idim-1)
        ! Check if particle has escaped from its parent grid
        ok(j)=ok(j).or.i.ne.1 
     end do
  end do
  if(error)then
     write(*,*)'Problem in rbd_check_tree at level ',ilevel
     write(*,*)'A particle has moved outside allowed boundaries'
     do idim=1,ndim
        do j=1,np
           i=int((rbd_xp(idim, ind_part(j))/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx/2.0D0)
           if(i<0.or.i>2)then
              write(*,*) rbd_xp(idim, ind_part(j)),x0(ind_grid_part(j),idim), i, ind_part(j)
           endif
        end do
     end do
     stop
  end if

  ! Compute neighboring grid index
  do j=1,np
     igrid_son(j)=son(nbors_father_cells(ind_grid_part(j),ind_son(j)))
  end do

  ! If escaped particle sits in unrefined cell, leave it to its parent grid.
  ! For ilevel=levelmin, this should never happen.
  do j=1,np
     if(igrid_son(j)==0)ok(j)=.false.
  end do

  ! Periodic box
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xxx=rbd_xp(idim, ind_part(j))/scale+skip_loc(idim)-xg(igrid_son(j),idim)
           if(xxx> xbound(idim)/2.0)then
              rbd_xp(idim, ind_part(j))=xp(ind_part(j),idim)-(xbound(idim)-skip_loc(idim))*scale
           endif
           if(xxx<-xbound(idim)/2.0)then
              rbd_xp(idim, ind_part(j))=xp(ind_part(j),idim)+(xbound(idim)-skip_loc(idim))*scale
           endif
        endif
     enddo
  enddo

  ! Switch particles linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do
  call rbd_remove_list(ind_part,list1,ok,np)
  call rbd_add_list(ind_part,list2,ok,np)

end subroutine rbd_check_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_kill_tree_fine(ilevel)
  use rbd_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine sorts particle between ilevel grids and their 
  ! ilevel+1 children grids. Particles are disconnected from their parent 
  ! grid linked list and connected to their corresponding child grid linked 
  ! list. If the  child grid does not exist, the particle is left to its 
  ! original parent grid. 
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::i,ig,ip,npart1,icpu
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer :: tot, part_tot, tmp_ilevel

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel+1)==0)return
  if(verbose)write(*,111)ilevel

  ! Reset all linked lists at level ilevel+1
  do i=1,active(ilevel+1)%ngrid
     rbd_headp(active(ilevel+1)%igrid(i))=0
     rbd_tailp(active(ilevel+1)%igrid(i))=0
     rbd_numbp(active(ilevel+1)%igrid(i))=0
  end do
    
  do icpu=1,ncpu
     do i=1,reception(icpu,ilevel+1)%ngrid
        rbd_headp(reception(icpu,ilevel+1)%igrid(i))=0
        rbd_tailp(reception(icpu,ilevel+1)%igrid(i))=0
        rbd_numbp(reception(icpu,ilevel+1)%igrid(i))=0
     end do
  end do
  
  ! Sort particles between ilevel and ilevel+1
  tot = 0
  part_tot = 0

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=rbd_numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then        
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=rbd_headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=rbd_nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   

              rbd_levelp(ind_part(ip)) = ilevel

              if(ip==nvector)then
                 call rbd_kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,tot)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0)call rbd_kill_tree(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,tot)
  end do 
  ! End loop over cpus

111 format('   Entering rbd_kill_tree_fine for level ',I2)

end subroutine rbd_kill_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_kill_tree(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,tot)
  use amr_commons
  use rbd_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine rbd_kill_tree_fine.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,xxx,scale
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle based arrays
  integer,dimension(1:nvector),save::igrid_son,ind_son
  integer,dimension(1:nvector),save::list1,list2
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:3)::skip_loc

  integer :: tot

  ! Mesh spacing in that level
  dx=0.5D0**ilevel   
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Compute lower left corner of grid
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-dx
     end do
  end do

  ! Select only particles within grid boundaries
  ok(1:np)=.true.
  do idim=1,ndim
     do j=1,np
        xxx=(rbd_xp(idim, ind_part(j))/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx
        ok(j)=ok(j) .and. (xxx >= 0.d0 .and. xxx < 2.0d0)
     end do
  end do

  ! Debug
  do j=1,np
     if (ok(j)) tot = tot + 1
  end do
  
  ! Determines in which son particles sit
  ind_son(1:np)=0
  do idim=1,ndim
     do j=1,np
        i=int((rbd_xp(idim, ind_part(j))/scale+skip_loc(idim)-x0(ind_grid_part(j),idim))/dx)
        ind_son(j)=ind_son(j)+i*2**(idim-1)
     end do
  end do 
  do j=1,np
     ind_son(j)=ncoarse+ind_son(j)*ngridmax+ind_grid(ind_grid_part(j))
  end do

  ! Determine which son cell is refined
  igrid_son(1:np)=0
  do j=1,np
     if(ok(j))igrid_son(j)=son(ind_son(j))
  end do
  do j=1,np
     ok(j)=igrid_son(j)>0
  end do

  ! Compute particle linked list
  do j=1,np
     if(ok(j))then
        list1(j)=ind_grid(ind_grid_part(j))
        list2(j)=igrid_son(j)
     end if
  end do
  
  ! Remove particles from their original linked lists
  call rbd_remove_list(ind_part,list1,ok,np)
  ! Add particles to their new linked lists
  call rbd_add_list(ind_part,list2,ok,np)

end subroutine rbd_kill_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_merge_tree_fine(ilevel)
  use rbd_commons
  use amr_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------
  ! This routine disconnects all particles contained in children grids
  ! and connects them to their parent grid linked list.
  !---------------------------------------------------------------
  integer::igrid,iskip,icpu
  integer::i,ind,ncache,ngrid
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_grid_son
  logical,dimension(1:nvector),save::ok

  if(numbtot(1,ilevel)==0)return
  if(ilevel==nlevelmax)return
  if(verbose)write(*,111)ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     if(icpu==myid)then
        ncache=active(ilevel)%ngrid
     else
        ncache=reception(icpu,ilevel)%ngrid
     end if
     ! Loop over grids by vector sweeps
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        if(icpu==myid)then
           do i=1,ngrid
              ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
           end do
        else
           do i=1,ngrid
              ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
           end do
        end if
        ! Loop over children grids
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ngrid
              ind_grid_son(i)=son(ind_cell(i))
           end do
           do i=1,ngrid
              ok(i)=ind_grid_son(i)>0
           end do
           do i=1,ngrid
           if(ok(i))then
              if(rbd_numbp(ind_grid_son(i))>0)then
                 if(rbd_numbp(ind_grid(i))>0)then
                    ! Connect son linked list at the tail of father linked list
                    rbd_nextp(rbd_tailp(ind_grid(i)))=rbd_headp(ind_grid_son(i))
                    rbd_prevp(rbd_headp(ind_grid_son(i)))=rbd_tailp(ind_grid(i))
                    rbd_numbp(ind_grid(i))=rbd_numbp(ind_grid(i))+rbd_numbp(ind_grid_son(i))
                    rbd_tailp(ind_grid(i))=rbd_tailp(ind_grid_son(i))
                 else
                    ! Initialize father linked list
                    rbd_headp(ind_grid(i))=rbd_headp(ind_grid_son(i))
                    rbd_tailp(ind_grid(i))=rbd_tailp(ind_grid_son(i))
                    rbd_numbp(ind_grid(i))=rbd_numbp(ind_grid_son(i))
                 end if
                 
              end if
           end if
           end do
        end do
        ! End loop over children
     end do
     ! End loop over grids
  end do
  ! End loop over cpus

111 format('   Entering rbd_merge_tree_fine for level ',I2)

end subroutine rbd_merge_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_virtual_tree_fine(ilevel)
  use rbd_commons
  use amr_commons
  use mpi_mod
  
  implicit none
  integer::ilevel
  !-----------------------------------------------------------------------
  ! This subroutine move particles across processors boundaries.
  !-----------------------------------------------------------------------
  integer::igrid,ipart,jpart,ncache_tot,next_part
  integer::ip,ipcom,npart1,icpu,ncache
  integer::info,buf_count,tag=101,tagf=102,tagu=102
  integer::countsend,countrecv
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
  integer,dimension(1:nvector),save::ind_part,ind_list,ind_com
  logical::ok_free,ok_all
  integer::particle_data_width

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Count particle sitting in virtual boundaries
  do icpu=1,ncpu
     reception(icpu,ilevel)%npart=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        
        reception(icpu,ilevel)%npart = reception(icpu,ilevel)%npart + rbd_numbp(reception(icpu,ilevel)%igrid(igrid))
     end do
     sendbuf(icpu)=reception(icpu,ilevel)%npart
  end do

  ! Calculate how many particle properties are being transferred
  particle_data_width = 7 ! 1 mass, 3 pos, 3 vel

  ! Allocate communication buffer in emission
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(reception(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(reception(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Gather particle in communication buffer
  do icpu=1,ncpu
     if(reception(icpu,ilevel)%npart>0)then
        ! Gather particles by vector sweeps
        ipcom=0
        ip=0
        do igrid=1,reception(icpu,ilevel)%ngrid
           npart1=rbd_numbp(reception(icpu,ilevel)%igrid(igrid))
           ipart =rbd_headp(reception(icpu,ilevel)%igrid(igrid))
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <--- Very important !!!
              next_part=rbd_nextp(ipart)
              ip=ip+1
              ipcom=ipcom+1
              ind_com (ip)=ipcom
              ind_part(ip)=ipart
              ind_list(ip)=reception(icpu,ilevel)%igrid(igrid)
              reception(icpu,ilevel)%fp(ipcom,1)=igrid
              if(ip==nvector)then
                 call rbd_fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
                 ip=0
              end if
              ipart=next_part  ! Go to next particle
           end do
        end do
        if(ip>0)call rbd_fill_comm(ind_part,ind_com,ind_list,ip,ilevel,icpu)
     end if
  end do
  
  ! Communicate virtual particle number to parent cpu
  call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_RAMSES,info)

  ! Allocate communication buffer in reception
  do icpu=1,ncpu
     emission(icpu,ilevel)%npart=recvbuf(icpu)
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        ! Allocate reception buffer
        allocate(emission(icpu,ilevel)%fp(1:ncache,1:3))
        allocate(emission(icpu,ilevel)%up(1:ncache,1:particle_data_width))
     end if
  end do

  ! Receive particles
  countrecv=0
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countrecv=countrecv+1
#ifndef LONGINT
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_RAMSES,reqrecv(countrecv),info)
#else
        call MPI_IRECV(emission(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_RAMSES,reqrecv(countrecv),info)
#endif
        buf_count=ncache*particle_data_width
        countrecv=countrecv+1
        call MPI_IRECV(emission(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_RAMSES,reqrecv(countrecv),info)
     end if
  end do

  ! Send particles
  countsend=0
  do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        buf_count=ncache*3
        countsend=countsend+1
#ifndef LONGINT
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER,icpu-1,&
             & tagf,MPI_COMM_RAMSES,reqsend(countsend),info)
#else
        call MPI_ISEND(reception(icpu,ilevel)%fp,buf_count, &
             & MPI_INTEGER8,icpu-1,&
             & tagf,MPI_COMM_RAMSES,reqsend(countsend),info)
#endif
        buf_count=ncache*particle_data_width
        countsend=countsend+1
        call MPI_ISEND(reception(icpu,ilevel)%up,buf_count, &
             & MPI_DOUBLE_PRECISION,icpu-1,&
             & tagu,MPI_COMM_RAMSES,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Compute total number of newly created particles
  ncache_tot=0
  do icpu=1,ncpu
     ncache_tot=ncache_tot+emission(icpu,ilevel)%npart
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  call MPI_ALLREDUCE(rbd_numbp_free, rbd_numbp_free_tot,1,MPI_INTEGER,MPI_MIN,&
       & MPI_COMM_RAMSES,info)
  ok_free=(rbd_numbp_free-ncache_tot)>=0
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase rbdpartmax'
     write(*,*)rbd_numbp_free,ncache_tot
     write(*,*)myid
     write(*,*)emission(1:ncpu,ilevel)%npart
     write(*,*)'============================'
     write(*,*)reception(1:ncpu,ilevel)%npart
     call MPI_ABORT(MPI_COMM_WORLD,1,info)
  end if

  ! Scatter new particles from communication buffer
  do icpu=1,ncpu
     ! Loop over particles by vector sweeps
     ncache=emission(icpu,ilevel)%npart
     do ipart=1,ncache,nvector
        npart1=min(nvector,ncache-ipart+1)
        do ip=1,npart1
           ind_com(ip)=ipart+ip-1
        end do
        call rbd_empty_comm(ind_com,npart1,ilevel,icpu)
     end do
  end do

  ! Deallocate temporary communication buffers
  do icpu=1,ncpu
     ncache=emission(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(emission(icpu,ilevel)%fp)
        deallocate(emission(icpu,ilevel)%up)
     end if
     ncache=reception(icpu,ilevel)%npart
     if(ncache>0)then
        deallocate(reception(icpu,ilevel)%fp)
        deallocate(reception(icpu,ilevel)%up)
     end if
  end do

111 format('   Entering virtual_tree_fine for level ',I2)
end subroutine rbd_virtual_tree_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_fill_comm(ind_part,ind_com,ind_list,np,ilevel,icpu)
  use rbd_commons
  use amr_commons
  implicit none
  integer::np,ilevel,icpu
  integer,dimension(1:nvector)::ind_part,ind_com,ind_list
  integer::current_property
  integer::i,idim
  logical,dimension(1:nvector),save::ok=.true.

  ! Gather particle level and identity
  do i=1,np
     reception(icpu,ilevel)%fp(ind_com(i),2)=rbd_levelp(ind_part(i))
     reception(icpu,ilevel)%fp(ind_com(i),3)=rbd_id(ind_part(i))
  end do
  
  ! Gather particle position, velocity and mass
  do idim=1,ndim
     do i=1,np
        reception(icpu,ilevel)%up(ind_com(i),idim     )=rbd_xp(idim, ind_part(i))
        reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=rbd_vp(idim, ind_part(i))
     end do
  end do

  do i=1, np
     reception(icpu,ilevel)%up(ind_com(i),7) = rbd_mp(ind_part(i))
  end do
  
  ! Remove particles from parent linked list
  call rbd_remove_list(ind_part,ind_list,ok,np)
  call rbd_add_free(ind_part,np)
  
end subroutine rbd_fill_comm
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_empty_comm(ind_com,np,ilevel,icpu)
  use rbd_commons
  use amr_commons
  implicit none
  include  'mpif.h'
  integer::np,icpu,ilevel
  integer,dimension(1:nvector)::ind_com
  
  integer::i,idim,igrid
  integer,dimension(1:nvector),save::ind_list,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  integer::current_property

  ! Compute parent grid index
  do i=1,np
     igrid=emission(icpu,ilevel)%fp(ind_com(i),1)
     ind_list(i)=emission(icpu,ilevel)%igrid(igrid)
  end do

  ! Add particle to parent linked list
  call rbd_remove_free(ind_part,np)
  call rbd_add_list(ind_part,ind_list,ok,np)

  ! Scatter particle level and identity
  do i=1,np
     rbd_levelp(ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),2)
     rbd_id    (ind_part(i))=emission(icpu,ilevel)%fp(ind_com(i),3)
     write(6,*) 'Rbd id(', ind_part(i), ') = ', rbd_id(ind_part(i))
     call flush(6)
  end do

  ! Scatter particle position and velocity
  do idim=1,ndim
     do i=1,np
        rbd_xp(idim, ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),idim     )
        rbd_vp(idim, ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),idim+ndim)
     end do
  end do

  current_property = twondim+1

  ! Scatter particle mass
  do i=1,np
     rbd_mp(ind_part(i))=emission(icpu,ilevel)%up(ind_com(i),current_property)
  end do

end subroutine rbd_empty_comm
!################################################################
!################################################################
!################################################################
!################################################################

subroutine rbd_compute_level_min
  use rbd_commons
  use amr_commons
  use mpi_mod
  implicit none
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale
  integer::igrid, jgrid, part_tot, ierr, jpart, tmp
  ! Counting rbd particles on every level of the grid
  rbd_level_min = nlevelmax
  rbd_p_per_level = 0
  do ilevel=nlevelmax, levelmin, -1
     part_tot = 0
     do icpu=1, ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           part_tot = part_tot + rbd_numbp(igrid)  ! Number of particles in the grid
           rbd_p_per_level(ilevel) = rbd_p_per_level(ilevel) + rbd_numbp(igrid)
           igrid=next(igrid)
        end do
     end do
     
     if (part_tot > 0 .and. ilevel < rbd_level_min) then
        rbd_level_min = ilevel
     end if
  end do

  tmp = rbd_level_min
  call MPI_Allreduce(tmp, rbd_level_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_RAMSES, ierr)
  ! Check, should be useless but better safe than sorry
  rbd_level_min = max(rbd_level_min, levelmin)

end subroutine rbd_compute_level_min


subroutine rbd_compute_gc_owner(add_virtual)
  use rbd_commons
  use amr_commons
  use pm_commons
  use mpi_mod
  
  implicit none
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale
  logical, intent(in) :: add_virtual
  integer::igrid, jgrid, part_tot, ierr, jpart, tmp
  ! Counting rbd particles on every level of the grid
  !rbd_level_min = nlevelmax
  rbd_p_per_level = 0
  rbd_gc_owner = -1
  rbd_gc_level = -1

  do ilevel=1, nlevelmax
     igrid = headl(myid, ilevel)

     do jgrid=1, numbl(myid, ilevel)
        npart1 = numbp(igrid)

        if (npart1 > 0) then
           ipart = headp(igrid)

           do jpart=1, npart1
              if (idp(ipart) == rbd_gc_id) then
                 rbd_gc_owner = myid
                 rbd_gc_level = ilevel
              end if

              ipart = nextp(ipart)
           end do
           if (rbd_gc_owner > -1) exit
        end if

        igrid = next(igrid)
     end do

     ! And also in the virtual boundaries
     if (add_virtual .and. rbd_gc_owner == -1) then
        do icpu=1, ncpu
           do jgrid=1, reception(icpu,ilevel)%ngrid
              igrid = reception(icpu, ilevel)%igrid(jgrid)
              
              npart1 = numbp(igrid)
              ipart = headp(igrid)
              
              do jpart=1, npart1
                 if (idp(ipart) == rbd_gc_id) then
                    rbd_gc_owner = myid
                    rbd_gc_level = ilevel
                    exit
                 end if
                 
                 ipart = nextp(ipart)
              end do
              if (rbd_gc_owner > -1) exit
           end do
           if (rbd_gc_owner > -1) exit
        end do
     end if

     tmp = rbd_gc_owner
     call MPI_Allreduce(tmp, rbd_gc_owner, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_RAMSES, ierr)

     if (rbd_gc_owner > -1) exit
  end do
     
  tmp = rbd_gc_level
  call MPI_Allreduce(tmp, rbd_gc_level, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_RAMSES, ierr)
  write(6,*) 'RBD : GC owner = process #', rbd_gc_owner

  if (myid == 1 .and. rbd_gc_owner == -1) then
     write(6,*) 'ERROR : Rambody guiding center has disappeared !', add_virtual
     call flush(6)
     call MPI_ABORT(MPI_COMM_WORLD, 111, ierr)
  end if

end subroutine rbd_compute_gc_owner
