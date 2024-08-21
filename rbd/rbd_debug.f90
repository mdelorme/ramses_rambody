!#################################################################
!# Helpers for Rambody development
!#################################################################

subroutine rbd_dbg_part_count(ilevel, name)
  use rbd_commons
  use amr_commons
  implicit none
  integer, intent(in) :: ilevel
  character(len=*), intent(in) :: name
  
  integer :: sum_part3, igrid, icpu, jgrid
  sum_part3 = 0
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     do jgrid=1,numbl(icpu,ilevel)
        sum_part3 = sum_part3 + rbd_numbp(igrid)
        igrid=next(igrid)
     end do
  end do
  
  write(6,*) 'RBD : ', name, '(', ilevel, ') :', sum_part3

end subroutine rbd_dbg_part_count

subroutine rbd_print_particles
  use rbd_commons
  use amr_commons
  implicit none
  include 'mpif.h'
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale
  logical :: switch
  integer::igrid, jgrid, part_tot, ierr, jpart
  do ilevel=1, nlevelmax
     part_tot = 0
     do icpu=1, ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           part_tot = part_tot + rbd_numbp(igrid)  ! Number of particles in the grid

           npart1 = rbd_numbp(igrid)
           if(npart1>0)then
              ipart=rbd_headp(igrid)
              ! Loop over particles
              do jpart=1,npart1
                 write(6,*) ilevel, ipart, rbd_id(ipart), rbd_xp(1, ipart), rbd_levelp(ipart)
                 
                 ipart=rbd_nextp(ipart)  ! Go to next particle
              end do
              ! End loop over particles
           end if
           
           igrid=next(igrid)
        end do
     end do
  end do
end subroutine rbd_print_particles

subroutine rbd_check_particles(ilevel)
  use rbd_commons
  use amr_commons
  implicit none
  include 'mpif.h'
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale, dx, x, x0
  logical :: switch
  integer::igrid, jgrid, part_tot, ierr, jpart

  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  write(6,*) 'RBD : Checking rbd particles'
  igrid=headl(myid,ilevel)
  ! Loop over grids
  do jgrid=1,numbl(myid,ilevel)
     npart1 = rbd_numbp(igrid)
     if(npart1>0)then
        ipart=rbd_headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           do idim=1,ndim
              x0 = xg(igrid, idim)-3.0D0*dx
              x = rbd_xp(idim, ipart) / scale + skip_loc(idim) - x0
              x = x / dx
              
              if (x < 0.5D0 .or. x > 5.5D0) then
                 write(6,*) 'Error at particle : ', ipart, x, idim
              end if
           end do
           
           ipart=rbd_nextp(ipart)  ! Go to next particle
        end do
        ! End loop over particles
     end if
     
     igrid=next(igrid)
  end do

end subroutine rbd_check_particles

subroutine rbd_ramses_part_per_level(name)
  use amr_commons
  use rbd_commons
  use pm_commons
  implicit none
  include 'mpif.h'
  character(len=*), intent(in) :: name
  integer :: counter, ilevel, igrid, jgrid, tot, ipart, jpart, icpu, part_tot
  write(6,*) '------------'
  tot = 0
  do ilevel=1, nlevelmax
     part_tot = 0
     do icpu=1, ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           part_tot = part_tot + numbp(igrid)
           igrid=next(igrid)
        end do
     end do

     write(6,*) name, ' : Level (', ilevel, ') : ', part_tot
  end do

end subroutine rbd_ramses_part_per_level

subroutine rbd_particles_per_level
  use rbd_commons
  use amr_commons
  implicit none
  include 'mpif.h'
  integer::ipart,idim,i,nxny,ilevel
  integer::npart1,info,icpu,nx_loc
  logical::error
  real(dp),dimension(1:3)::xbound
  integer,dimension(1:nvector),save::ix,iy,iz
  integer,dimension(1:nvector),save::ind_grid,ind_part
  logical,dimension(1:nvector),save::ok=.true.
  real(dp),dimension(1:3)::skip_loc
  real(dp)::scale
  logical :: switch
  integer::igrid, jgrid, part_tot, ierr, jpart
  ! Counting rbd particles on every level of the grid
  rbd_p_per_level = 0
  switch = .false.
  do ilevel=1, nlevelmax
     part_tot = 0
     do icpu=1, ncpu
        igrid=headl(icpu,ilevel)
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           rbd_p_per_level(ilevel) = rbd_p_per_level(ilevel) + rbd_numbp(igrid)
           igrid=next(igrid)
        end do
     end do

     write(6,*) 'rbd_p_per_level (', ilevel, ') : ', rbd_p_per_level(ilevel)
  end do
end subroutine rbd_particles_per_level

subroutine rbd_dbg_gc_owner(tag, ilevel)
  use rbd_commons 
  implicit none

  character(len=*), intent(in) :: tag
  integer, intent(in) :: ilevel

  write(6,*) 'RBD DBG GC :', tag, ilevel, rbd_gc_id
  call flush(6)
  call rbd_sync_gc(.false.)
end subroutine rbd_dbg_gc_owner
