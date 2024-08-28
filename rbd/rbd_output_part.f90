subroutine rbd_output_part(filename)

  use amr_commons
  use rbd_commons
  use pm_commons
  use mpi_mod  
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ilun,ipart, idpart, pid
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  integer(i8b),allocatable,dimension(:)::ii8
  integer,allocatable,dimension(:)::ll
  logical,allocatable,dimension(:)::nb
  integer, allocatable, dimension(:)::ind_part
  integer,parameter::tag=1122
  integer::dummy_io,info2, ierr, tpart

  integer:: ilevel, icpu, npart1, igrid, jgrid, jpart, cur_part, tot_part, tot_esc

  real(dp), dimension(1:3) :: dbg

  if(verbose)write(*,*)'  Entering rbd_output_part'

  !call rbd_sync_gc(.false.)

  call MPI_Reduce(npart, tpart, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_RAMSES, ierr)
  
  if (myid == 1) then
     ilun=2*ncpu+myid+10
     fileloc=TRIM(filename)
     open(unit=ilun,file=TRIM(fileloc),form='unformatted')
     rewind(ilun)
     ! Write header
     write(ilun)ncpu
     write(ilun)ndim
     write(ilun)nb6_npart
     write(ilun)rbd_mesh_last_scale

     write(6,*) 'RBD : Writing ', nb6_npart, ' particles to disk'
     
     ! Cluster guiding center
     write(ilun) rbd_xc
     write(ilun) rbd_vc
     write(ilun) rbd_gc_owner
     allocate(xdp(nb6_npart))

     !if (rbd_dbg_nesc > 0) then
     !   dbg = 0
        ! Debug escapers
        !do i=1, rbd_dbg_nesc
        !   dbg = dbg + rbd_dbg_esc(:, i)
        !end do
        
       !dbg = dbg / rbd_dbg_nesc
     !end if

     do i=1,3
        xdp = nb6_xp(i,1:nb6_npart) + rbd_xc(i)
        write(ilun) xdp
     end do

     do i=1,3
        xdp = nb6_vp(i,1:nb6_npart)
        write(ilun) xdp
     end do
     xdp = nb6_mp(1:nb6_npart)
     write(ilun)xdp
     
     deallocate(xdp)
     close(ilun)
  end if
end subroutine rbd_output_part
