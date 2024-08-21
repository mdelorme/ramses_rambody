subroutine rbd_output_mesh(filename)

  use amr_commons
  use rbd_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer :: ilun, i
  
  if(verbose) write(*,*)' Entering rbd_output_mesh'

  ! Maybe unnecessary since rbd_output_part is called right before and also does the sync
  !call rbd_sync_gc(.false.)

  ! Only writing on master process
  if (myid == 1) then
     ilun = 2*ncpu+myid+10
     fileloc=TRIM(filename)
     open(unit=ilun, file=TRIM(fileloc), form='unformatted')
     rewind(ilun)

     write(ilun) rbd_xc
     write(ilun) rbd_vc
     write(ilun) rbd_gc_id
     write(ilun) rbd_mesh_np, rbd_mesh_Nx
     
     do i=1, 3
        write(ilun) rbd_mesh_pos(i, 1:rbd_mesh_np)
     end do

     do i=1, 3
        write(ilun) rbd_mesh_force(i, 1:rbd_mesh_np)
     end do
     close(ilun)
     
  end if
end subroutine rbd_output_mesh
