subroutine rbd_clean
  use rbd_commons
  implicit none
  
  deallocate(rbd_recv_buf)
  deallocate(rbd_send_buf)
  deallocate(rbd_xp)
  deallocate(rbd_vp)
  deallocate(rbd_fp)
  deallocate(rbd_mp)
  deallocate(rbd_id)
  deallocate(rbd_escapers)

  deallocate(nb6_xp)
  deallocate(nb6_vp)
  deallocate(nb6_mp)

  deallocate(rbd_nextp)
  deallocate(rbd_prevp)
  !deallocate(rbd_levelp)

  deallocate(rbd_headp)
  deallocate(rbd_tailp)
  deallocate(rbd_numbp)

  deallocate(rbd_mesh_pos)
  deallocate(rbd_mesh_force)
end subroutine rbd_clean
  
