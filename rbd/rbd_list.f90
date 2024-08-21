!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_add_list(ind_part,list2,ok,np)
  use amr_commons
  use rbd_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part,list2
  logical,dimension(1:nvector)::ok
  !
  ! Add particles to their new linked lists
  !
  integer::j

  do j=1,np
     if(ok(j))then
        if(rbd_numbp(list2(j))>0)then
           ! Add particle at the tail of its linked list
           rbd_nextp(rbd_tailp(list2(j)))=ind_part(j)
           rbd_prevp(ind_part(j))=rbd_tailp(list2(j))
           rbd_nextp(ind_part(j))=0
           rbd_tailp(list2(j))=ind_part(j)
           rbd_numbp(list2(j))=rbd_numbp(list2(j))+1
        else
           ! Initialise linked list
           rbd_headp(list2(j))=ind_part(j)
           rbd_tailp(list2(j))=ind_part(j)
           rbd_prevp(ind_part(j))=0
           rbd_nextp(ind_part(j))=0
           rbd_numbp(list2(j))=1
        end if
        rbd_npart = rbd_npart + 1
        rbd_total_added = rbd_total_added + 1
     end if
  end do

end subroutine rbd_add_list

!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_remove_list(ind_part,list1,ok,np)
  use amr_commons
  use rbd_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part,list1
  logical,dimension(1:nvector)::ok
  !----------------------------------------------------
  ! Remove particles from their original linked lists
  !----------------------------------------------------
  integer::j
  do j=1,np
     if(ok(j))then
        if(rbd_prevp(ind_part(j)) .ne. 0) then
           if(rbd_nextp(ind_part(j)) .ne. 0 )then
              rbd_nextp(rbd_prevp(ind_part(j)))=rbd_nextp(ind_part(j))
              rbd_prevp(rbd_nextp(ind_part(j)))=rbd_prevp(ind_part(j))
           else
              rbd_nextp(rbd_prevp(ind_part(j)))=0
              rbd_tailp(list1(j))=rbd_prevp(ind_part(j))
           end if
        else
           if(rbd_nextp(ind_part(j)) .ne. 0)then
              rbd_prevp(rbd_nextp(ind_part(j)))=0
              rbd_headp(list1(j))=rbd_nextp(ind_part(j))
           else
              rbd_headp(list1(j))=0
              rbd_tailp(list1(j))=0
           end if
        end if
        rbd_numbp(list1(j))=rbd_numbp(list1(j))-1
        rbd_npart = rbd_npart - 1
        rbd_total_removed = rbd_total_removed + 1
     end if
  end do
end subroutine rbd_remove_list

!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_add_free(ind_part,np)
  use amr_commons
  use rbd_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        rbd_xp(idim, ind_part(j))=0.0
        rbd_vp(idim, ind_part(j))=0.0
     end do
  end do
  do j=1,np
     rbd_mp(ind_part(j))=0.0
     rbd_id(ind_part(j))=0
     rbd_levelp(ind_part(j))=0
  end do

  do j=1,np
     if(rbd_numbp_free>0)then
        ! Add particle at the tail of its linked list
        rbd_nextp(rbd_tailp_free)=ind_part(j)
        rbd_prevp(ind_part(j))=rbd_tailp_free
        rbd_nextp(ind_part(j))=0
        rbd_tailp_free=ind_part(j)
        rbd_numbp_free=rbd_numbp_free+1
     else
        ! Initialise linked list
        rbd_headp_free=ind_part(j)
        rbd_tailp_free=ind_part(j)
        rbd_prevp(ind_part(j))=0
        rbd_nextp(ind_part(j))=0
        rbd_numbp_free=1
     end if
  end do
  rbd_npart=rbdpartmax-rbd_numbp_free

end subroutine rbd_add_free
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_add_free_cond(ind_part,ok,np)
  use amr_commons
  use rbd_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  logical,dimension(1:nvector)::ok
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           rbd_xp(idim, ind_part(j))=0.0
           rbd_vp(idim, ind_part(j))=0.0
        endif
     end do
  end do
  do j=1,np
     if(ok(j))then
        rbd_mp(ind_part(j))=0.0
        rbd_id(ind_part(j))=0
        rbd_levelp(ind_part(j))=0
     endif
  end do
  
  do j=1,np
     if(ok(j))then
        if(rbd_numbp_free>0)then
           ! Add particle at the tail of its linked list
           rbd_nextp(rbd_tailp_free)=ind_part(j)
           rbd_prevp(ind_part(j))=rbd_tailp_free
           rbd_nextp(ind_part(j))=0
           rbd_tailp_free=ind_part(j)
           rbd_numbp_free=rbd_numbp_free+1
        else
           ! Initialise linked list
           rbd_headp_free=ind_part(j)
           rbd_tailp_free=ind_part(j)
           rbd_prevp(ind_part(j))=0
           rbd_nextp(ind_part(j))=0
           rbd_numbp_free=1
        end if
     endif
  end do
  rbd_npart=rbdpartmax-rbd_numbp_free

end subroutine rbd_add_free_cond
!################################################################
!################################################################
!################################################################
!################################################################
subroutine rbd_remove_free(ind_part,np)
  use amr_commons
  use rbd_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  !-----------------------------------------------
  ! Get np particle from free memory linked list
  !-----------------------------------------------
  integer::j,ipart
  do j=1,np
     ipart=rbd_headp_free
     ind_part(j)=ipart
     rbd_numbp_free=rbd_numbp_free-1
     if(rbd_numbp_free<0)then
        write(*,*)'No more free memory'
        write(*,*)'in PE ',myid
        write(*,*)'Increase npartmax'
        call clean_stop
     end if
     rbd_headp_free=rbd_nextp(rbd_headp_free)
  end do
  !rbd_npart=rbdpartmax-rbd_numbp_free
end subroutine rbd_remove_free
