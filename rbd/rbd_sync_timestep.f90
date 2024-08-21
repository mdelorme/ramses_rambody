subroutine rbd_get_nb6dt
  use amr_commons
  use rbd_commons
  use rbd_parameters

  implicit none
  include 'mpif.h'

  integer :: ierr

  call timer('rbd_wait_timestep', 'start')
  
  ! The first process recovers NBody6's timestep
  if (myid == 1) then
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 3 (<- NB6 dt)'
        call flush(6)
     end if

     call MPI_Recv(rbd_recv_buf, 1, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)
  end if

  ! Then broadcasts it to everyone
  call MPI_Bcast(rbd_recv_buf, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
  rbd_nb6_dt = rbd_recv_buf(1, 1)

  if (myid == 1) then
     write(6,*) 'RBD : NBody6 timestep is ', rbd_nb6_dt, ' Myrs'
     call flush(6)
  end if

  call timer('rbd_sync', 'start')
end subroutine rbd_get_nb6dt

subroutine rbd_sync_timestep(ilevel)
  use amr_commons
  use rbd_commons
  use rbd_parameters
  
  implicit none
  include 'mpif.h'

  integer, intent(in) :: ilevel
  integer  :: ierr, istatus, i, k
  real(dp) :: dt_r
  real(dp) :: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(dp) :: rbd_new_nb6_dt, limited_dt
  real(dp), parameter :: myrs_to_s = 3.15576D13
  ! Setting the timestep of the level so that it's compliant with NBody6's

  call timer('rbd_wait_timestep', 'start')
  if (myid==1) then
     ! Conversion factor from user units to cgs units                              
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     
     ! Ratio between the dt_n and dt_r
     dt_r = dtnew(ilevel) * scale_t / myrs_to_s
     write(6,*) '---- RBD timestep synchronization ----'
     write(6,*) 'Ramses DT = ', dtnew(ilevel), '; In Physical units : ', dt_r, 'Myrs'
     rbd_dt_k = dt_r / rbd_nb6_dt

     k      = ceiling(rbd_dt_k)
     write(6,*) 'NBody6 Repetitions, k = ', k
     if (k > max_nb6_steps) then
        k = max_nb6_steps
        write(6,*) 'Limiting to the maximum number of steps. k = ', k
     end if
     rbd_next_sync = t + (k * rbd_nb6_dt * myrs_to_s / scale_t)

     ! Limiting dt if over the limit
     if (rbd_next_sync - t < dtnew(ilevel)) then
        limited_dt = rbd_next_sync - t
        write(6,*) 'Limiting ramses timestep. Old value=', dtnew(ilevel), '; New value=', limited_dt
        dtnew(ilevel) = limited_dt
     end if

     if (rbd_dbg_comm) then
        write(6,*) 'Comm tag 11 (-> k steps)'
        call flush(6)
     end if

     call MPI_Send(k, 1, MPI_INTEGER, 0, 11, MPI_COMM_RAMBODY, ierr)
     
     write(6,*) 'Sync time : ', t*scale_t/myrs_to_s, 'Myrs'
     write(6,*) 'Next sync time : ', rbd_next_sync
     write(6,*) '--------------------------------------'
  end if

  ! Sending info to the rest of Ramses
  call MPI_Bcast(rbd_next_sync, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
  call timer('rbd_sync', 'start')
  
end subroutine rbd_sync_timestep
