! RAMBODY Version : No WITHOUTMPI flag !
#ifdef WITHOUTMPI
#error "This Ramses version is a special RAMBODY version, you have to compile it with MPI activated !"
#endif

program ramses
  use rbd_parameters
  implicit none

  ! Read run parameters
  call read_params

#ifndef CRAY
  ! Set signal handler
  call set_signal_handler
#endif

  ! if we are in rambody mode we set it up
  if (rambody) then
     write(6,*) 'RBD : Initialising Rambody !'
     call init_rambody
  end if

  ! Start time integration
  call adaptive_loop

end program ramses

#ifndef CRAY
! sets the hook to catch signal 10, doesn't work with CRAY
subroutine set_signal_handler
  implicit none
  external output_signal
  integer::istatus=-1
#ifdef NOSYSTEM
  integer::jsigact,jhandle
#endif

#ifndef NOSYSTEM
  call SIGNAL(10,output_signal,istatus)
#else
  call PXFSTRUCTCREATE("sigaction",jsigact,istatus)
  call PXFGETSUBHANDLE(output_signal,jhandle,istatus)
  call PXFINTSET(jsigact,"sa_handler",jhandle,istatus)
  call PXFSIGACTION(10,jsigact,0,istatus)
#endif
end subroutine set_signal_handler

! signal handler subroutine
subroutine output_signal
  use amr_commons
  implicit none

  if (myid==1) write (*,*) 'SIGNAL 10: Output will be written to disk during next main step.'

  ! output will be written to disk at next main step
  output_now = .true.

end subroutine output_signal

! Initialisation of Rambody
subroutine init_rambody
  use amr_commons
  use pm_commons
  use rbd_commons
  implicit none
  include 'mpif.h'

  integer::ierr, mpi_st
  integer:: greetings, k
  real(dp) :: dt_ram, dt_nb6

#ifdef WITHOUTMPI
  ! ERROR !
  write(*,*) ' ERROR ! Using parameter rambody=.true. but the code has not been'
  write(*,*) ' compiled with MPI ! Recompile the code, removing the -DWITHOUTMPI flag'
  call clean_stop
#endif
  if (myid .eq. 1) then
     ! Sending greetings message to NBody6
     greetings = 94111535

     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 1 (-> Greetings)'
        call flush(6)
     end if
     
     call MPI_Send(greetings, 1, MPI_INTEGER, 0, 1, MPI_COMM_RAMBODY, ierr)

     ! Waiting for acknowledgment
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 2 (<- Greetings)'
        call flush(6)
     end if
     
     call MPI_Recv(greetings, 1, MPI_INTEGER, 0, 2, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)

     ! Checking greetings 94111535
     if (greetings .ne. 1180276) then
        write(6,*) 'RMS : Message received is not proper greeting. Something might be wrong ...'
        write(6,*) '   - Received from NBody6 : ', greetings
        call MPI_Abort(MPI_COMM_WORLD, 125)
     end if
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

end subroutine init_rambody

