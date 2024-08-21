subroutine rbd_build_mesh
  use amr_commons
  use rbd_commons
  use pm_commons

  implicit none
  
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2, scale_M
  real(dp), parameter :: cm2pc  = 3.240779270005395D-19
  real(dp) :: dx, min_dx, half_scale, px, py
  real(dp), dimension(1:3) :: gx0
  real(dp) :: max_dist = 0.0
  integer :: pid, i, j, k, ierr, counter

  include 'mpif.h'

  write(6,*) 'RBD : Generating force mesh'
  rbd_mesh_np = rbd_mesh_Nx**3+1

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)
  if (myid==1) then
     ! Calculating spacing between mesh points
     dx = rbd_mesh_scale / (rbd_mesh_Nx - 1)
     half_scale = rbd_mesh_scale * 0.5

     ! And we loop
     rbd_mesh_pos(:, 1) = 0.0 !rbd_xc 
     pid = 2
     do i=1, rbd_mesh_Nx
        px = (i-1)*dx - half_scale
        !px = (i-1)*dx + gx0(1)
        do j=1, rbd_mesh_Nx
           py = (j-1)*dx - half_scale
           !py = (j-1)*dx + gx0(2)
           do k=1, rbd_mesh_Nx
              rbd_mesh_pos(1, pid) = px
              rbd_mesh_pos(2, pid) = py
              rbd_mesh_pos(3, pid) = (k-1)*dx - half_scale
              !rbd_mesh_pos(3, pid) = (k-1)*dx + gx0(3)
              pid = pid + 1
           end do
        end do
     end do

     ! Debug
     ! Counting how many are outside of the mesh scale
     !counter = 0
     !do i=1, nb6_npart
     !   do k=1,3
     !      if (nb6_xp(k, i) < rbd_mesh_pos(k, 2) .or. nb6_xp(k, i) > rbd_mesh_pos(k, pid-1)) then
     !         counter = counter + 1
     !      end if
     !   end do
     !end do
     !write(6,*) 'MESH BUILDING : ', counter, ' nbody particles are outside !'
     !call flush(6)
  end if

  ! Broadcasting the info to everyone
  call MPI_Bcast(rbd_mesh_scale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr) ! Necessary ? Not sure ...
  call MPI_Bcast(rbd_mesh_np, 1, MPI_INTEGER, 0, MPI_COMM_RAMSES, ierr)
  call MPI_Bcast(rbd_mesh_pos, 3*rbd_mesh_np, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)
end subroutine rbd_build_mesh

subroutine rbd_sync_mesh
  use amr_commons
  use rbd_commons
  use rbd_parameters

  implicit none

  include 'mpif.h'

  integer :: ierr, i, k, j
  integer, dimension(1:nvector) :: cc
  real(dp), dimension(1:3, 1:nvector) :: xx_dp
  real(dp) :: scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2
  real(dp), parameter :: pc2cm   = 3.086e18
  real(dp), parameter :: kms2cms = 1.0e5
  real(dp), parameter :: M2g     = 1.988e33
  real(dp), dimension(1:3) :: cur_xp

  if (rbd_gc_owner == -1) then
     call MPI_Abort(MPI_COMM_RAMSES, 110, ierr)
  end if

  call units(scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2)

  call timer('rbd_wait_mesh', 'start')
  if (myid==1) then
     if (rbd_dbg_comm) then
        write(6,*) 'RBD : Comm tag 8 (<- Mesh scale)'
        call flush(6)
     end if

     call MPI_Recv(rbd_mesh_scale, 1, MPI_DOUBLE_PRECISION, 0, 8, MPI_COMM_RAMBODY, MPI_STATUS_IGNORE, ierr)
     rbd_mesh_scale = rbd_mesh_scale * pc2cm / scale_l
  end if

  ! Broadcasting mesh scale to all ramses processes
  call MPI_Bcast(rbd_mesh_scale, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)

  ! We build the mesh on process 1
  ! TODO : Maybe test myid here ?
  call rbd_build_mesh

  ! Transmission of the coords
  call MPI_Bcast(rbd_mesh_pos, rbd_mesh_np*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_RAMSES, ierr)

  ! Init of mesh points
  rbd_xp = 0
  rbd_vp = 0
  rbd_mp = 0
  rbd_id = 0
  rbd_levelp = levelmin

  rbd_npart = 0
  
  do i=1, rbd_mesh_np
     cur_xp = rbd_mesh_pos(:,i) + rbd_xc(:)
     xx_dp(:,1) = cur_xp / boxlen

     ! TODO Vectorize this !!!
     call cmp_cpumap(xx_dp, cc, 1)
     if (cc(1) == myid) then
        rbd_npart = rbd_npart + 1
        rbd_xp(:, rbd_npart) = cur_xp
        rbd_id(rbd_npart) = i
     end if
  end do

  call timer('rbd_sync', 'start')
  ! And we push the points on the rbd grid
  call rbd_init_tree

end subroutine rbd_sync_mesh

! Debug/diag routine
subroutine rbd_build_force_profile
  use amr_commons
  use rbd_commons
  use rbd_parameters

  implicit none

  include 'mpif.h'
  
  integer::ierr, i, j, k, pm
  real(dp) :: min_x, max_x, dx
  real(dp), dimension(1:3) :: x, fc
  integer::Np=100

  min_x = -0.02 ! 1 kpc
  max_x =  0.02
  dx = (max_x - min_x) / (Np-1)

  open(unit=92548, file="rbd_force_data")
  do i=1, Np
     x(1) = min_x + i*dx
     do j=1, Np
        x(2) = min_x + j*dx
        do k=1, Np
           x(3) = min_x + k*dx

           call rbd_get_force_contribution(rbd_xc + x, fc)
           pm=0
           if (rbd_dbg_pm) pm=1
           
           write(92548, '(7e20.9, i5)') x, fc, norm2(fc), pm
        end do
     end do
  end do

  do i=1, rbd_mesh_np
     write(92548, '(7e20.9, i5)') rbd_mesh_pos(:,i), rbd_mesh_force(:,i), norm2(rbd_mesh_force(:,i)), 2
  end do

  close(92548)

  
  !call MPI_Abort(MPI_COMM_RAMSES, 3, ierr)
end subroutine rbd_build_force_profile

subroutine rbd_get_force_contribution(xi, fc)
  use rbd_commons
  use pm_commons
  use amr_commons
  implicit none

  include 'mpif.h'

  real(dp), dimension(1:3), intent(out) :: fc
  real(dp), dimension(1:3), intent(in)  :: xi

  real(dp), dimension(1:3) :: rij
  real(dp) :: dist
  logical :: inside_box
  integer :: idim
  ! Interpolation
  integer ii, ij, ik, ngx, ierr, i
  integer, dimension(1:8) :: iid
  real(dp) :: xd, yd, zd, half_box, grid_dx, mlim
  real(dp), dimension(1:3) :: c00, c01, c10, c11, c0, c1, p1, p2, xj
  integer, dimension(1:3) :: idx
  real(dp), dimension(1:8) :: NF
  
  ! Checking if point is inside the box grid
  inside_box = .true.
  xj = xi - rbd_xc ! xj is xi in the frame of the guiding center

  do idim=1,ndim
     if (xj(idim) < rbd_mesh_pos(idim, 2) .or. xj(idim) > rbd_mesh_pos(idim, rbd_mesh_np)) then
        inside_box = .false.
     end if
  end do

  ! If we're inside, we interpolate the force
  p1 = rbd_mesh_pos(:,2)
  p2 = rbd_mesh_pos(:,rbd_mesh_np)

  if (inside_box) then
     rbd_dbg_pm = .false.
     
     ngx = rbd_mesh_Nx

     idx = (xj-p1)*(ngx-1) / (p2-p1)
     ii = idx(1)
     ij = idx(2)
     ik = idx(3)

     ! Ward : this should NEVER happen !
     if (ii < 0 .or. ii >= ngx-1 .or. ij < 0 .or. ij >= ngx-1 .or. ik < 0 .or. ik >= ngx-1) then
        write(6,*) 'ERROR: trying to ask for interpolated force outside of the grid mesh !'
        write(6,*) 'X  = ', xi(:)
        write(6,*) 'XJ = ', xj(:)
        write(6,*) 'XC = ', rbd_xc(:)
        write(6,*) 'ii, ij, ik = ', ii, ij, ik
        write(6,*) 'ngx = ', ngx
        write(6,*) 'Min mesh point = ', rbd_mesh_pos(:,2) 
        write(6,*) 'Max mesh point = ', rbd_mesh_pos(:,rbd_mesh_np) 
        call MPI_Abort(MPI_COMM_RAMSES, 112, ierr)
     end if

     ! Copied from the nbody6 interpolation
     iid(1) = ii*ngx**2 + ij*ngx + ik + 2 
     iid(2) = iid(1)+1
     iid(3) = iid(1)+ngx
     iid(4) = iid(1)+ngx+1
     iid(5) = iid(1)+ngx**2
     iid(6) = iid(1)+ngx**2+1
     iid(7) = iid(1)+ngx**2+ngx
     iid(8) = iid(1)+ngx**2+ngx+1

     xd = (xj(1)-rbd_mesh_pos(1,iid(1))) / (rbd_mesh_pos(1,iid(8))-rbd_mesh_pos(1,iid(1)))

     yd = (xj(2)-rbd_mesh_pos(2,iid(1))) / (rbd_mesh_pos(2,iid(8))-rbd_mesh_pos(2,iid(1)))

     zd = (xj(3)-rbd_mesh_pos(3,iid(1))) / (rbd_mesh_pos(3,iid(8))-rbd_mesh_pos(3,iid(1)))

     ! Trilinear interpolation
     do idim=1, 3
        c00(idim) = rbd_mesh_force(idim, iid(1))*(1.0-xd) + rbd_mesh_force(idim, iid(5))*xd
        c01(idim) = rbd_mesh_force(idim, iid(2))*(1.0-xd) + rbd_mesh_force(idim, iid(6))*xd
        c10(idim) = rbd_mesh_force(idim, iid(3))*(1.0-xd) + rbd_mesh_force(idim, iid(7))*xd
        c11(idim) = rbd_mesh_force(idim, iid(4))*(1.0-xd) + rbd_mesh_force(idim, iid(8))*xd

        c0(idim) = c00(idim)*(1.0-yd) + c10(idim)*yd
        c1(idim) = c01(idim)*(1.0-yd) + c11(idim)*yd

        fc(idim) = c0(idim)*(1.0-zd)+c1(idim)*zd
     end do

     !if (rbd_debug_mesh) then
     !   do i=1, 8
     !      NF(i) = norm2(rbd_mesh_force(:,iid(i)))
     !   end do
     !   write(6,*) 'T  = ', t
     !   write(6,*) 'X  = ', xi(:)
     !   write(6,*) 'XJ = ', xj(:)
     !   write(6,*) 'Min mesh point = ', rbd_mesh_pos(:,2) 
     !   write(6,*) 'Max mesh point = ', rbd_mesh_pos(:,rbd_mesh_np)
     !   write(6,*) 'ii, ij, ik = ', ii, ij, ik
     !   write(6,*) 'iid = ', iid
     !   write(6,*) 'xd, yd, zd = ', xd, yd, zd
     !   write(6,*) 'NF = ', NF
     !   write(6,*) 'fc = ', fc
     !   write(6,*) '||fc|| = ', norm2(fc)
     !end if
     
  else ! Else we consider the point mass case
     dist = sqrt(xj(1)**2.0+xj(2)**2.0+xj(3)**2.0)
     fc = -xj(:) * nb6_tot_mass / dist**3.0
     
     !if (rbd_debug_mesh) then         
     !   write(6,*) 'Dist=', dist
     !   write(6,*) 'Mass=', nb6_tot_mass
     !   write(6,*) 'FC=', fc
     !end if
  end if
  
end subroutine rbd_get_force_contribution
