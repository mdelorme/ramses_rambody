module rbd_parameters
  use amr_parameters
  
  integer :: rbdpartmax      = 1000000           ! Maximum number of particles in the NBody simulation
  logical :: rambody         = .false.           ! Is Rambody on ?
  logical :: rbd_restart     = .false.
  character(len=15)::sync_method = 'coarse_dt' ! Sychronisation method for Rambody
  real(dp), dimension(1:3) :: rbd_xc0          ! Initial position for the cluster
  real(dp), dimension(1:3) :: rbd_vc0          ! Initial velocity for the cluster
  logical :: refine_on_rambody = .true.       ! Should there be some refinement on rambody particles
  integer :: rbd_refine_npart = 1              

  integer :: rbd_mesh_Nx=10
  integer :: rbd_escapers_buf=5000
  
  logical :: rbd_dbg_comm = .true. ! Messaging the communication tags
  logical :: rbd_dbg_barriers = .false.
  logical :: rbd_nbody6_force = .false.

  integer :: max_nb6_steps = 100 ! To avoid having too many nbody6 steps in a single Ramses step
  logical :: rbd_limit_dt = .false.
  real(dp) :: rbd_max_rms_dt = 1.0

  logical :: rbd_store_escapers = .true.
  real(dp) :: rbd_epsilon = 1.0e-2 ! Smoothing length

end module rbd_parameters
