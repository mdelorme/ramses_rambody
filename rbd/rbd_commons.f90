module rbd_commons
  use amr_parameters
  use rbd_parameters
  
  real(dp) :: rbd_nb6_dt                                ! Nbody6 timestep
  logical  :: rbd_sync_done                             ! Has the synchronisation already been made on a finer level ?
  integer  :: rbd_nbody6_k                              ! Number of Nbody6 integrations
  real(dp), allocatable, dimension(:,:) :: rbd_recv_buf ! Communication buffer (only on process 0)
  real(dp), allocatable, dimension(:,:) :: rbd_send_buf ! Communication buffer (only on process 0)
  integer, parameter :: n_dbg_buf = 100
  real(dp), dimension(3,n_dbg_buf) :: rbd_dbg_buf  ! Debug buffer (only on process 0)

  ! Cluster particles : TODO : cleanup these only for mesh
  integer(i8b),  allocatable, dimension(:) :: rbd_id  
  real(dp), allocatable, dimension(:,:)    :: rbd_xp
  real(dp), allocatable, dimension(:,:)    :: rbd_vp
  real(dp), allocatable, dimension(:)      :: rbd_mp
  real(dp), allocatable, dimension(:,:)    :: rbd_fp

  ! NB6 particles
  real(dp), allocatable, dimension(:,:) :: nb6_xp
  real(dp), allocatable, dimension(:,:) :: nb6_vp
  real(dp), allocatable, dimension(:)   :: nb6_mp
  real(dp) :: nb6_tot_mass
  

  integer ,allocatable,dimension(:)  ::rbd_headp    ! Head rambody particle in grid
  integer ,allocatable,dimension(:)  ::rbd_tailp    ! Tail rambody particle in grid
  integer ,allocatable,dimension(:)  ::rbd_numbp    ! Number of rambody particles in grid
  integer ,allocatable,dimension(:)  ::rbd_nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::rbd_prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::rbd_levelp   ! Current level of particle
  integer :: rbd_level_min                          ! Minimum level where particles are put
  integer::rbd_headp_free, rbd_tailp_free, rbd_numbp_free=0, rbd_numbp_free_tot=0

  integer, dimension(1:20) :: rbd_p_per_level   ! DEBUG
  real(dp)                 :: rbd_dt_k          ! How many nbody6 timestep do we put in one ramses dt
  real(dp)                 :: rbd_next_sync     ! When is the next synchronisation time
  logical                  :: rbd_skip_subcycle ! Should we stop subcycling now
  integer, parameter       :: rbd_margin = 2    ! What margin on DT do we put. Because of floating point precision

  ! Cluster globals 
  real(dp), dimension(1:3) :: rbd_xc     ! Guiding center of the cluster
  real(dp), dimension(1:3) :: rbd_last_xc
  real(dp), dimension(1:3) :: rbd_vc     ! Velocity of the guiding center of the cluster
  real(dp), dimension(1:3) :: rbd_fc     ! Force applied on the guiding center

  real(dp)                 :: rbd_dt
  integer                  :: nb6_npart  ! How many particles on this process
  integer                  :: rbd_npart  ! Token to check list integrity. TODO : Remove this
  integer                  :: rbd_sync_state = 0
  logical                  :: rbd_sync_needed = .false.
  integer                  :: rbd_gc_owner, rbd_gc_level
  logical                  :: rbd_gc_force_synced

  integer                  :: rbd_gc_id

  ! Mesh
  integer  :: rbd_mesh_np
  real(dp) :: rbd_mesh_last_scale, rbd_mesh_scale
  real(dp), allocatable, dimension(:,:) :: rbd_mesh_pos
  real(dp), allocatable, dimension(:,:) :: rbd_mesh_force
  
  real(dp), dimension(1:3, 1:100) :: rbd_dbg_esc
  integer :: rbd_dbg_nesc

  ! Debug
  integer, dimension(1:3) :: bar_sync_count
  integer                 :: dbg_counter = 0
  logical                 :: rbd_debug_mesh = .false.
  logical                 :: rbd_dbg_pmass = .true.

  ! Book keeping
  integer                 :: rbd_total_added = 0
  integer                 :: rbd_total_removed = 0

  ! Escapers
  integer                                  :: rbd_tot_escapers = 0 ! In general
  integer                                  :: rbd_n_escapers   = 0 ! At this coarse step
  real(dp), allocatable, dimension(:,:)    :: rbd_escapers

  logical :: rbd_dbg_pm


end module rbd_commons
