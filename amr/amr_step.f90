recursive subroutine amr_step(ilevel,icount)
  use amr_commons
  use rbd_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use tracer_utils, only: reset_tracer_move_flag
#ifdef RT
  use rt_hydro_commons
  use SED_module
  use UV_module
  use coolrates_module, only: update_coolrates_tables
  use rt_cooling_module, only: update_UVrates
#endif
  use sink_feedback_parameters, only: sn_feedback_sink
#if USE_TURB==1
  use turb_commons
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::mpi_err
#endif
  integer, intent(in)::ilevel,icount
  !-------------------------------------------------------------------!
  ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
  ! Each routine is called using a specific order, don't change it,   !
  ! unless you check all consequences first.                          !
  !-------------------------------------------------------------------!
  integer::i,idim,ivar
  logical::ok_defrag,output_now_all
  logical,save::first_step=.true.
  logical::is_sync_time

  if(numbtot(1,ilevel)==0)return

  if(verbose)write(*,999)icount,ilevel

  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  call timer('refine','start')
  if(levelmin.lt.nlevelmax .and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
     if(ilevel==levelmin.or.icount>1)then
        do i=ilevel,nlevelmax
           if(i>levelmin)then

              !--------------------------
              ! Build communicators
              !--------------------------
              call build_comm(i)

              !--------------------------
              ! Update boundaries
              !--------------------------
              call make_virtual_fine_int(cpu_map(1),i)
              if(hydro)then
#ifdef SOLVERmhd
                 do ivar=1,nvar+3
#else
                 do ivar=1,nvar
#endif
                    call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
                 end do
#else
                 end do
#endif
                 if(momentum_feedback>0)call make_virtual_fine_dp(pstarold(1),i)
                 if(strict_equilibrium>0)call make_virtual_fine_dp(rho_eq(1),i)
                 if(strict_equilibrium>0)call make_virtual_fine_dp(p_eq(1),i)
                 if(simple_boundary)call make_boundary_hydro(i)
              end if
#ifdef RT
              if(rt)then
                 do ivar=1,nrtvar
                    call make_virtual_fine_dp(rtuold(1,ivar),i)
                 end do
                 if(simple_boundary)call rt_make_boundary_hydro(i)
              end if
#endif
              if(poisson)then
                 call make_virtual_fine_dp(phi(1),i)
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),i)
                 end do
                 if(simple_boundary)call make_boundary_force(i)
              end if
           end if

           !--------------------------
           ! Refine grids
           !--------------------------
           call refine_fine(i)
        end do
     end if
  end if

  !--------------------------
  ! Load balance
  !--------------------------
  call timer('loadbalance','start')
  ok_defrag=.false.

  if(levelmin.lt.nlevelmax)then
     if(ilevel==levelmin)then
        if(nremap>0)then
           ! Skip first load balance because it has been performed before file dump
           if(nrestart>0.and.first_step)then
              if(nrestart.eq.nrestart_quad) restart_remap=.true.
              if(restart_remap) then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
              first_step=.false.
           else
              if(MOD(nstep_coarse,nremap)==0)then
                 call load_balance
                 call defrag
                 ok_defrag=.true.
              endif
           end if
        end if
     endif
  end if

  !-----------------
  ! Rambody, 1st sync
  !-----------------
  if (rambody) then
     call timer('rbd_sync', 'start')
     if (t > rbd_next_sync .and. t > 0 .and. rbd_sync_state > 0) then
        write(6,'(a, f15.7, a, f15.7)') 'WARNING, Current time (', t, ') > next sync (', rbd_next_sync, ') but no sync in view ... : '
     end if
     
     if (ilevel == levelmin .and. rbd_sync_state == 0) then
        rbd_sync_state    = 1
        rbd_skip_subcycle = .false.
        rbd_sync_needed   = .true.

        ! Syncing cluster, and putting the particles on ilevel        
        call rbd_sync_cluster ! Getting cluster positions and velocities
        call rbd_sync_mesh    ! Getting the scale of the cluster
        call rbd_get_nb6dt    ! Getting original timestep
     end if
  end if

  !-----------------
  ! Update sink cloud particle properties
  !-----------------
#if NDIM==3
  call timer('sinks','start')
  if(sink)call update_cloud(ilevel)
#endif

  !-----------------
  ! Particle leakage
  !-----------------
  call timer('particles','start')
  if(pic)call make_tree_fine(ilevel)

  !------------------------
  ! Output results to files
  !------------------------
  if(ilevel==levelmin)then

#ifdef WITHOUTMPI
     output_now_all = output_now
#else
     ! check if any of the processes received a signal for output
     call MPI_BARRIER(MPI_COMM_RAMSES,mpi_err)
     call MPI_ALLREDUCE(output_now,output_now_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_RAMSES,mpi_err)
#endif
     if(mod(nstep_coarse,foutput)==0.or.aexp>=aout(iout).or.t>=tout(iout) &
        &.or.aexp>=aout_next.or.t>=tout_next.or.output_now_all.EQV..true.)then

        call timer('io','start')
        if(.not.ok_defrag)then
           call defrag
        endif

        ! Run the clumpfinder, (produce output, don't keep arrays alive on output)
        ! CAREFUL: create_output is used to destinguish between the case where
        ! the clumpfinder is called from create_sink or directly from amr_step.
#if NDIM==3
        if(clumpfind .and. ndim==3) call clump_finder(.true.,.false.)
#endif
        
        call dump_all
        if (output_now_all.EQV..true.) then
          output_now=.false.
          if (finish_run) then
            ! trick to stop the code after walltime triggered output
            tout(iout)=t
            aout(iout)=aexp
          endif
        endif

     endif
  
     ! Dump lightcone
     if(lightcone .and. ndim==3) call output_cone()

  endif

  !----------------------------
  ! Output frame to movie dump (without synced levels)
  !----------------------------
  if(movie) then
     if(imov.le.imovout)then
        if(aexp>=amovout(imov).or.t>=tmovout(imov))then
           call timer('io','start')

           write(6,*) 'Starting movie output'
           call flush(6)

           call output_frame()

           write(6,*) 'End movie output'
           call flush(6)

        endif
     endif
  end if

  !-----------------------------------------------------------
  ! Put here all stuffs that are done only at coarse time step
  !-----------------------------------------------------------
  if(ilevel==levelmin)then
     !----------------------------------------------------
     ! Kinetic feedback from giant molecular clouds
     !----------------------------------------------------
     call timer('feedback','start')
     if(hydro.and.star.and.eta_sn>0.and.f_w>0)call kinetic_feedback
  endif

  !----------------------------------------------------
  ! Feedback on sink particles
  !----------------------------------------------------
  if(stellar) then
     call make_stellar_from_sinks
  endif
  if (sn_feedback_sink) then
     call make_sn_stellar
  endif

  !--------------------
  ! Poisson source term
  !--------------------
  if(poisson)then
     call timer('poisson','start')
     !save old potential for time-extrapolation at level boundaries
     call save_phi_old(ilevel)
     call timer('rho','start')
     call rho_fine(ilevel,icount)
  endif

  !-------------------------------------------
  ! Sort particles between ilevel and ilevel+1
  !-------------------------------------------
  if(pic)then
     ! Remove particles to finer levels
     call timer('particles','start')
     call kill_tree_fine(ilevel)
     ! Update boundary conditions for remaining particles
     call virtual_tree_fine(ilevel)
  end if

  !---------------
  ! Gravity update
  !---------------
  if(poisson)then
     call timer('poisson','start')

     ! Remove gravity source term with half time step and old force
     if(hydro)then
        call synchro_hydro_fine(ilevel,-0.5*dtnew(ilevel),1)
     endif

     ! Compute gravitational potential
     if(ilevel>levelmin)then
        if(ilevel .ge. cg_levelmin) then
           call phi_fine_cg(ilevel,icount)
        else
           call multigrid_fine(ilevel,icount)
        end if
     else
        call multigrid_fine(levelmin,icount)
     end if
     !when there is no old potential...
     if (nstep==0)call save_phi_old(ilevel)

     ! Compute gravitational acceleration
     call force_fine(ilevel,icount)

     ! Synchronize remaining particles for gravity
     if(pic)then
        call timer('particles','start')
        if(static_dm.or.static_stars)then
           call synchro_fine_static(ilevel)
        else
           call synchro_fine(ilevel)

           if (rambody) then
             call timer('rbd_synchro_fine', 'start')
             call rbd_synchro_fine(ilevel)
           end if
        end if
     end if

     if(hydro)then
        call timer('poisson','start')

        ! Add gravity source term with half time step and new force
        call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel),1)

        ! Update boundaries
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
        end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)

        ! Compute Bondi-Hoyle accretion parameters
#if NDIM==3
        call timer('sinks','start')
        if(sink.and.hydro)call collect_acczone_avg(ilevel)
#endif
     end if
  end if

#ifdef RT
  ! Turn on RT in case of rt_stars and first stars just created:
  ! Update photon packages according to star particles and sink particles
  call timer('radiative transfer','start')
  if(rt .and. rt_star) call update_star_RT_feedback(ilevel)
#if NDIM==3
  if(rt .and. rt_sink) call update_sink_RT_feedback
#endif
#endif

#if USE_TURB==1
  ! Compute turbulent forcing
                               call timer('turb','start')
  if (turb .and. turb_type/=3) then
     ! Calculate turbulent acceleration on each cell in this level
     call calc_turb_forcing(ilevel)
  end if
#endif

  !----------------------
  ! Compute new time step
  !----------------------
  call timer('courant','start')
  call newdt_fine(ilevel)
  if(ilevel>levelmin)then
     dtnew(ilevel)=MIN(dtnew(ilevel-1)/real(nsubcycle(ilevel-1)),dtnew(ilevel))
  end if

  ! Set unew equal to uold
  call timer('hydro - set unew','start')
  if(hydro)call set_unew(ilevel)

#ifdef RT
  ! Set rtunew equal to rtuold
  call timer('radiative transfer','start')
  if(rt)call rt_set_unew(ilevel)
#endif

  ! Syncing rambody now that dt has been computed
  if (rambody) then
     call timer('rbd_sync_timestep', 'start')
     if (ilevel == rbd_level_min .and. rbd_sync_state == 1) then
        ! Synchronising timesteps
        call rbd_sync_timestep(ilevel)
        if (ilevel == levelmin) write(6,*) 'RBD : Timestep synchronization; T = ', t, 'Next sync = ', rbd_next_sync
        rbd_sync_state = 2
     end if
  end if

  !---------------------------
  ! Recursive call to amr_step
  !---------------------------
  if(ilevel<nlevelmax)then
     ! Are there sub-levels ?
     if(numbtot(1,ilevel+1)>0)then
        if(nsubcycle(ilevel)==2)then
           call amr_step(ilevel+1,1)

           if (rambody .and. rbd_skip_subcycle) then
              dtnew(ilevel) = rbd_next_sync - t
           else
              call amr_step(ilevel+1,2)
           end if
        else
           call amr_step(ilevel+1,1)
        endif
     else
        ! Rambody, checking timestep does not go further than the next sync time
        if (rambody .and. (t + dtnew(ilevel) - dtnew(ilevel)**rbd_margin >= rbd_next_sync) .and. rbd_sync_state == 3) then
           write(6,'(a, E15.7, a, E15.7, a, E15.7)') ' RBD : Going over next sync time. Current time : ', t, '; Next sync : ', rbd_next_sync, '; dt = ', dtnew(ilevel)
           write(6,*) 'RBD : Changing dt from ', dtnew(ilevel), ' to ', rbd_next_sync - t
           dtnew(ilevel) = rbd_next_sync - t

           rbd_sync_state    = 0
           rbd_skip_subcycle = .true.
        end if
        
        ! Otherwise, update time and finer level time-step
        dtold(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        dtnew(ilevel+1)=dtnew(ilevel)/dble(nsubcycle(ilevel))
        call update_time(ilevel)
#if NDIM==3
        if(sink)call update_sink(ilevel)
#endif
     end if
  else
     ! Rambody, checking timestep does not go further than the next sync time
     if (rambody .and. (t + dtnew(ilevel) - dtnew(ilevel)**rbd_margin >= rbd_next_sync) .and. rbd_sync_state == 3) then
        write(6,'(a, E15.7, a, E15.7, a, E15.7)') ' RBD : Going over next sync time. Current time : ', t, '; Next sync : ', rbd_next_sync, '; dt = ', dtnew(ilevel)
        write(6,*) 'RBD : Changing dt from ', dtnew(ilevel), ' to ', rbd_next_sync - t
        dtnew(ilevel) = rbd_next_sync - t
        rbd_sync_state = 0
        rbd_skip_subcycle = .true.
     end if
     call update_time(ilevel)
     if(sink)call update_sink(ilevel)
  end if

  ! Thermal feedback from stars
  call timer('feedback','start')
  if(hydro.and.star.and.eta_sn>0)call thermal_feedback(ilevel)

  ! Density threshold or Bondi accretion onto sink particle
#if NDIM==3
  if(sink.and.hydro)then
     call timer('sinks','start')
     call grow_sink(ilevel,.false.)
  end if
#endif
  !-----------
  ! Hydro step
  !-----------
  if ((hydro).and.(.not.static_gas)) then

     ! Hyperbolic solver
     call timer('hydro - godunov','start')
     call godunov_fine(ilevel)

     ! Reverse update boundaries
     call timer('hydro - rev ghostzones','start')
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_reverse_dp(unew(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     ! MC Tracer
     ! Communicate fluxes accross boundaries
     if(MC_tracer)then
        call timer('tracer','start')
        do ivar=1,twondim
           call make_virtual_reverse_dp(fluxes(1,ivar),ilevel-1)
           call make_virtual_fine_dp(fluxes(1,ivar),ilevel-1)
        end do
     end if

     if(momentum_feedback>0)then
        call make_virtual_reverse_dp(pstarnew(1),ilevel)
     endif
     if(pressure_fix)then
        call make_virtual_reverse_dp(enew(1),ilevel)
        call make_virtual_reverse_dp(divu(1),ilevel)
     endif

     ! Set uold equal to unew
     call timer('hydro - set uold','start')
     call set_uold(ilevel)

     ! Add gravity source term with half time step and old force
     ! in order to complete the time step 
     call timer('poisson','start')
     if(poisson)call synchro_hydro_fine(ilevel,+0.5*dtnew(ilevel),1)

#if USE_TURB==1
     ! Compute turbulent forcing
     call timer('turb','start')
     if (turb .AND. turb_type/=3) then
        ! Euler step, adding turbulent acceleration
        call synchro_hydro_fine(ilevel,dtnew(ilevel),2)
     end if
#endif

     ! Restriction operator
     call timer('hydro upload fine','start')
     call upload_fine(ilevel)

  endif

  !---------------------
  ! Do RT/Chemistry step
  !---------------------
#ifdef RT
  if(rt .and. rt_advect) then
     call timer('radiative transfer','start')
     call rt_step(ilevel)
  else
     ! Still need a chemistry call if RT is defined but not
     ! actually doing radiative transfer (i.e. rt==false):
     call timer('cooling','start')
     if(hydro .and. (neq_chem.or.cooling.or.T2_star>0.0.or.barotropic_eos))call cooling_fine(ilevel)
  endif
  ! Regular updates and book-keeping:
  if(ilevel==levelmin) then
     call timer('radiative transfer','start')
     if(cosmo) call update_rt_c
     if(cosmo .and. haardt_madau) call update_UVrates(aexp)
     if(cosmo .and. rt_isDiffuseUVsrc) call update_UVsrc
     call timer('cooling','start')
     if(cosmo) call update_coolrates_tables(dble(aexp))
     call timer('radiative transfer','start')
     if(ilevel==levelmin) call output_rt_stats
  endif
#else
  call timer('cooling','start')
  if((hydro).and.(.not.static_gas)) then
    if(neq_chem.or.cooling.or.T2_star>0.0.or.barotropic_eos)call cooling_fine(ilevel)
  endif
#endif

  !---------------
  ! Move particles
  !---------------
  if(pic)then
     call timer('particles','start')
     if(static_dm.or.static_stars)then
        call move_fine_static(ilevel) ! Only remaining particles
     else
        call move_fine(ilevel) ! Only remaining particles
     end if

     call timer('rbd_move_part', 'start')
     if (rambody .and. ilevel == rbd_gc_level) then
        call rbd_sync_gc(.false.)
     end if

  end if

  !----------------------------------
  ! Star formation in leaf cells only
  !----------------------------------
#if NDIM==3
  call timer('feedback','start')
  if(hydro.and.star.and.(.not.static_gas))call star_formation(ilevel)
#endif
  !---------------------------------------
  ! Update physical and virtual boundaries
  !---------------------------------------
  if((hydro).and.(.not.static_gas))then
     call timer('hydro - ghostzones','start')
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
     end do
#else
     end do
#endif
     if(momentum_feedback>0)call make_virtual_fine_dp(pstarold(1),ilevel)
     if(strict_equilibrium>0)call make_virtual_fine_dp(rho_eq(1),ilevel)
     if(strict_equilibrium>0)call make_virtual_fine_dp(p_eq(1),ilevel)
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

#ifdef SOLVERmhd
  ! Magnetic diffusion step
  if((hydro).and.(.not.static_gas))then
     if(eta_mag>0d0.and.ilevel==levelmin)then
        call timer('hydro - diffusion','start')
        call diffusion
     endif
  end if
#endif


  ! Sending the force back to the cluster, and moving the guiding center
  if (rambody) then
     call timer('rbd_sync_forces', 'start')
     ! Level max to make sure we send this at the beginning of integration and not after the timestep !
     is_sync_time = .false. ! Why order priority does not work and am I forced to do this ?
     if (ilevel == nlevelmax) then
      is_sync_time = .true.
     else if (numbtot(1, ilevel+1) == 0) then
      is_sync_time = .true.
     endif
     if (is_sync_time .and. rbd_sync_state == 2) then
        call rbd_sync_forces
        rbd_sync_state = 3
     end if
  end if
  
  !-----------------------
  ! Compute refinement map
  !-----------------------
  call timer('flag','start')
  if(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)) call flag_fine(ilevel,icount)

  !----------------------------
  ! Merge finer level particles
  !----------------------------
  call timer('particles','start')
  if(pic)call merge_tree_fine(ilevel)

  !---------------
  ! Radiation step
  !---------------
#ifdef ATON
  if(aton.and.ilevel==levelmin)then
     call timer('aton','start')
     call rad_step(dtnew(ilevel))
  endif
#endif

  if(sink)then
     call timer('sinks','start')
     !-------------------------------
     ! Update coarser level sink velocity
     !-------------------------------
     if(ilevel>levelmin)then
        vsold(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel-1)
        if(nsubcycle(ilevel-1)==1)vsnew(1:nsink,1:ndim,ilevel-1)=vsnew(1:nsink,1:ndim,ilevel)
        if(icount==2)vsnew(1:nsink,1:ndim,ilevel-1)= &
             (vsold(1:nsink,1:ndim,ilevel)*dtold(ilevel)+vsnew(1:nsink,1:ndim,ilevel)*dtnew(ilevel))/ &
             (dtold(ilevel)+dtnew(ilevel))
     end if
     !---------------
     ! Sink production
     !---------------
#if NDIM==3
     if(ilevel==levelmin)call create_sink
#endif
  end if

  !-------------------------------
  ! Update coarser level time-step
  !-------------------------------
  if(ilevel>levelmin)then
     if(nsubcycle(ilevel-1)==1)dtnew(ilevel-1)=dtnew(ilevel)
     if(icount==2)dtnew(ilevel-1)=dtold(ilevel)+dtnew(ilevel)
  end if

  ! Reset move flag flag
  if(MC_tracer) then
                                call timer('tracer','start')
     ! Decrease the move flag by 1
     call reset_tracer_move_flag(ilevel)
  end if

999 format(' Entering amr_step(',i1,') for level',i2)

end subroutine amr_step

!##########################################################################
!##########################################################################
!##########################################################################
!##########################################################################

#ifdef RT
subroutine rt_step(ilevel)
  use amr_parameters, only: dp
  use amr_commons,    only: levelmin, t, dtnew, myid
  use rt_cooling_module, only: update_UVrates
  use rt_hydro_commons
  use UV_module
  use SED_module,     only: star_RT_feedback
  use mpi_mod
  implicit none
  integer, intent(in) :: ilevel

!--------------------------------------------------------------------------
!  Radiative transfer and chemistry step. Either do one step on ilevel,
!  with radiation field updates in coarser level neighbours, or, if
!  rt_nsubsteps>1, do many substeps in ilevel only, using Dirichlet
!  boundary conditions for the level boundaries.
!--------------------------------------------------------------------------

  real(dp) :: dt_hydro, t_left, dt_rt, t_save
  integer  :: i_substep, ivar

  dt_hydro = dtnew(ilevel)                   ! Store hydro timestep length
  t_left = dt_hydro
  ! We shift the time backwards one hydro-dt, to get evolution of stellar
  ! ages within the hydro timestep, in the case of rt subcycling:
  t_save=t ; t=t-t_left

  i_substep = 0
  do while (t_left > 0)                      !                RT sub-cycle
     i_substep = i_substep + 1
     call get_rt_courant_coarse(dt_rt)
     ! Temporarily change timestep length to rt step:
     dtnew(ilevel) = MIN(t_left, dt_rt/2**(ilevel-levelmin))
     t = t + dtnew(ilevel) ! Shift the time forwards one dt_rt

     ! If (myid==1) write(*,900) dt_hydro, dtnew(ilevel), i_substep, ilevel
     if (i_substep > 1) call rt_set_unew(ilevel)

     if(rt_star) call star_RT_feedback(ilevel,dtnew(ilevel))
#if NDIM==3
     if(rt_sink) call sink_RT_feedback(ilevel,dtnew(ilevel))
#endif

     ! Hyperbolic solver
     if(rt_advect) call rt_godunov_fine(ilevel,dtnew(ilevel))

     call add_rt_sources(ilevel,dtnew(ilevel))

     ! Reverse update boundaries
     do ivar=1,nrtvar
        call make_virtual_reverse_dp(rtunew(1,ivar),ilevel)
     end do

     ! Set rtuold equal to rtunew
     call rt_set_uold(ilevel)

     call timer('cooling', 'start')
     if(neq_chem.or.cooling.or.T2_star>0.0.or.barotropic_eos)call cooling_fine(ilevel)
     call timer('radiativbe transfer', 'start')
     do ivar=1,nrtvar
        call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
     end do
     if(simple_boundary)call rt_make_boundary_hydro(ilevel)

     t_left = t_left - dtnew(ilevel)
  end do                                   !          End RT subcycle loop
  dtnew(ilevel) = dt_hydro                 ! Restore hydro timestep length
  t = t_save       ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update coarser level split cells
  call rt_upload_fine(ilevel)

  if (myid==1 .and. rt_nsubcycle .gt. 1) write(*,901) ilevel, i_substep

  !900 format (' dt_hydro=', 1pe12.3, ' dt_rt=', 1pe12.3, ' i_sub=', I5, ' level=', I5)
901 format (' Performed level', I3, ' RT-step with ', I5, ' subcycles')

end subroutine rt_step
#endif
