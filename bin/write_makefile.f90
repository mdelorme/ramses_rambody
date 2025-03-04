subroutine output_makefile(filename)
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=30)::format
  integer::ilun

  ilun=11

  fileloc=TRIM(filename)
  format="(A)"
  open(unit=ilun,file=fileloc,form='formatted')
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# If you have problems with this makefile, contact Romain.Teyssier@gmail.com"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# Compilation time parameters"
  write(ilun,format)""
  write(ilun,format)"# Do we want a debug build? 1=Yes, 0=No"
  write(ilun,format)"DEBUG = 0"
  write(ilun,format)"# Compiler flavor: GNU or INTEL"
  write(ilun,format)"COMPILER = GNU"
  write(ilun,format)"# Size of vector cache"
  write(ilun,format)"NVECTOR = 32"
  write(ilun,format)"# Number of dimensions"
  write(ilun,format)"NDIM = 1"
  write(ilun,format)"# Float precision size"
  write(ilun,format)"NPRE = 8"
  write(ilun,format)"# hydro/mhd/rhd solver"
  write(ilun,format)"SOLVER = hydro"
  write(ilun,format)"# Patch"
  write(ilun,format)"PATCH ="
  write(ilun,format)"# Use RT? 1=Yes, 0=No"
  write(ilun,format)"RT = 0"
  write(ilun,format)"# Use turbulence? 1=Yes, 0=No (requires fftw3)"
  write(ilun,format)"USE_TURB = 0"
  write(ilun,format)"# Use MPI? 1=Yes, 0=No"
  write(ilun,format)"MPI = 0"
  write(ilun,format)"MPIF90 = mpif90"
  write(ilun,format)"# Root name of executable"
  write(ilun,format)"EXEC = ramses"
  write(ilun,format)"# Number of additional energies"
  write(ilun,format)"NENER = 0"
  write(ilun,format)"# Use Grackle cooling? 1=Yes, 0=No"
  write(ilun,format)"GRACKLE = 0"
  write(ilun,format)"# Number of metal species"
  write(ilun,format)"NMETALS = 0"
  write(ilun,format)"# Use ATON? 1=Yes, 0=No"
  write(ilun,format)"ATON = 0"
  write(ilun,format)"# Number of ions for RT"
  write(ilun,format)"# use 1 for HI+HII, +1 for added H2, +2 for added HeI, HeII, HeIII"
  write(ilun,format)"NIONS = 0"
  write(ilun,format)"# Number of photon groups for RT"
  write(ilun,format)"NGROUPS = 0"
  write(ilun,format)"# Number of passive scalars"
  write(ilun,format)"NPSCAL = 0"
  write(ilun,format)"# Light MPI communicator structure (for > 10k MPI processes)"
  write(ilun,format)"LIGHT_MPI_COMM = 0"
  write(ilun,format)""
  write(ilun,format)"# Compute NVAR"
  write(ilun,format)"NVAR = 2+$(NDIM)"
  write(ilun,format)"ifeq ($(SOLVER),mhd)"
  write(ilun,format)"   NVAR = 8"
  write(ilun,format)"endif"
  write(ilun,format)"NVAR := $(NVAR)+$(NENER)+$(NPSCAL)+$(NMETALS)"
  write(ilun,format)"ifeq ($(RT),1)"
  write(ilun,format)"   NVAR := $(NVAR)+$(NIONS)"
  write(ilun,format)"endif"
  write(ilun,format)""
  write(ilun,format)"# Set to one to use 'include ""mpif.h""' instead of more recent ""use mpi"""
  write(ilun,format)"OLD_MPI_SUPPORT = 0"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)"
  write(ilun,format)"GITHASH = $(shell git log --pretty=format:'%H' -n 1)"
  write(ilun,format)"GITREMOTE = $(shell git config --get branch.$(GITBRANCH).remote)"
  write(ilun,format)"GITREPO = $(shell git config --get remote.$(GITREMOTE).url)"
  write(ilun,format)"BUILDDATE = $(shell date +""%D-%T"")"
  write(ilun,format)"DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) \"
  write(ilun,format)"          -DNVAR=$(NVAR) -DSOLVER$(SOLVER)"
  write(ilun,format)"ifeq ($(ATON),1)"
  write(ilun,format)"   DEFINES += -DATON"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(GRACKLE),1)"
  write(ilun,format)"   DEFINES += -Dgrackle"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(OLD_MPI_SUPPORT),1)"
  write(ilun,format)"   DEFINES += -DMPI_OLD"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(RT),1)"
  write(ilun,format)"   DEFINES += -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(USE_TURB),1)"
  write(ilun,format)"   DEFINES += -DUSE_TURB"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(LIGHT_MPI_COMM), 1)"
  write(ilun,format)"   DEFINES += -DLIGHT_MPI_COMM"
  write(ilun,format)"endif"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# Fortran compiler options and directives"
  write(ilun,format)""
  write(ilun,format)"# GNU compiler (gfortran)"
  write(ilun,format)"ifeq ($(COMPILER),GNU)"
  write(ilun,format)"   FFLAGS = -cpp $(DEFINES)"
  write(ilun,format)"   ifeq ($(MPI),1)"
  write(ilun,format)"      F90 = $(MPIF90)"
  write(ilun,format)"      FC = gfortran"
  write(ilun,format)"   else"
  write(ilun,format)"      F90 = gfortran"
  write(ilun,format)"      FC = gfortran"
  write(ilun,format)"      FFLAGS += -DWITHOUTMPI"
  write(ilun,format)"   endif"
  write(ilun,format)"   F90 += -ffree-line-length-none -fimplicit-none"
  write(ilun,format)"   FC += -ffree-line-length-none -fimplicit-none"
  write(ilun,format)"   ifeq ($(DEBUG),1)"
  write(ilun,format)"      F90 += -g -O0 -fbacktrace -fbounds-check -Wuninitialized -Wall"
  write(ilun,format)"      FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan"
  write(ilun,format)"   else"
  write(ilun,format)"      F90 += -O3"
  write(ilun,format)"   endif"
  write(ilun,format)"endif"
  write(ilun,format)""
  write(ilun,format)"# Intel compiler"
  write(ilun,format)"ifeq ($(COMPILER),INTEL)"
  write(ilun,format)"   FFLAGS = -cpp $(DEFINES)"
  write(ilun,format)"   ifeq ($(MPI),1)"
  write(ilun,format)"      F90 = $(MPIF90)"
  write(ilun,format)"      FC = ifort"
  write(ilun,format)"      FFLAGS += -DNOSYSTEM"
  write(ilun,format)"   else"
  write(ilun,format)"      F90 = ifort"
  write(ilun,format)"      FC = ifort"
  write(ilun,format)"      FFLAGS += -DWITHOUTMPI"
  write(ilun,format)"   endif"
  write(ilun,format)"   F90 += -fp-model source"
  write(ilun,format)"   FC += -fp-model source"
  write(ilun,format)"   ifeq ($(DEBUG),1)"
  write(ilun,format)"      F90 += -warn all -O0 -g -traceback -check bounds"
  write(ilun,format)"      FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays"
  write(ilun,format)"   else"
  write(ilun,format)"      F90 += -O3"
  write(ilun,format)"   endif"
  write(ilun,format)"endif"
  write(ilun,format)""
  write(ilun,format)"#############################################################################"
  write(ilun,format)"MOD = mod"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# MPI librairies"
  write(ilun,format)"LIBMPI ="
  write(ilun,format)"#LIBMPI = -lfmpi -lmpi -lelan"
  write(ilun,format)""
  write(ilun,format)"# --- CUDA libraries, for Titane ---"
  write(ilun,format)"LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart"
  write(ilun,format)""
  write(ilun,format)"ifeq ($(GRACKLE),1)"
  write(ilun,format)"   # Add include and library install path for grackle and hdf5 here"
  write(ilun,format)"   LIBS_GRACKLE = -L$(HOME)/local/lib -lgrackle -lhdf5 -lz -lgfortran -ldl"
  write(ilun,format)"   LIBS_OBJ     = -I$(HOME)/local/include -DCONFIG_BFLOAT_8 -DH5_USE_16_API -fPIC"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(USE_TURB),1)"
  write(ilun,format)"   LIBS_TURB = -L$(HOME)/local/lib -lfftw3"
  write(ilun,format)"   LIBS_OBJ_TURB = -I$(HOME)/local/include -lfftw3"
  write(ilun,format)"endif"
  write(ilun,format)"LIBS = $(LIBMPI) $(LIBS_GRACKLE) $(LIBS_TURB)"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# Sources directories are searched in this exact order"
  write(ilun,format)"VPATH = $(PATCH):../$(SOLVER):../aton:"
  write(ilun,format)"ifeq ($(RT),1)"
  write(ilun,format)"   VPATH += ../rt:"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(USE_TURB),1)"
  write(ilun,format)"   VPATH += ../turb:"
  write(ilun,format)"endif"
  write(ilun,format)"VPATH += ../hydro:../pm:../poisson:../amr:../io:../rbd"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"# All objects"
  write(ilun,format)"MODOBJ = mpi_mod.o amr_parameters.o amr_commons.o random.o pm_parameters.o sink_feedback_parameters.o \"
  write(ilun,format)"         pm_commons.o poisson_parameters.o dump_utils.o constants.o file_module.o rbd_parameters.o rbd_commons.o"
  write(ilun,format)"ifeq ($(GRACKLE),1)"
  write(ilun,format)"   MODOBJ += grackle_parameters.o"
  write(ilun,format)"endif"
  write(ilun,format)"MODOBJ += poisson_commons.o hydro_parameters.o hydro_commons.o \"
  write(ilun,format)"          cooling_module.o bisection.o sparse_mat.o clfind_commons.o \"
  write(ilun,format)"          gadgetreadfile.o write_makefile.o write_patch.o write_gitinfo.o \"
  write(ilun,format)"          sink_sn_feedback.o"
  write(ilun,format)"ifeq ($(RT),1)"
  write(ilun,format)"   MODOBJ += rt_parameters.o rt_hydro_commons.o coolrates_module.o \"
  write(ilun,format)"             rt_spectra.o rt_cooling_module.o rt_flux_module.o"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(USE_TURB),1)"
  write(ilun,format)"   MODOBJ += turb_parameters.o turb_commons.o"
  write(ilun,format)"endif"
  write(ilun,format)""
  write(ilun,format)"AMROBJ = read_params.o init_amr.o init_time.o init_refine.o tracer_utils.o adaptive_loop.o \"
  write(ilun,format)"         amr_step.o update_time.o output_amr.o flag_utils.o \"
  write(ilun,format)"         physical_boundaries.o virtual_boundaries.o refine_utils.o \"
  write(ilun,format)"         nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o \"
  write(ilun,format)"         eos.o units.o light_cone.o movie.o memory.o end.o"
  write(ilun,format)"# Particle-Mesh objects"
  write(ilun,format)"PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o \"
  write(ilun,format)"        newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o \"
  write(ilun,format)"        sink_particle.o feedback.o clump_finder.o clump_merger.o output_clump.o \"
  write(ilun,format)"        flag_formation_sites.o init_sink.o output_sink.o \"
  write(ilun,format)"        unbinding.o merger_tree.o move_tracer.o init_tracer.o \"
  write(ilun,format)"        read_sink_feedback_params.o sink_rt_feedback.o \"
  write(ilun,format)"        stellar_particle.o init_stellar.o output_stellar.o"
  write(ilun,format)""
  write(ilun,format)"# Rambody"
  write(ilun,format)"RBDOBJ = rbd_clean.o rbd_debug.o rbd_escapers.o rbd_force_fine.o \"
  write(ilun,format)"	 rbd_init.o rbd_list.o rbd_mesh.o rbd_output_mesh.o rbd_output_part.o \"
  write(ilun,format)"	 rbd_ptree.o rbd_sync_cluster.o rbd_sync_forces.o \"
  write(ilun,format)"	 rbd_synchro_fine.o rbd_sync_timestep.o"
  write(ilun,format)"# Poisson solver objects"
  write(ilun,format)"POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o \"
  write(ilun,format)"             multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o \"
  write(ilun,format)"             multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o \"
  write(ilun,format)"             output_poisson.o"
  write(ilun,format)"# Hydro objects"
  write(ilun,format)"HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o \"
  write(ilun,format)"           courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o \"
  write(ilun,format)"           godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o \"
  write(ilun,format)"           read_hydro_params.o synchro_hydro_fine.o cooling_module_ism.o"
  write(ilun,format)"# RT objects"
  write(ilun,format)"RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o \"
  write(ilun,format)"        rt_output_hydro.o rt_godunov_fine.o rt_interpol_hydro.o \"
  write(ilun,format)"        rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o \"
  write(ilun,format)"        rt_boundana.o rt_units.o"
  write(ilun,format)"# Turbulence objects"
  write(ilun,format)"TURBOBJ = turb_force_utils.o read_turb_params.o init_turb.o mpi_share_turb_fields.o \"
  write(ilun,format)"          turb_next_field.o turb_check_time.o write_turb_fields.o read_turb_fields.o \"
  write(ilun,format)"          add_turb_forcing.o"
  write(ilun,format)"# Patch objects"
  write(ilun,format)"sinclude $(PATCH)/Makefile"
  write(ilun,format)""
  write(ilun,format)"# All objects"
  write(ilun,format)"AMRLIB = $(AMROBJ) $(RBDOBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) "
  write(ilun,format)""
  write(ilun,format)"ifeq ($(RT),1)"
  write(ilun,format)"   AMRLIB += $(RTOBJ)"
  write(ilun,format)"endif"
  write(ilun,format)"ifeq ($(USE_TURB),1)"
  write(ilun,format)"   AMRLIB += $(TURBOBJ)"
  write(ilun,format)"endif"
  write(ilun,format)"# ATON objects"
  write(ilun,format)"ATON_MODOBJ = timing.o radiation_commons.o rad_step.o"
  write(ilun,format)"ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o \"
  write(ilun,format)"           rad_backup.o ../aton/atonlib/libaton.a"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"ramses:	$(MODOBJ) $(AMRLIB) ramses.o"
  write(ilun,format)"	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)"
  write(ilun,format)"	rm write_makefile.f90"
  write(ilun,format)"	rm write_patch.f90"
  write(ilun,format)"ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o"
  write(ilun,format)"	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o \"
  write(ilun,format)"		-o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)"
  write(ilun,format)"	rm write_makefile.f90"
  write(ilun,format)"	rm write_patch.f90"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"write_gitinfo.o: FORCE"
  write(ilun,format)"	$(FC) -O0 -cpp -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' \"
  write(ilun,format)"		-DGITHASH=\'""$(GITHASH)""\' -DGITREPO=\'$(GITREPO)\' \"
  write(ilun,format)"		-DBUILDDATE=\'""$(BUILDDATE)""\' -c ../amr/write_gitinfo.f90 -o $@"
  write(ilun,format)"write_makefile.o: FORCE"
  write(ilun,format)"	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)"
  write(ilun,format)"	$(FC) -O0 -c write_makefile.f90 -o $@"
  write(ilun,format)"write_patch.o: FORCE"
  write(ilun,format)"	../utils/scripts/cr_write_patch.sh $(PATCH)"
  write(ilun,format)"	$(FC) -O0 -c write_patch.f90 -o $@"
  write(ilun,format)"%.o:%.F"
  write(ilun,format)"	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ) $(LIBS_OBJ_TURB)"
  write(ilun,format)"%.o:%.f90"
  write(ilun,format)"	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ) $(LIBS_OBJ_TURB)"
  write(ilun,format)"FORCE:"
  write(ilun,format)"#############################################################################"
  write(ilun,format)"clean:"
  write(ilun,format)"	rm -f *.o *.$(MOD) *.i"
  write(ilun,format)"#############################################################################"

  close(ilun)
end subroutine output_makefile
