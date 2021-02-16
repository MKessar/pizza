!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the 1-d advection diffusion problem on a periodic domain

!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
module main_LIBPFASST_CALL 
!   use pf_mod_mpi
use parallel_mod
  use iso_c_binding  
  use pf_mod_dtype
  use pf_mod_stop
!   use pf_mod_zndarray
!   use pizza_zdnarray
  use pizza_zndsysarray
  integer ::  ierror


contains
  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run_pfasst_pizza(time_comm,timers)  
    use pfasst  !< This module has include statements for the main pfasst routines
    use pf_my_sweeper  !< Local module for sweeper and function evaluations
    use pf_my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
!     use probin  !< Local module reading/parsing problem parameters
   use namelists, only: imex_stat, n_modes_m,n_points_phi,n_points_r,nsteps,dt, l_chem
    use blocking, only: nMstart, nMstop!,nRstart, nRstop
       use timers_mod, only: timers_type

    implicit none

    type(timers_type),   intent(inout) :: timers
    integer           ::  time_comm   

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
!     type(pf_zndarray_t):: y_0      !<  the initial condition
!     type(pf_zndarray_t):: y_end    !<  the solution at the final time
    type(pizza_zndsysarray_t):: y_0      !<  the initial condition
    type(pizza_zndsysarray_t):: y_end    !<  the solution at the final time
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters

    integer           ::  l   !  loop variable over levels
    integer           ::  mpibuflen  !  Length of MPI buffer
    integer           ::  grid_shape(3)   !  size of the spatial discretization
    real(8) :: runStart, runStop


!   namelist /params/   imex_stat
    pf_fname="input.nml"
     

    !>  Set up communicator
    call pf_mpi_create(comm, time_comm)

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    !> Loop over levels and set some level specific parameters

    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
       allocate(pizza_zndsysarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)

       !>  Allocate the shape array for level (here just one dimension)
       if ( l_chem ) then

          grid_shape=[nMstop-nMstart+1,n_points_r(l),5]
       else
          grid_shape=[nMstop-nMstart+1,n_points_r(l),4]
       endif
       
       mpibuflen=(product(grid_shape))*2 ! The two is because data is complex
       
       !>  Set the size of the data on this level (here just one)
       call pf_level_set_size(pf,l,grid_shape,mpibuflen)

    end do

    
    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)

    
    !> Add some hooks for output
!     call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)
!     call pf_add_hook(pf, -1, PF_POST_CONVERGENCE, time_series)
    call pf_add_hook(pf, -1, PF_POST_BLOCK, time_series)
!     call pf_add_hook(pf, -1, PF_POST_SWEEP, time_series)
!     call pf_add_hook(pf, -1, PF_POST_ITERATION, time_series)
!     print*, "Before pf_print_options"

    !>  Output the run options 
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
!     call print_loc_options(pf,un_opt=6)

    !> Allocate initial and final solutions


    call pizza_zndsysarray_build(y_0, [ pf%levels(pf%nlevels)%lev_shape])
    call pizza_zndsysarray_build(y_end, [ pf%levels(pf%nlevels)%lev_shape])
       
    !> compute initial condition
    
    call initial_condition(y_0)


    !> Do the PFASST stepping
    runStart = MPI_Wtime()
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
    runStop = MPI_Wtime()
    
    timers%pfasst_run   = timers%pfasst_run   + (runStop-runStart)
    timers%n_pfasst_run = timers%n_pfasst_run + nsteps

    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

    !>  Deallocate initial condition and final solution
    call pizza_zndsysarray_destroy(y_0)
    

    call pizza_zndsysarray_destroy(y_end)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run_pfasst_pizza

end module
