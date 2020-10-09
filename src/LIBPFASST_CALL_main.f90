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
  subroutine run_pfasst_pizza(time_comm)  
    use pfasst  !< This module has include statements for the main pfasst routines
    use pf_my_sweeper  !< Local module for sweeper and function evaluations
    use pf_my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
!     use probin  !< Local module reading/parsing problem parameters
   use namelists, only: imex_stat, n_modes_m,n_points_phi,n_points_r,nsteps,dt
    use blocking, only: nMstart, nMstop!,nRstart, nRstop
    implicit none

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator (here it is mpi)
!     type(pf_zndarray_t):: y_0      !<  the initial condition
!     type(pf_zndarray_t):: y_end    !<  the solution at the final time
    type(pizza_zndsysarray_t):: y_0      !<  the initial condition
    type(pizza_zndsysarray_t):: y_end    !<  the solution at the final time
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters

    integer           ::  l   !  loop variable over levels
    integer           ::  time_comm   !  loop variable over levels
    integer           ::  mpibuflen  !  Length of MPI buffer
    integer           ::  grid_shape(3)   !  size of the spatial discretization
!     complex(pfdp),         pointer :: omega(:,:), psi(:,:), theta(:,:), uphi0(:,:)
    complex(pfdp),         pointer :: omega(:,:), u_phi(:,:), theta(:,:), u_s(:,:)

!     integer           ::  un   !  loop variable over levels
!   integer, save :: n_modes_m(PF_MAXLEVS)     ! number of grid points
!   integer, save :: n_points_phi(PF_MAXLEVS)     ! number of grid points
!   integer, save :: n_points_r(PF_MAXLEVS)     ! number of grid points
!   integer, save :: nsteps          ! number of time steps
!   integer, save :: imex_stat       ! type of imex splitting
!   integer, save :: nprob           ! which problem
! ! n_modes_m,n_points_phi,n_points_r
! !   character(len=128), save :: pfasst_nml
! 
!   character(len=64), save :: output ! directory name for output
!   CHARACTER(LEN=255) :: istring  ! stores command line argument
!   CHARACTER(LEN=255) :: message           ! use for I/O error messages
!   real(pfdp), save :: dt     ! time step
!   real(pfdp), save :: Tfin   ! Final time

!   namelist /params/  n_modes_m,n_points_r,n_points_phi,nprob, nsteps, dt, Tfin
!   namelist /params/   imex_stat
    pf_fname="input.nml"
     
!     print*, "Before pf_mpi_create"

    !>  Set up communicator
    call pf_mpi_create(comm, time_comm)
!     print*, "Before pf_pfasst_create"

    !>  Create the pfasst structure
    call pf_pfasst_create(pf, comm, fname=pf_fname)

    !> Loop over levels and set some level specific parameters

    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(ad_level_t::pf%levels(l)%ulevel)
!        allocate(pf_zndarray_factory_t::pf%levels(l)%ulevel%factory)
!        allocate(pizza_zndarray_factory_t::pf%levels(l)%ulevel%factory)
       allocate(pizza_zndsysarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(ad_sweeper_t::pf%levels(l)%ulevel%sweeper)
!     print*, "Before pf_level_set_size"
       !>  Allocate the shape array for level (here just one dimension)
       grid_shape=[nMstop-nMstart+1,n_points_r(l),4]
       mpibuflen=(product(grid_shape))*2 ! The two is because data is complex
       
       !>  Set the size of the data on this level (here just one)
!        call pf_level_set_size(pf,l,[nMstop-nMstart+1,n_points_r(l)])
       call pf_level_set_size(pf,l,grid_shape,mpibuflen)

    end do
!     print*, "Before pf_pfasst_setup"
    
    !>  Set up some pfasst stuff
    call pf_pfasst_setup(pf)
!     print*, "Before pf_add_hook"

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
!     print*, "Before zndarray_build"

!     call pizza_zndarray_build(y_0, [ pf%levels(pf%nlevels)%lev_shape])
!     call pizza_zndarray_build(y_end, [ pf%levels(pf%nlevels)%lev_shape])
    call pizza_zndsysarray_build(y_0, [ pf%levels(pf%nlevels)%lev_shape])
    call pizza_zndsysarray_build(y_end, [ pf%levels(pf%nlevels)%lev_shape])

!     select type (y_0)
!     class is (pizza_zndsysarray_t)
!               print*,"size y_0%arr_shape",y_0%arr_shape
!               print*,"size y_0%flatarray",SIZE(y_0%flatarray)
!               
!               print*,"size y_end%arr_shape",y_end%arr_shape
!               print*,"size y_end%flatarray",SIZE(y_end%flatarray)
! 
!        omega  => get_array2d(y_0,1)
!        u_phi  => get_array2d(y_0,2)
!        u_s    => get_array2d(y_0,3)
!        theta  => get_array2d(y_0,4)
!    
!               
!     print*,"size omega",SIZE(omega)!,"this%omega_Mloc   =",SIZE(this%omega_Mloc)
!     print*,"size u_phi",SIZE(u_phi)!,"this%u_phi_m_Mloc =",SIZE(this%u_phi_m_Mloc)
!     print*,"size u_s"  ,SIZE(u_s)  !,"this%u_s_m_Mloc   =",SIZE(this%u_s_m_Mloc)
!     print*,"size theta",SIZE(theta)!,"this%theta_Mloc   =",SIZE(this%theta_Mloc)              
!     end select        
       
!     print*, "[ pf%levels(pf%nlevels)%lev_shape])=", [ pf%levels(pf%nlevels)%lev_shape]
    !> compute initial condition
!     print*, "pf%levels(pf%nlevels)%lev_shape",pf%levels(pf%nlevels)%lev_shape
!     print*, "Before initial_condition"
    
    call initial_condition(y_0)

!     print*, "Before pf_pfasst_run"

    !> Do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps,y_end)
!     print*, "Before mpi_barrier"

    !>  Wait for everyone to be done
    call mpi_barrier(pf%comm%comm, ierror)

!     print*, "Before zndarray_destroy a"

    !>  Deallocate initial condition and final solution
!     call pizza_zndarray_destroy(y_0)
    call pizza_zndsysarray_destroy(y_0)
    
!     print*, "Before zndarray_destroy b"

!     call pizza_zndarray_destroy(y_end)
    call pizza_zndsysarray_destroy(y_end)
!     print*, "Before pf_pfasst_destroy"

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)
!     print*, "after pf_pfasst_destroy"

  end subroutine run_pfasst_pizza

! 
!   function pizza_zndarray_norm(this,flags) result (norm)
!     class(pizza_zndarray_t), intent(in) :: this
!     integer,     intent(in   ), optional :: flags
!     real(pfdp) :: norm
!     real(pfdp) :: norm_local
! 
!     norm_local = maxval(abs(this%flatarray))
!     
!       call MPI_ALLREDUCE(norm, norm_local, 1, MPI_DEF_REAL, &
!            &          MPI_MAX, Comm_Pizza, ierr)
!   end function pizza_zndarray_norm

end module
