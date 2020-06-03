!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

!   real(pfdp), save :: v      ! advection velocity
!   real(pfdp), save :: nu     ! viscosity
!   real(pfdp), save :: Pr     ! Prandl number
!   integer,    save :: kfreq  ! initial condition parameter
  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time

  integer, save :: n_modes_m(PF_MAXLEVS)     ! number of grid points
  integer, save :: n_points_phi(PF_MAXLEVS)     ! number of grid points
  integer, save :: n_points_r(PF_MAXLEVS)     ! number of grid points
  integer, save :: nsteps          ! number of time steps
  integer, save :: imex_stat       ! type of imex splitting
  integer, save :: nprob           ! which problem
! n_modes_m,n_points_phi,n_points_r
  character(len=128), save :: pfasst_nml

  character(len=64), save :: output ! directory name for output
  CHARACTER(LEN=255) :: istring  ! stores command line argument
  CHARACTER(LEN=255) :: message           ! use for I/O error messages

  integer :: ios,iostat
  namelist /params/  n_modes_m,n_points_r,n_points_phi,nprob, nsteps, dt, Tfin
  namelist /params/  pfasst_nml, a,b, nu, t00, sigma, kfreq,imex_stat

contains

  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    character(len=128) :: arg  !  command line argument
    character(128)    :: probin_fname   !<  file name for input parameters
    
    !> Set the name of the input file
!     probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    probin_fname = "input.nml" ! default file name - can be overwritten on the command line
!     if (command_argument_count() >= 1) &
!          call get_command_argument(1, value=probin_fname)
    
    !> set defaults
!     nsteps  = -1
!     v       = 1.0_pfdp
    a       = 1.0_pfdp
    b       = 1.0_pfdp
    nu      = 0.01_pfdp
    kfreq   = 1.0_pfdp
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    nsteps = 1
    imex_stat=2    !  Default is full IMEX
    pfasst_nml=probin_fname
!     print*,"probin_fname=",probin_fname

    !>  Read in stuff from input file
    un = 9
    write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
    open(unit=un, file = probin_fname, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    !>  Read the command line
    i = 0
    do
       call get_command_argument(i, arg)
       if (LEN_TRIM(arg) == 0) return !EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          read(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/(1.0_8*nsteps)

    !  Return the name of the file from which to read PFASST parameters
    pf_fname=pfasst_nml
    print*,"pf_fname=",pf_fname
  end subroutine probin_init

  subroutine print_loc_options(pf, un_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    integer :: un = 6

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt

    !  Print out the local parameters
    write(un,*) '=================================================='
    write(un,*) ' '
    write(un,*) 'Local Variables'
    write(un,*) '----------------'
    write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    write(un,*) 'Dt:     ', Dt, '! Time step size'
    write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
    write(un,*) 'nr:     ',  n_points_r(1:pf%nlevels), '! grid size per level'
    write(un,*) 'nm:     ',  n_modes_m(1:pf%nlevels), '! grid size per level'
    write(un,*) 'nphi:   ',  n_points_phi(1:pf%nlevels), '! grid size per level'
!     write(un,*) 'v:      ',  v, '! advection constant'
!     write(un,*) 'nu:     ', nu, '! diffusion constant'
    select case (imex_stat)
    case (0)  
       write(un,*) 'imex_stat:', imex_stat, '! Fully explicit'
    case (1)  
       write(un,*) 'imex_stat:', imex_stat, '! Fully implicit'
    case (2)  
       write(un,*) 'imex_stat:', imex_stat, '! Implicit/Explicit'
    case DEFAULT
       print *,'Bad case for imex_stat in probin ', imex_stat
!        call exit(0)
    end select

    
    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
