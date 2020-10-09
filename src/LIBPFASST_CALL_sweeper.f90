!
! This file is part of LIBPFASST.
!
!
!> Sweeper and RHS routines for 1-D advection/diffusion example.
!>     u_t + v*u_x = nu*u_xx
module pf_my_sweeper
  use pf_mod_dtype
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
!   use pf_mod_zndarray
  use radial_scheme, only: type_rscheme
  use chebyshev, only: type_cheb
  use algebra, only: prepare_full_mat, solve_full_mat
  use useful, only: abortRun
  use truncation, only: minc
  use communications, only: transp_m2r, m2r_fields, transp_r2m, r2m_fields

  use parallel_mod
!   use pizza_zdnarray
  use pizza_zndsysarray
  use update_temp_coll, only: update_temp_co, get_temp_rhs_imp_coll,finish_exp_temp_coll
   use update_psi_coll_smat, only: update_om_coll_smat, finish_exp_psi_coll_smat
  use update_psi_coll_smat, only: get_psi_rhs_imp_coll_smat
  
  use rloop, only:radial_loop
  implicit none
  logical :: lMat=.true.
  logical :: l_log_next=.false.
  integer :: i_substep=1

!   integer ::  ierror

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: ad_sweeper_t
     integer ::     nr   !  Grid size
     integer ::     nm   !  Grid size

     !>  FFT and Spectral derivatives
!      type(pf_fft_t), pointer :: fft_tool
     complex(pfdp), allocatable :: u_s_m_Mloc(:,:)      !u_s velocity in complex space
     complex(pfdp), allocatable :: u_phi_m_Mloc(:,:)    !u_phi velocity in complex space
     complex(pfdp), allocatable :: u_s_m_Rloc(:,:)      !u_s velocity in complex space
     complex(pfdp), allocatable :: u_phi_m_Rloc(:,:)    !u_phi velocity in complex space
     complex(pfdp), allocatable :: theta_m_Rloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: omega_Mloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: theta_Mloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: psi_Mloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: xi_m_Rloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: omega_m_Rloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: uphi0_Mloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: fvec_omega_Mloc(:,:)       !u_phi velocity in complex space
     complex(pfdp), allocatable :: fvec_theta_Mloc(:,:)       !u_phi velocity in complex space
!      complex(pfdp), allocatable :: fvec_psi_Mloc(:,:)       !u_phi velocity in complex space
     complex(pfdp), allocatable :: fvec_u_s_Mloc(:,:)       !u_phi velocity in complex space
     complex(pfdp), allocatable :: fvec_u_phi_Mloc(:,:)       !u_phi velocity in complex space
!      complex(pfdp), allocatable :: fvec_uphi0_Mloc(:,:)       !u_phi velocity in complex space
     complex(pfdp), allocatable :: rhsvec_theta_Mloc(:,:)     !u_phi velocity in complex space
!      complex(pfdp), allocatable :: rhsvec_psi_Mloc(:,:)     !u_phi velocity in complex space
     complex(pfdp), allocatable :: rhsvec_u_s_Mloc(:,:)     !u_phi velocity in complex space
     complex(pfdp), allocatable :: rhsvec_u_phi_Mloc(:,:)     !u_phi velocity in complex space
     complex(pfdp), allocatable :: rhsvec_omega_Mloc(:,:)     !u_phi velocity in complex space
!      complex(pfdp), allocatable :: rhsvec_uphi0_Mloc(:,:)     !u_phi velocity in complex space
     real(pfdp), allocatable :: u_s_Rloc(:,:)           ! u_s velocity in real space
     real(pfdp), allocatable :: u_phi_Rloc(:,:)         !u_phi velocity in real space
     real(pfdp), allocatable :: om_Rloc(:,:)            !u_phi velocity in real space
     real(pfdp), allocatable :: theta_r_Rloc(:,:)       !temperature in real space
     real(pfdp), allocatable :: xi_Rloc(:,:)            !u_phi temperature in real space
     real(pfdp), allocatable :: tmp_Mloc(:,:)           ! used for operator in real space
     real(pfdp), allocatable :: tmp_bis_Mloc(:,:)       ! used for operator in real space
     real(pfdp), allocatable :: tmp_Rloc(:,:)           ! used for operator in real space
     real(pfdp), allocatable :: tmp_bis_Rloc(:,:)       ! used for operator in real space
     real(pfdp), allocatable :: tmp_ter_Rloc(:,:)       ! used for operator in real space
     complex(pfdp), allocatable :: dtempdt_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dVsT_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dxidt_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dVsXi_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dpsidt_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dVsOm_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dVsOm_Mloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: dVsT_Mloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: buo_Mloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: tmphat_Mloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: tmphat_bis_Mloc(:,:) ! used for operator in spectral space
     complex(pfdp), allocatable :: tmphat_Rloc(:,:)     ! used for operator in spectral space
     complex(pfdp), allocatable :: tmphat_bis_Rloc(:,:) ! used for operator in spectral space
     real(pfdp), allocatable :: dtr_Rloc(:)  ! used for operator in real space
     real(pfdp), allocatable :: dth_Rloc(:)  ! used for operator in real space
     real(pfdp), allocatable :: r(:)  ! used for operator in real space
     real(pfdp), allocatable :: or1(:)  ! used for operator in real space
     real(pfdp), allocatable :: or2(:)  ! used for operator in real space
!      class(type_rscheme), pointer :: rscheme
     real(pfdp), allocatable :: tMat_temp(:,:,:)
     integer, allocatable :: tPivot_temp(:,:)
     real(pfdp), allocatable :: tMat_fac(:,:)
     complex(pfdp), allocatable :: rhs_mat(:)
     complex(pfdp), allocatable :: rhs_bounds(:,:)
   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves
     procedure :: initialize  !  Bypasses base sweeper initialize
     procedure :: destroy     !  Bypasses base sweeper destroy

  end type ad_sweeper_t


contains

  !>  Helper function to return sweeper pointer
  function as_ad_sweeper(sweeper) result(r)
    class(pf_sweeper_t), intent(inout), target :: sweeper
    class(ad_sweeper_t), pointer :: r
    select type(sweeper)
    type is (ad_sweeper_t)
       r => sweeper
    class default
       stop
    end select
  end function as_ad_sweeper

  

  !>  Routine to initialize sweeper (bypasses imex sweeper initialize)
  subroutine initialize(this, pf,level_index)
!     use probin, only:  imex_stat,n_modes_m,n_points_phi,n_points_r!,u_s_0,u_phi_0
    use namelists, only: tag, alph1, alph2, l_newmap,&
                       & n_modes_m,n_points_phi,n_points_r
!     use truncation, only: n_cheb_max, m_max
    use fourier, only: fft, ifft
    use radial_functions, only: r, or1, or2!, dtcond, rgrav,rscheme
    use fields, only: us_Mloc, up_Mloc
    use blocking, only: nMstart, nMstop,nRstart, nRstop
  
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index
!    use truncation, only: n_m_max, n_phi_max, idx2m, m2idx
!     integer             :: nMstart,nMstop,n_m_max,n_phi_max,nRstart,nRstop,n_r_max

 
!     integer     :: i,nx,ny
    integer     :: n_r,n_m

    integer             :: n_m_max,n_phi_max,n_r_max
    real(pfdp) :: error
    
!     integer             :: nMstart,nMstop
!     integer             :: nRstart,nRstop
!     nx=pf%levels(level_index)%lev_shape(1)
!     ny=pf%levels(level_index)%lev_shape(2)
!     print*,"nm=",nx,"nr=",ny
     !  Call the imex sweeper initialize
!     nMstart   = 1  ! nMstart_lvl(level_index)
!     nMstop    = n_modes_m(level_index) !nMstop_lvl(level_index)
    n_m_max   = n_modes_m(level_index)
    
    n_phi_max = n_points_phi(level_index)

    n_r_max   = n_points_r(level_index)
!     nRstart = 1   ! nRstart_lvl(level_index)
!     nRstop   = n_r_max !nRstop_lvl(level_index)
    
        call this%imex_initialize(pf,level_index)    


    this%implicit=.TRUE.
    this%explicit=.TRUE.


    allocate(this%u_s_m_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%u_phi_m_Mloc (nMstart:nMstop,n_r_max))
    allocate(this%theta_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%psi_Mloc     (nMstart:nMstop,n_r_max))
    allocate(this%omega_Mloc     (nMstart:nMstop,n_r_max))
!     allocate(this%uphi0_Mloc   (1,n_r_max))
    
    allocate(this%u_s_m_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%u_phi_m_Rloc (n_m_max,nRstart:nRstop))
    allocate(this%omega_m_Rloc (n_m_max,nRstart:nRstop))
    allocate(this%theta_m_Rloc   (n_m_max,nRstart:nRstop))
    
    allocate(this%u_s_Rloc     (n_phi_max,nRstart:nRstop))
    allocate(this%u_phi_Rloc   (n_phi_max,nRstart:nRstop))
    allocate(this%theta_r_Rloc (n_phi_max,nRstart:nRstop))
    allocate(this%om_Rloc      (n_phi_max,nRstart:nRstop))
    allocate(this%xi_Rloc      (n_phi_max,nRstart:nRstop))

    allocate(this%fvec_omega_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%fvec_theta_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%fvec_u_phi_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%fvec_u_s_Mloc   (nMstart:nMstop,n_r_max))
!     allocate(this%fvec_psi_Mloc   (nMstart:nMstop,n_r_max))
!     allocate(this%fvec_uphi0_Mloc   (1,n_r_max))
    
!     allocate(this%rhsvec_Mloc (nMstart:nMstop,n_r_max))
    allocate(this%rhsvec_omega_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%rhsvec_theta_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%rhsvec_u_phi_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%rhsvec_u_s_Mloc   (nMstart:nMstop,n_r_max))
!     allocate(this%rhsvec_psi_Mloc   (nMstart:nMstop,n_r_max))
!     allocate(this%rhsvec_uphi0_Mloc   (1,n_r_max))

!         this%u_phi_m_Mloc = u_phi
!        this%u_s_m_Mloc   = u_s


    this%u_s_m_Mloc(:,:)    = us_Mloc(:,:)
    this%u_phi_m_Mloc(:,:)  = up_Mloc(:,:)
    this%u_s_m_Rloc   = 0.0
    this%u_phi_m_Rloc = 0.0
    this%theta_Mloc   = 0.0
    this%omega_Mloc   = 0.0
!     this%theta_Rloc   = 0.0
    this%u_s_Rloc     = 0.0
    this%u_phi_Rloc   = 0.0
    this%theta_r_Rloc = 0.0

    this%fvec_omega_Mloc    = 0.0
    this%fvec_u_phi_Mloc = 0.0_pfdp
    this%fvec_theta_Mloc    = 0.0
    this%rhsvec_u_s_Mloc = 0.0_pfdp
!     this%rhsvec_Mloc  = 0.0

!     this%fvec_psi_Mloc   = 0.0_pfdp
!     this%fvec_uphi0_Mloc = 0.0_pfdp
    
    allocate(this%dtempdt_Rloc (n_m_max,nRstart:nRstop))
    allocate(this%dVsT_Rloc    (n_m_max,nRstart:nRstop))
    allocate(this%dxidt_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dVsXi_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dpsidt_Rloc  (n_m_max,nRstart:nRstop))
    allocate(this%dVsOm_Rloc   (n_m_max,nRstart:nRstop))
    this%dtempdt_Rloc = 0.0
    this%dVsT_Rloc    = 0.0
    this%dxidt_Rloc   = 0.0
    this%dVsXi_Rloc   = 0.0
    this%dpsidt_Rloc  = 0.0
    this%dVsOm_Rloc   = 0.0
    allocate(this%dVsOm_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%dVsT_Mloc    (nMstart:nMstop,n_r_max))
    allocate(this%buo_Mloc     (nMstart:nMstop,n_r_max))
    this%dVsOm_Mloc   = 0.0
    this%dVsT_Mloc    = 0.0
    this%buo_Mloc     = 0.0
    
    call transp_m2r(m2r_fields, this%u_s_m_Mloc, this%u_s_m_Rloc)

    call transp_m2r(m2r_fields, this%u_phi_m_Mloc, this%u_phi_m_Rloc)
       
       do n_r=nRstart,nRstop
         call ifft(this%u_s_m_Rloc(:,n_r), this%u_s_Rloc(:,n_r))
         call ifft(this%u_phi_m_Rloc(:,n_r), this%u_phi_Rloc(:,n_r))
      enddo

    allocate(this%tmp_Rloc(n_phi_max,nRstart:nRstop))
    allocate(this%tmp_bis_Rloc(n_phi_max,nRstart:nRstop))
    allocate(this%tmp_ter_Rloc(n_phi_max,nRstart:nRstop))
    allocate(this%tmphat_Rloc(n_m_max,nRstart:nRstop))
    allocate(this%tmphat_bis_Rloc(n_m_max,nRstart:nRstop))
    allocate(this%tmphat_Mloc(nMstart:nMstop,n_r_max))
    allocate(this%tmphat_bis_Mloc(nMstart:nMstop,n_r_max))
    this%tmp_Rloc        = 0.0
    this%tmp_bis_Rloc    = 0.0
    this%tmp_ter_Rloc    = 0.0
    this%tmphat_Rloc     = 0.0
    this%tmphat_bis_Rloc = 0.0
    this%tmphat_Mloc     = 0.0
    this%tmphat_bis_Mloc = 0.0

    
    allocate(this%dtr_Rloc(nRstart:nRstop))
    allocate(this%dth_Rloc(nRstart:nRstop))    
    
    allocate(this%r(n_r_max))
    allocate(this%or1(n_r_max))
    allocate(this%or2(n_r_max))
    this%r   = r
    this%or1 = or1
    this%or2 = or2
    
!     allocate(this%rhs_mat(n_r_max))
!     allocate(this%rhs_bounds(nMstart:nMstop,n_r_max))
!     allocate(this%tMat_fac(n_r_max,nMstart:nMstop))
!     allocate(this%tPivot_temp(n_r_max,nMstart:nMstop))
!     allocate(this%tMat_temp(n_r_max,n_r_max,nMstart:nMstop))
!     this%rhs_mat     = 0.0
!     this%rhs_bounds  = 0.0
!     this%tMat_fac    = 0.0
!     this%tPivot_temp = 0
!     this%tMat_temp   = 0.0
!     allocate ( type_cheb :: this%rscheme )

     
     !      real(pfdp), allocatable :: rhs_mat(:)
! 
!     n_in = n_cheb_max
!     if ( l_newmap ) then
!     call this%rscheme%initialize(n_points_r(level_index),n_cheb_max,1,l_cheb_coll)
!     else
!     call this%rscheme%initialize(n_points_r(level_index),n_cheb_max,0,l_cheb_coll)
!     end if
! 
!     call rscheme%initialize(n_points_r(level_index),n_cheb_max,n_in_2,l_cheb_coll)

  end subroutine initialize
  

  !>  Destroy sweeper (bypasses base sweeper destroy)
  subroutine destroy(this,pf,level_index)
    class(ad_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),  target, intent(inout) :: pf
    integer,              intent(in)    :: level_index

    !>  Call base sweeper destroy
    call this%imex_destroy(pf,level_index)


    deallocate(this%u_s_m_Mloc   )
    deallocate(this%u_phi_m_Mloc )
    deallocate(this%theta_Mloc   )

    deallocate(this%u_s_m_Rloc   )
    deallocate(this%u_phi_m_Rloc )
    deallocate(this%omega_m_Rloc )
    deallocate(this%theta_m_Rloc )
    
    deallocate(this%u_s_Rloc     )
    deallocate(this%u_phi_Rloc   )
    deallocate(this%theta_r_Rloc )
    deallocate(this%om_Rloc )
    deallocate(this%xi_Rloc      )

    deallocate(this%fvec_omega_Mloc   )
    deallocate(this%fvec_theta_Mloc   )
!     deallocate(this%rhsvec_Mloc )


    deallocate(this%dtempdt_Rloc )
    deallocate(this%dVsT_Rloc    )
    deallocate(this%dxidt_Rloc   )
    deallocate(this%dVsXi_Rloc   )
    deallocate(this%dpsidt_Rloc  )
    deallocate(this%dVsOm_Rloc   )

    deallocate(this%dVsT_Mloc)
    deallocate(this%buo_Mloc)

    deallocate(this%tmp_Rloc)
    deallocate(this%tmp_bis_Rloc)
    deallocate(this%tmp_ter_Rloc)
    deallocate(this%tmphat_Rloc)
    deallocate(this%tmphat_bis_Rloc)
    deallocate(this%tmphat_Mloc)
    deallocate(this%tmphat_bis_Mloc)
   
    deallocate(this%dtr_Rloc)
    deallocate(this%dth_Rloc)    
    
    deallocate(this%r)
    deallocate(this%or1)
    deallocate(this%or2)


  end subroutine destroy
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
!     use probin, only:  imex_stat ,nu, v

    use namelists, only:  imex_stat , n_modes_m,n_points_phi,n_points_r,&
    &TdiffFac, amp_t,amp_u, r_cmb, r_icb,l_non_rot,l_vort!
    use radial_functions, only: rscheme, r, or1, or2!, dtcond, rgrav
    use radial_der, only: get_ddr, get_dr
    use fourier, only: fft, ifft
    use dct_fftw 
    use truncation, only: idx2m!, m2idx
!     use horizontal, only: hdif_T
    use constants, only: zero, one, two, three, ci, pi, half
    use blocking, only: nMstart, nMstop,nRstart, nRstop
!     use fields, only: us_Mloc,up_Mloc
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece  !  Which piece to solve for
    
    complex(pfdp),      pointer :: theta(:,:),omega(:,:)
    complex(pfdp),      pointer :: u_phi(:,:),u_s(:,:)
    complex(pfdp),      pointer :: fvec_theta(:,:),fvec_omega(:,:)
    complex(pfdp),      pointer :: fvec_u_phi(:,:),fvec_u_s(:,:)
!     integer :: nMstart, nMstop
!     integer :: n_r_max
!     integer :: n_phi_max
    integer :: n_r,n_m,n_phi,m
    integer :: multilevel=1
    real(pfdp) ::h2,  dm2
    real(pfdp) :: phi
!     real(pfdp) :: error
!     integer             :: n_m_max,n_phi_max,n_r_max
    integer             :: n_phi_max,n_r_max
    logical :: dummylogical=.true.
!     integer             :: nMstart,nMstop
!     integer             :: nRstart,nRstop
!     nMstart   = 1  ! nMstart_lvl(level_index)
!     nMstop    = n_modes_m(level_index) !nMstop_lvl(level_index)
   
    n_r_max   = n_points_r(level_index)
    n_phi_max = n_points_phi(level_index)
    
!     nRstart = 1   ! nRstart_lvl(level_index)
!     nRstop   = n_r_max !nRstop_lvl(level_index)
    !  Grab the arrays from the encap
!     theta_Mloc  => get_array2d(y)

       omega  => get_array2d(y,1)
!        psi    => get_array2d(y,2)
       u_phi  => get_array2d(y,2)
       u_s    => get_array2d(y,3)
       theta  => get_array2d(y,4)
!        uphi0  => get_array2d(y,4)

! 
!     select type (y)
!     class is (pizza_zndsysarray_t)
!               print*,"size y%arr_shape",y%arr_shape
!               print*,"size y%flatarray",SIZE(y%flatarray)
!               
!     end select           
  
       fvec_omega => get_array2d(f,1)
!        fvec_psi   => get_array2d(f,2)
       fvec_u_phi => get_array2d(f,2)
       fvec_u_s   => get_array2d(f,3)
       fvec_theta => get_array2d(f,4)
!     fvec_uphi0 => get_array2d(f,4)


    this%fvec_omega_Mloc = 0.0_pfdp
    this%fvec_u_phi_Mloc = 0.0_pfdp
    this%fvec_u_s_Mloc   = 0.0_pfdp
    this%fvec_theta_Mloc = 0.0_pfdp


    this%theta_Mloc   = theta 
    this%omega_Mloc   = omega
    this%u_phi_m_Mloc = u_phi
    this%u_s_m_Mloc   = u_s

    select case (piece)
    case (1)  ! Explicit piece

 



    call transp_m2r(m2r_fields, this%omega_Mloc, this%omega_m_Rloc)
    call transp_m2r(m2r_fields, this%theta_Mloc, this%theta_m_Rloc)
    call transp_m2r(m2r_fields, this%u_s_m_Mloc, this%u_s_m_Rloc)
    call transp_m2r(m2r_fields, this%u_phi_m_Mloc, this%u_phi_m_Rloc)

    call radial_loop(this%u_s_m_Rloc, this%u_phi_m_Rloc, this%omega_m_Rloc, this%theta_m_Rloc, this%xi_m_Rloc,   &
                 &           this%dtempdt_Rloc, this%dVsT_Rloc, this%dxidt_Rloc, this%dVsXi_Rloc, &
                 &           this%dpsidt_Rloc, this%dVsOm_Rloc, this%dtr_Rloc, this%dth_Rloc)

!      if ( l_heat ) then
        call transp_r2m(r2m_fields, this%dtempdt_Rloc, &
                       &          this%fvec_theta_Mloc)

        call transp_r2m(r2m_fields, this%dVsT_Rloc, this%dVsT_Mloc)
!     end if
    
!     if ( l_chem ) then
!        call transp_r2m(r2m_fields, dxidt_Rloc, &
!                       &          dxidt%expl(:,:,tscheme%istage))
!        call transp_r2m(r2m_fields, dVsXi_Rloc, dVsXi_Mloc)
!     end if
    if (l_vort) then

    call transp_r2m(r2m_fields,  this%dpsidt_Rloc, &
         &          this%fvec_omega_Mloc)
    call transp_r2m(r2m_fields, this%dVsOm_Rloc, this%dVsOm_Mloc)
    endif    
  
    !--------------------
    !-- Finish assembling the explicit terms
    !--------------------
!     runStart = MPI_Wtime()


    call finish_exp_temp_coll(this%theta_Mloc, this%u_s_m_Mloc,    &
         &                    this%dVsT_Mloc, this%buo_Mloc,    &
         &                    this%fvec_theta_Mloc)
    if (l_vort) then
       call finish_exp_psi_coll_smat(this%u_s_m_Mloc, this%dVsOm_Mloc, this%buo_Mloc,  &
            &                        this%fvec_omega_Mloc)     
    endif     

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
            
               this%fvec_u_phi_Mloc(n_m,n_r)   = this%fvec_omega_Mloc(n_m,n_r)

            end if
!             if ( m == 3 ) then
!             print*,"this%fvec_omega_Mloc(m=3,n_r=",n_r,")=",this%fvec_omega_Mloc(n_m,n_r)
!             endif
         end do
      end do 
      
    case (2)  ! Implicit piece
          call get_temp_rhs_imp_coll(this%theta_Mloc, this%tmphat_Mloc, this%tmphat_bis_Mloc, &
              &                     this%fvec_theta_Mloc, .true.)

          if (l_vort) then
             call get_psi_rhs_imp_coll_smat(this%u_s_m_Mloc, this%u_phi_m_Mloc, this%omega_Mloc, this%tmphat_Mloc,    &
                  &                         this%tmphat_bis_Mloc,         &
                  &                         this%fvec_omega_Mloc,l_calc_lin_rhs=.true.)
          endif

!       do n_r=1,n_r_max
!          do n_m=nMstart,nMstop
!             m = idx2m(n_m)
!             if ( m == 3 ) then
! !             print*,"this%fvec_omega_Mloc(m=3,n_r=",n_r,")=",this%fvec_omega_Mloc(n_m,n_r)
! !             print*,"this%fvec_omega_Mloc(m=3,n_r=",n_r,")=",this%fvec_omega_Mloc(n_m,n_r)
! !             print*,"this%fvec_omega_Mloc(m=3,n_r=",n_r,")=",this%fvec_omega_Mloc(n_m,n_r)
!             endif
!          end do
!       end do 
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
            
               this%fvec_u_phi_Mloc(n_m,n_r)   = this%fvec_omega_Mloc(n_m,n_r)
!                print*,"this%fvec_omega_Mloc(m=3,n_r=",n_r,")=",this%fvec_omega_Mloc(n_m,n_r)
               
!                this%fvec_omega_Mloc(n_m,n_r) = 0.0_8
!                this%fvec_omega_Mloc(n_m,n_r) = 0.0_8
!                this%fvec_omega_Mloc(n_m,n_r) = 0.0_8
               
            end if
         end do
      end do 

    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       return
    end select

!     fvec =this%fvec_Mloc

   


    fvec_theta = this%fvec_theta_Mloc
    fvec_omega = this%fvec_omega_Mloc
    fvec_u_phi = this%fvec_u_phi_Mloc
    fvec_u_s   = this%fvec_u_s_Mloc


  end subroutine f_eval

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
  use fieldsLast, only: dTdt
!     use probin, only:  imex_stat ,nu,v
    use namelists, only:  imex_stat , n_modes_m,n_points_phi,n_points_r,&
    &TdiffFac, amp_t,amp_u, r_cmb, r_icb,l_non_rot,nnodes,l_vort!    use radial_functions, only: rscheme, or1, or2!, dtcond, rgrav
    use truncation, only: idx2m!, m2idx
    use horizontal, only: bott_Mloc, topt_Mloc
    use blocking, only: nMstart, nMstop!,nRstart, nRstop
    use radial_functions, only: rscheme, r, or1, or2!, dtcond, rgrav
    use constants, only: zero, one, two, three, ci, pi, half
    use radial_der, only: get_ddr, get_dr

    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)

    complex(pfdp),         pointer :: omega(:,:), u_phi(:,:), theta(:,:), u_s(:,:)
    complex(pfdp),         pointer :: fvec_omega(:,:), fvec_u_phi(:,:),fvec_u_s(:,:), fvec_theta(:,:)
    complex(pfdp),         pointer :: rhs_omega(:,:), rhs_u_phi(:,:),rhs_u_s(:,:), rhs_theta(:,:)
    integer :: m,n_m,n_r,n_r_out
    real(pfdp) ::h2,  dm2
    integer             :: n_r_max

!
!     integer :: i_substep

!     integer::nMstart, nMstop
! 
!     nMstart=1
!     nMstop =n_modes_m(level_index)
    
    n_r_max   = n_points_r(level_index)
   
    if (piece == 2) then

       omega  => get_array2d(y,1)
       u_phi  => get_array2d(y,2)
       u_s    => get_array2d(y,3)
       theta  => get_array2d(y,4)

!     print*,"size y",SIZE(y%flatarray),"size f   =",SIZE(f%flatarray)
!     select type (y)
!     class is (pizza_zndsysarray_t)
!               print*,"size y%arr_shape",y%arr_shape
!               print*,"size y%flatarray",SIZE(y%flatarray)
!               
!     end select           
!               
!     print*,"size omega",SIZE(omega),"this%omega_Mloc   =",SIZE(this%omega_Mloc)
!     print*,"size u_phi",SIZE(u_phi),"this%u_phi_m_Mloc =",SIZE(this%u_phi_m_Mloc)
!     print*,"size u_s"  ,SIZE(u_s)  ,"this%u_s_m_Mloc   =",SIZE(this%u_s_m_Mloc)
!     print*,"size theta",SIZE(theta),"this%theta_Mloc   =",SIZE(this%theta_Mloc)
       
       fvec_omega => get_array2d(f,1)
       fvec_u_phi => get_array2d(f,2)
       fvec_u_s   => get_array2d(f,3)
       fvec_theta => get_array2d(f,4)
       
       rhs_omega => get_array2d(rhs,1)
       rhs_u_phi => get_array2d(rhs,2)
       rhs_u_s   => get_array2d(rhs,3)
       rhs_theta => get_array2d(rhs,4)
       

       this%omega_Mloc   = omega
       this%u_phi_m_Mloc = u_phi
       this%u_s_m_Mloc   = u_s
       this%theta_Mloc   = theta 
       

       this%rhsvec_theta_Mloc = rhs_theta 
       this%rhsvec_u_phi_Mloc = rhs_u_phi 
       this%rhsvec_u_s_Mloc   = rhs_u_s 
       this%rhsvec_omega_Mloc = rhs_omega 

       this%fvec_omega_Mloc   = fvec_omega 
       this%fvec_u_phi_Mloc   = fvec_u_phi 
       this%fvec_u_s_Mloc     = fvec_u_s 
       this%fvec_theta_Mloc   = fvec_theta 
   
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
            
!                this%fvec_uphi0_Mloc(n_m,n_r)=this%fvec_omega_Mloc(n_m,n_r)
               
!                this%fvec_omega_Mloc(n_m,n_r)   = this%fvec_u_phi_Mloc(n_m,n_r)
               this%rhsvec_omega_Mloc(n_m,n_r) = this%rhsvec_u_phi_Mloc(n_m,n_r)
               
            end if
         end do
      end do    

       call update_temp_co(this%theta_Mloc, this%tmphat_Mloc, this%tmphat_bis_Mloc, dTdt, &
            &              lMat, l_log_next,dtq=dtq,work_Mloc_pfasst=this%rhsvec_theta_Mloc,int_mat=i_substep )
       if (l_vort) then

          call update_om_coll_smat(this%psi_Mloc, this%omega_Mloc, this%tmphat_Mloc, this%u_s_m_Mloc, this%u_phi_m_Mloc, &
                 &                   this%buo_Mloc, &
                 &                   lMat, dtq=dtq,work_Mloc_pfasst=this%rhsvec_omega_Mloc,int_mat=i_substep )
       endif
       !  The function is easy to derive
       this%fvec_theta_Mloc = (this%theta_Mloc - this%rhsvec_theta_Mloc) / dtq
       this%fvec_u_s_Mloc = 0.0
!        this%fvec_psi_Mloc = 0.0
       
       this%fvec_omega_Mloc = (this%omega_Mloc - this%rhsvec_omega_Mloc) / dtq

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            if ( m == 0 ) then
            
!                this%uphi0_Mloc(1,n_r) = this%u_phi_m_Mloc(n_m,n_r)
               this%fvec_u_phi_Mloc(n_m,n_r)   = (this%u_phi_m_Mloc(n_m,n_r) &
                                               & - this%rhsvec_u_phi_Mloc(n_m,n_r)) / dtq
!                this%rhsvec_u_phi_Mloc(n_m,n_r) = this%rhsvec_omega_Mloc(n_m,n_r)
               
               this%fvec_omega_Mloc(n_m,n_r)   = 0.0_8
!                this%rhsvec_omega_Mloc(n_m,n_r) = 0.0_8
               
               
            else
               this%fvec_u_phi_Mloc(n_m,n_r)   = 0.0
!                this%rhsvec_u_phi_Mloc(n_m,n_r) = 0.0
               
            endif
         end do
      end do    

!        this%fvec_uphi0_Mloc = (this%uphi0_Mloc - this%rhsvec_uphi0_Mloc) / dtq


  
       fvec_theta = this%fvec_theta_Mloc
       fvec_omega = this%fvec_omega_Mloc
       fvec_u_s   = this%fvec_u_s_Mloc    
       fvec_u_phi = this%fvec_u_phi_Mloc    
!        fvec_psi  = this%fvec_psi_Mloc    
!        fvec_uphi0 =  this%fvec_uphi0_Mloc
         
       theta = this%theta_Mloc 
       omega = this%omega_Mloc
       u_s   = this%u_s_m_Mloc 
       u_phi = this%u_phi_m_Mloc 
!        psi   = this%psi_Mloc 
!        uphi0 = this%uphi0_Mloc
   
    else
       print *,'Bad piece in f_comp ',piece
    end if
    if (i_substep.LT.(nnodes)) then
       i_substep=i_substep+1
    end if
    if (i_substep.EQ.(nnodes)) then
       i_substep = 1
       lMat=.false.
    end if
 
  end subroutine f_comp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  subroutine initial_condition(y_0)
    use fields, only: temp_Mloc,om_Mloc,us_Mloc,up_Mloc
!     use namelists, only:  imex_stat , n_modes_m,n_points_phi,n_points_r,r_icb
!     use radial_functions, only: rscheme, r, or1, or2!, dtcond, rgrav
!        use constants, only: zero, one, two, three, ci, pi, half
!     use blocking, only: nMstart, nMstop,nRstart, nRstop
!     use namelists, only:  imex_stat ,n_modes_m,n_points_phi,n_points_r,TdiffFac
!     use radial_functions, only: r, or1, or2!, dtcond, rgrav

!     type(pf_zndarray_t), intent(inout) :: y_0
    type(pizza_zndsysarray_t), intent(inout) :: y_0
!     class(pf_encap_t),   intent(in   ) :: y_0
!     integer,             intent(in   ) :: level_index !  Which level this is
    complex(pfdp), pointer :: theta_Mloc(:,:)
    complex(pfdp), pointer :: omega_Mloc(:,:)
    complex(pfdp), pointer :: psi_Mloc(:,:)
    complex(pfdp), pointer :: u_phi_Mloc(:,:)
    complex(pfdp), pointer :: u_s_Mloc(:,:)
    complex(pfdp), pointer :: uphi0_Mloc(:,:)
!     integer :: n_r,n_m
!     integer :: idx

    omega_Mloc => get_array2d(y_0,1)
!     psi_Mloc   => get_array2d(y_0,2)
    u_phi_Mloc => get_array2d(y_0,2)
    u_s_Mloc   => get_array2d(y_0,3)
    theta_Mloc => get_array2d(y_0,4)
!     uphi0_Mloc => get_array2d(y_0,4)
    
    omega_Mloc = om_Mloc
    u_phi_Mloc = up_Mloc
    u_s_Mloc   = us_Mloc
    theta_Mloc = temp_Mloc
       


  end subroutine initial_condition

  !> Routine to return the exact solution
!   subroutine exact(t, yex)
! !     use probin, only: nprob,nu, a,b, kfreq, t00
!     real(pfdp), intent(in)  :: t
!     
!     complex(pfdp), intent(out) :: yex(:,:)
! 
! !     integer    :: nx,ny!, i, j
! !     real(pfdp) :: tol, x,y, r2, omega
! ! 
! !     nx = size(yex,1)
! !     ny = size(yex,2)
! 
!     !  Using sin wave initial condition
! !     if (nprob .eq. 1) then
! ! !!$       omega = two_pi*kfreq
! ! !!$       do j = 1, ny
! ! !!$          y = dble(j-1)/dble(ny)-0.5_pfdp - t*b 
! ! !!$          do i = 1, nx
! ! !!$             x = dble(i-1)/dble(nx)-0.5_pfdp - t*a 
! ! !!$             yex(i,j) = sin(omega*x)*sin(omega*y)*exp(-2.0_pfdp*omega*omega*nu*t)
! ! !!$          end do
! ! !!$       end do
! ! !        call exact_ad_cos(t,yex,nu,[a,b],[kfreq,kfreq],[1.0_pfdp,1.0_pfdp])
! !     else
! !        do j = 1, ny
! !           y = (1.0_8*(j-1))/(1.0_8*(ny))-0.5_pfdp - t*b 
! !           do i = 1, nx
! !              x = (1.0_8*(i-1))/(1.0_8*(nx))-0.5_pfdp - t*a
! !              r2=x*x+y*y
! !              yex(i,j) = t00/(t00+t)*exp(-r2/(4.0*nu*(t00+t)))             
! !           end do
! !        end do
! !     endif
!   end subroutine exact
! 
! #ifdef WITH_PRECOND_S
!    subroutine get_tempMat_pfasst(this,dtq,rscheme, m, tMat, tPivot, tMat_fac,n_r_max)
! #else
!    subroutine get_tempMat_pfasst(this,dtq,rscheme, m, tMat, tPivot,n_r_max)
! #endif
!        use truncation, only: m2idx
!       use namelists, only: kbott, ktopt, TdiffFac
! !       use horizontal, only: hdif_T
! !     use radial_functions, only: or1, or2!, dtcond, rgrav
! 
!       !-- Input variables
! !       class(type_tscheme), intent(in) :: tscheme        ! time step
!       integer,             intent(in) :: m
!       integer :: nR_out, nR, info, n_m,n_r_max
!       class(ad_sweeper_t), intent(inout) :: this
! 
!       !-- Output variables
!       real(pfdp), intent(out) :: tMat(n_r_max,n_r_max)
!       integer,  intent(out) :: tPivot(n_r_max)
! #ifdef WITH_PRECOND_S
!       real(pfdp),intent(out) :: tMat_fac(n_r_max)
! #endif
!      class(type_rscheme), pointer :: rscheme
! 
!       !-- Local variables
!       real(pfdp) :: dm2
!       real(pfdp) :: dtq
! 
!       dm2 = real(m,pfdp)*real(m,pfdp)
!       n_m = m2idx(m)
! 
!       !----- Boundary coditions:
!       do nR_out=1,rscheme%n_max
!          if ( ktopt == 1 ) then
!             tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
!          else
!             tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
!          end if
!          if ( kbott == 1 ) then
!             tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
!          else
!             tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
!          end if
!       end do
! 
!       if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
!          do nR_out=rscheme%n_max+1,n_r_max
!             tMat(1,nR_out)      =0.0_pfdp
!             tMat(n_r_max,nR_out)=0.0_pfdp
!          end do
!       end if
! 
!       !----- Other points:
!       do nR_out=1,n_r_max
!          do nR=2,n_r_max-1
!             tMat(nR,nR_out)= rscheme%rnorm * (                          &
!             &                                 rscheme%rMat(nR,nR_out) - &
!             &                                          dtq*TdiffFac*(   &
!             &                               rscheme%d2rMat(nR,nR_out) + &
!             &          this%or1(nR)*         rscheme%drMat(nR,nR_out) - &
!             &      dm2*this%or2(nR)*          rscheme%rMat(nR,nR_out) ) )
!          end do
!       end do
! 
!       !----- Factor for highest and lowest cheb:
!       do nR=1,n_r_max
!          tMat(nR,1)      =rscheme%boundary_fac*tMat(nR,1)
!          tMat(nR,n_r_max)=rscheme%boundary_fac*tMat(nR,n_r_max)
!       end do
! 
! #ifdef WITH_PRECOND_S
!       ! compute the linesum of each line
!       do nR=1,n_r_max
!          tMat_fac(nR)=one/maxval(abs(tMat(nR,:)))
!       end do
!       ! now divide each line by the linesum to regularize the matrix
!       do nR=1,n_r_max
!          tMat(nR,:) = tMat(nR,:)*tMat_fac(nR)
!       end do
! 
! #endif
! 
!       !----- LU decomposition:
!       call prepare_full_mat(tMat,n_r_max,n_r_max,tPivot,info)
!       if ( info /= 0 ) then
!          call abortRun('Singular matrix tMat!')
!       end if
! 
!    end subroutine get_tempMat_pfasst
   
end module pf_my_sweeper
