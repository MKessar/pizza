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
  use pizza_zdnarray
  use update_temp_coll, only: update_temp_co, get_temp_rhs_imp_coll,finish_exp_temp_coll
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
     complex(pfdp), allocatable :: theta_Mloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: xi_m_Rloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: om_m_Rloc(:,:)      !u_phi velocity in complex space
     complex(pfdp), allocatable :: fvec_Mloc(:,:)       !u_phi velocity in complex space
     complex(pfdp), allocatable :: rhsvec_Mloc(:,:)     !u_phi velocity in complex space
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
! 
! subroutine initit_sweeper_pfasst()
! 
! lMat=.false.
! l_log_next=.false.
! 
! end subroutine initit_sweeper_pfasst


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

!     complex(pfdp), allocatable :: lap(:,:)         ! Lapclacian operators
!     complex(pfdp), allocatable :: ddx(:,:) ! First derivative operators
!     complex(pfdp), allocatable :: ddy(:,:) ! First derivative operators
 
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

! print*,"lMat=",lMat
! print*,"l_log_next=",l_log_next
! lMat=.true.
! l_log_next=.true.
! print*,"lMat=",lMat
! print*,"l_log_next=",l_log_next

    allocate(this%u_s_m_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%u_phi_m_Mloc (nMstart:nMstop,n_r_max))
    allocate(this%theta_Mloc   (nMstart:nMstop,n_r_max))

    allocate(this%u_s_m_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%u_phi_m_Rloc (n_m_max,nRstart:nRstop))
    allocate(this%om_m_Rloc (n_m_max,nRstart:nRstop))
    allocate(this%theta_m_Rloc   (n_m_max,nRstart:nRstop))
    
    allocate(this%u_s_Rloc     (n_phi_max,nRstart:nRstop))
    allocate(this%u_phi_Rloc   (n_phi_max,nRstart:nRstop))
    allocate(this%theta_r_Rloc (n_phi_max,nRstart:nRstop))
    allocate(this%om_Rloc      (n_phi_max,nRstart:nRstop))
    allocate(this%xi_Rloc      (n_phi_max,nRstart:nRstop))

    allocate(this%fvec_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%rhsvec_Mloc (nMstart:nMstop,n_r_max))

    this%u_s_m_Mloc(:,:)    = us_Mloc(:,:)
    this%u_phi_m_Mloc(:,:)  = up_Mloc(:,:)
    this%u_s_m_Rloc   = 0.0
    this%u_phi_m_Rloc = 0.0
    this%theta_Mloc   = 0.0
!     this%theta_Rloc   = 0.0
    this%u_s_Rloc     = 0.0
    this%u_phi_Rloc   = 0.0
    this%theta_r_Rloc = 0.0
    this%fvec_Mloc    = 0.0
    this%rhsvec_Mloc  = 0.0
    
    allocate(this%dtempdt_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dVsT_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dxidt_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dVsXi_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dpsidt_Rloc   (n_m_max,nRstart:nRstop))
    allocate(this%dVsOm_Rloc   (n_m_max,nRstart:nRstop))
    this%dtempdt_Rloc = 0.0
    this%dVsT_Rloc    = 0.0
    this%dxidt_Rloc   = 0.0
    this%dVsXi_Rloc   = 0.0
    this%dpsidt_Rloc  = 0.0
    this%dVsOm_Rloc   = 0.0
    allocate(this%dVsT_Mloc   (nMstart:nMstop,n_r_max))
    allocate(this%buo_Mloc   (nMstart:nMstop,n_r_max))
    this%dVsT_Mloc   = 0.0
    this%buo_Mloc   = 0.0
    
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
    deallocate(this%om_m_Rloc    )
    deallocate(this%theta_m_Rloc )
    
    deallocate(this%u_s_Rloc     )
    deallocate(this%u_phi_Rloc   )
    deallocate(this%theta_r_Rloc )
    deallocate(this%om_Rloc      )
    deallocate(this%xi_Rloc      )

    deallocate(this%fvec_Mloc   )
    deallocate(this%rhsvec_Mloc )

    
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
    &TdiffFac, amp_t,amp_u, r_cmb, r_icb!
    use radial_functions, only: rscheme!, r, or1, or2!, dtcond, rgrav
    use radial_der, only: get_ddr, get_dr
    use fourier, only: fft, ifft
    use dct_fftw 
    use truncation, only: idx2m!, m2idx
!     use horizontal, only: hdif_T
    use constants, only: zero, one, two, three, ci, pi, half
    use blocking, only: nMstart, nMstop,nRstart, nRstop
    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece  !  Which piece to solve for
    
    complex(pfdp),      pointer :: theta(:,:), fvec(:,:)
!     integer :: nMstart, nMstop
!     integer :: n_r_max
!     integer :: n_phi_max
    integer :: n_r,n_m,n_phi,m
    integer :: multilevel=1
    real(pfdp) :: dm2
    real(pfdp) :: phi
!     real(pfdp) :: error
!     integer             :: n_m_max,n_phi_max,n_r_max
    integer             :: n_phi_max,n_r_max
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

    theta  => get_array2d(y)
    fvec => get_array2d(f)
    this%fvec_Mloc=0.0_pfdp
 
 
    this%theta_Mloc = theta 

    select case (piece)
    case (1)  ! Explicit piece

    call transp_m2r(m2r_fields, this%theta_Mloc, this%theta_m_Rloc)

    call radial_loop(this%u_s_m_Rloc, this%u_phi_m_Rloc, this%om_m_Rloc, this%theta_m_Rloc, this%xi_m_Rloc,   &
                 &           this%dtempdt_Rloc, this%dVsT_Rloc, this%dxidt_Rloc, this%dVsXi_Rloc, &
                 &           this%dpsidt_Rloc, this%dVsOm_Rloc, this%dtr_Rloc, this%dth_Rloc)

!      if ( l_heat ) then
        call transp_r2m(r2m_fields, this%dtempdt_Rloc, &
                       &          this%fvec_Mloc)

        call transp_r2m(r2m_fields, this%dVsT_Rloc, this%dVsT_Mloc)
!     end if
    
!     if ( l_chem ) then
!        call transp_r2m(r2m_fields, dxidt_Rloc, &
!                       &          dxidt%expl(:,:,tscheme%istage))
!        call transp_r2m(r2m_fields, dVsXi_Rloc, dVsXi_Mloc)
!     end if
    
!     call transp_r2m(r2m_fields, dpsidt_Rloc, &
!          &          dpsidt%expl(:,:,tscheme%istage))
!     call transp_r2m(r2m_fields, this%dVsOm_Rloc, this%dVsOm_Mloc)
    
    
    !--------------------
    !-- Finish assembling the explicit terms
    !--------------------
!     runStart = MPI_Wtime()
!     call finish_explicit_assembly(this%theta_Mloc, xi_Mloc, psi_Mloc,      &
!     &                        us_Mloc, up_Mloc, om_Mloc,         &
!     &                        this%dVsT_Mloc, dVsXi_Mloc, dVsOm_Mloc, &
!     &                        buo_Mloc, dTdt, dxidt, dpsidt,     &
!     &                        tscheme, vp_bal, vort_bal)
    call finish_exp_temp_coll(this%theta_Mloc, this%u_s_m_Mloc,    &
         &                    this%dVsT_Mloc, this%buo_Mloc,    &
         &                    this%fvec_Mloc)
!                  
!       do n_r=nRstart,nRstop
!          call ifft(this%theta_Rloc(:,n_r), this%theta_r_Rloc(:,n_r))
!          
!          do n_phi=1,n_phi_max
!             this%tmp_Rloc(n_phi,n_r)= this%r(n_r)*this%u_s_Rloc(n_phi,n_r)*this%theta_r_Rloc(n_phi,n_r)
!          enddo
!       enddo
! 
! 
!       do n_r=nRstart,nRstop
!          call fft(this%tmp_Rloc(:,n_r), this%tmphat_Rloc(:,n_r))
!       enddo
! 
!             !!! tranpose r2m
!       call transp_r2m(r2m_fields,  this%tmphat_Rloc, this%tmphat_Mloc)
!     
!       call get_dr(this%tmphat_Mloc, this%fvec_Mloc, nMstart, nMstop, n_r_max, rscheme,multi_level=multilevel)
! 
!       do n_r=1,n_r_max
!          do n_m=nMstart,nMstop
!          this%fvec_Mloc(n_m,n_r)= -this%or1(n_r)*this%fvec_Mloc(n_m,n_r)
!          enddo
!       enddo
!  
!     
!       do n_r=nRstart,nRstop
!          do n_phi=1,n_phi_max
!             this%tmp_Rloc(n_phi,n_r)= this%u_phi_Rloc(n_phi,n_r)*this%theta_r_Rloc(n_phi,n_r)
!          enddo
!       enddo
! 
! 
!       do n_r=nRstart,nRstop
!          call fft(this%tmp_Rloc(:,n_r), this%tmphat_Rloc(:,n_r))
!       enddo
! 
!       !!! tranpose r2m
!       call transp_r2m(r2m_fields,  this%tmphat_Rloc, this%tmphat_Mloc)
! 
!       do n_r=1,n_r_max
!          do n_m=nMstart,nMstop
!             m = idx2m(n_m)
!             this%fvec_Mloc(n_m,n_r) = this%fvec_Mloc(n_m,n_r) &
!                         & - this%or1(n_r)*ci*m*this%tmphat_Mloc(n_m,n_r)
!          enddo
!       enddo

     
         
    case (2)  ! Implicit piece
          call get_temp_rhs_imp_coll(this%theta_Mloc, this%tmphat_Mloc, this%tmphat_bis_Mloc, &
              &                     this%fvec_Mloc, .true.)
              
! !          call get_ddr(theta, this%tmphat,this%tmphat_bis, nMstart, nMstop, n_points_r(level_index), rscheme,multi_level=multilevel)
!          call get_ddr(this%theta_Mloc, this%tmphat_Mloc,this%tmphat_bis_Mloc, nMstart, nMstop, n_r_max, rscheme)
!  
!          do n_r=1,n_r_max
!             do n_m=nMstart,nMstop
!                m = idx2m(n_m)
!                dm2 = real(m,pfdp)*real(m,pfdp)
!                this%fvec_Mloc(n_m,n_r) = TdiffFac*(             this%tmphat_bis_Mloc(n_m,n_r) &
!                                   & + this%or1(n_r)*                this%tmphat_Mloc(n_m,n_r) &
!                                   & - dm2*this%or2(n_r)*           this%theta_Mloc(n_m,n_r) )
!                                
!             enddo
!          enddo
!         print*,"test =",1

    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       return
    end select

    fvec =this%fvec_Mloc
          
  end subroutine f_eval

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
  use fieldsLast, only: dTdt
!     use probin, only:  imex_stat ,nu,v
    use namelists, only:  n_points_r,nnodes
    use radial_functions, only: rscheme!, or1, or2!, dtcond, rgrav
    use truncation, only: idx2m!, m2idx
    use horizontal, only: bott_Mloc, topt_Mloc
    use blocking, only: nMstart, nMstop!,nRstart, nRstop

    class(ad_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)

    complex(pfdp),         pointer :: yvec(:,:), rhsvec(:,:), fvec(:,:)
    integer :: m,n_m,n_r,n_r_out
!     integer :: i_substep

!     integer::nMstart, nMstop
! 
!     nMstart=1
!     nMstop =n_modes_m(level_index)
    
    
    
    if (piece == 2) then
       yvec  => get_array2d(y)
       rhsvec => get_array2d(rhs)
       fvec => get_array2d(f)

       this%fvec_Mloc   = fvec 
       this%theta_Mloc  = yvec 
       this%rhsvec_Mloc = rhsvec 

       call update_temp_co(this%theta_Mloc, this%tmphat_Mloc, this%tmphat_bis_Mloc, dTdt, &
            &              lMat, l_log_next,dtq=dtq,work_Mloc_pfasst=this%rhsvec_Mloc,int_mat=i_substep )
!               
!       do n_m=nMstart, nMstop
! 
!          m = idx2m(n_m)
!      
! #ifdef WITH_PRECOND_S
!             call get_tempMat_pfasst(this,dtq,rscheme, m, this%tMat_temp(:,:,n_m), this%tPivot_temp(:,n_m), &
!                  &           this%tMat_fac(:,n_m),n_points_r(level_index))
! #else
!             call get_tempMat_pfasst(this,dtq,rscheme, m, this%tMat_temp(:,:,n_m), this%tPivot_temp(:,n_m),n_points_r(level_index))
! #endif
! !             lTmat(n_m)=.true.
! 
!          !-- Inhomogeneous B.Cs (if not zero)
!          this%rhs_mat(1)                      =topt_Mloc(n_m)
!          this%rhs_mat(n_points_r(level_index))=bott_Mloc(n_m)
!          do n_r=2,n_points_r(level_index)-1
!             this%rhs_mat(n_r)=this%rhsvec_Mloc(n_m,n_r)
!          end do
!          
! !          this%rhs_bounds(n_m, :)=this%rhs_mat(:)
! 
! #ifdef WITH_PRECOND_S
!          do n_r=1,n_points_r(level_index)
!             this%rhs_mat(n_r) = this%tMat_fac(n_r,n_m)*this%rhs_mat(n_r)
!          end do
! #endif
!               
!          call solve_full_mat(this%tMat_temp(:,:,n_m), n_points_r(level_index), n_points_r(level_index), this%tPivot_temp(:, n_m), &
!               &              this%rhs_mat(:))
! 
!          do n_r_out=1,rscheme%n_max
! !             temp_Mloc(n_m, n_r_out)=this%rhs(n_r_out)
!             this%theta_Mloc(n_m, n_r_out)=this%rhs_mat(n_r_out)
!          end do
!       end do
!       
!     !-- set cheb modes > rscheme%n_max to zero (dealiazing)
!       do n_r_out=rscheme%n_max+1,n_points_r(level_index)
!          do n_m=nMstart,nMstop
!             this%theta_Mloc(n_m,n_r_out)=zero
!          end do
!       end do
! 
!       !-- Bring temperature back to physical space
!       call rscheme%costf1(this%theta_Mloc, nMstart, nMstop, n_points_r(level_index))

      !-- Assemble buoyancy in case this is treated implicitly
!       if ( l_buo_imp ) then
!          call tscheme%assemble_implicit_buo(buo_Mloc, temp_Mloc, dTdt,      &
!               &                             BuoFac, rgrav, nMstart, nMstop, &
!               &                             n_r_max, .true.)
!       end if

      !-- Roll the arrays before filling again the first block
!       call tscheme%rotate_imex(dTdt, nMstart, nMstop, n_r_max)

      !-- In case log is needed on the next iteration, recalculate dT/dr
!       if ( l_log_next ) then
!          call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
!       end if

       !  The function is easy to derive
       this%fvec_Mloc = (this%theta_Mloc - this%rhsvec_Mloc) / dtq
  
       fvec =this%fvec_Mloc 
       yvec =this%theta_Mloc   
    else
       print *,'Bad piece in f_comp ',piece
    end if
    i_substep=i_substep+1
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
    use fields, only: temp_Mloc
!     use namelists, only:  imex_stat , n_modes_m,n_points_phi,n_points_r,r_icb
!     use radial_functions, only: rscheme, r, or1, or2!, dtcond, rgrav
!        use constants, only: zero, one, two, three, ci, pi, half
!     use blocking, only: nMstart, nMstop,nRstart, nRstop
!     use namelists, only:  imex_stat ,n_modes_m,n_points_phi,n_points_r,TdiffFac
!     use radial_functions, only: r, or1, or2!, dtcond, rgrav

!     type(pf_zndarray_t), intent(inout) :: y_0
    type(pizza_zndarray_t), intent(inout) :: y_0
!     class(pf_encap_t),   intent(in   ) :: y_0
!     integer,             intent(in   ) :: level_index !  Which level this is
    complex(pfdp), pointer :: yvec(:,:)
!     integer :: n_r,n_m
!     integer :: idx
    
!  print*,"Barrier a init condition Key_Pizza =",Key_Pizza 
!     call mpi_barrier(Comm_Pizza, ierror)
! 
!  print*,"Barrier b init condition Key_Pizza =",Key_Pizza     
    yvec => get_array2d(y_0)  
    print*,"shape yvec", shape(yvec)
    print*,"shape temp_Mloc", shape(temp_Mloc)
!      yvec=0.0!temp_Mloc
     yvec=temp_Mloc
!  print*,"Barrier c init condition Key_Pizza =",Key_Pizza 
!     call mpi_barrier(Comm_Pizza, ierror)
! 
!  print*,"Barrier d init condition Key_Pizza =",Key_Pizza 
       
     
     
         
     
!       do n_r=1,n_points_r(1)
! !          E_temp_radial(n_r) = 0.0_pfdp
! !          do n_m=1,n_modes_m(1)
! !             m = idx2m(n_m)
!              yvec(4,n_r) = 100.0_8*exp(-r(n_r)/r_icb)*(tan(pi*r_icb)*cos(pi*r(n_r))-sin(pi*r(n_r)))
! !                E_temp_radial(n_r)   =E_temp_radial(n_r)+cc2real(y_end(n_m,n_r), m)
! !          print*,"temp_Mloc(nm=",n_m,",nr=",n_r,")=",temp_Mloc(n_m,n_r)
! !          end do
!       end do


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
