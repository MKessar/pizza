!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
!   use pf_mod_zndarray
!   use pizza_zdnarray
  use pizza_zndsysarray
  use pf_mod_utils
  
  implicit none
  integer :: frame_counter
contains

  !>  Output the error and residual in the solution
!   subroutine echo_error(pf, level_index)
! !   use pf_my_sweeper, only: exact
!     type(pf_pfasst_t), intent(inout) :: pf
!     integer, intent(in) :: level_index
! 
!     complex(pfdp), pointer :: y_end(:,:),y_ex(:,:)
!     real(pfdp) :: maxerr
!     type(pf_zndarray_t), target :: y_exact      !<  the initial condition
!     y_end => get_array2d(pf%levels(level_index)%qend)
!     call zndarray_build(y_exact, [ pf%levels(level_index)%lev_shape])
!     y_ex => get_array2d(y_exact)    
!     
!     !>  compute the exact solution
! !     call exact(pf%state%t0+pf%state%dt, y_ex)
!     !>  compute error
!     maxerr = maxval(abs(y_end-y_ex))
!     
!     print '("error: step: ",i3.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
!          pf%state%step+1, pf%state%iter,level_index, maxerr,pf%levels(level_index)%residual
! !     call flush(6)
!     call zndarray_destroy(y_exact)    
!   end subroutine echo_error
  
  
 subroutine time_series(pf, level_index)
!     use pf_my_sweeper, only: exact
    use namelists, only: imex_stat, n_modes_m,n_points_phi,n_points_r,n_frame_step
    use communications, only: reduce_radial_on_rank
    use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
    use integration, only: rInt_R
    use radial_functions, only: r, rscheme
    use truncation, only: idx2m
    use parallel_mod
    use blocking, only: nMstart, nMstop!,nRstart, nRstop
    use output_frames, only: write_snapshot_mloc
    use outputs, only: get_time_series
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

!     complex(pfdp), pointer :: y_end(:,:)!,y_ex(:,:)
    complex(pfdp),      pointer :: theta(:,:),omega(:,:)
    complex(pfdp),      pointer :: psi(:,:),uphi0(:,:)
    complex(pfdp),      pointer :: u_phi(:,:),u_s(:,:)
    complex(pfdp), allocatable  :: Temp_Mloc(:,:)
    complex(pfdp), allocatable  :: u_phi_Mloc(:,:)
    complex(pfdp), allocatable  :: u_s_Mloc(:,:)
    complex(pfdp), allocatable  :: omega_Mloc(:,:)

!     real(pfdp) :: maxerr
!     type(pf_zndarray_t), target :: y_exact      !<  the initial condition
    real(pfdp) :: E_temp_radial(n_points_r(level_index))
    real(pfdp) :: E_temp
    real(pfdp) :: E_kin_radial(n_points_r(level_index))
    real(pfdp) :: E_kin
    real(pfdp) :: us2_r(n_points_r(level_index))
    real(pfdp) :: up2_r(n_points_r(level_index))
    real(pfdp) :: enstrophy_r(n_points_r(level_index))
    real(pfdp) :: flux_r(n_points_r(level_index))
    
    integer :: n_r,m,n_m
    character(len=144) :: file_name
    integer :: n_temp_file_2D
    integer ::  ierror
    character(len=144) :: frame_name
    integer :: step,mod_step
    allocate(Temp_Mloc(nMstart:nMstop,n_points_r(level_index)))
    allocate(u_phi_Mloc(nMstart:nMstop,n_points_r(level_index)))
    allocate(u_s_Mloc(nMstart:nMstop,n_points_r(level_index)))
    allocate(omega_Mloc(nMstart:nMstop,n_points_r(level_index)))

!     y_end => get_array2d(pf%levels(level_index)%qend)
!     print*,"time_series -1"
!     print*,"size omega",SIZE(omega)
!     print*,"size psi",SIZE(psi)
!     print*,"size theta",SIZE(theta)
!     print*,"size uphi0",SIZE(uphi0)
    omega  => get_array2d(pf%levels(level_index)%qend,1)
    u_phi  => get_array2d(pf%levels(level_index)%qend,2)
    u_s    => get_array2d(pf%levels(level_index)%qend,3)
!     psi    => get_array2d(pf%levels(level_index)%qend,2)
    theta  => get_array2d(pf%levels(level_index)%qend,4)
!     uphi0  => get_array2d(pf%levels(level_index)%qend,4)

!     print*,"time_series 0"
!     print*,"size omega",SIZE(omega)
!     print*,"size psi",SIZE(psi)
!     print*,"size theta",SIZE(theta)
!     print*,"size uphi0",SIZE(uphi0)
!     print*,"size pf%levels(level_index)%qend%arr_shape(1)",pf%levels(level_index)%qend%arr_shape
!     print*,"size pf%levels(level_index)%qend%arr_shape(2)",pf%levels(level_index)%qend%arr_shape
    
    Temp_Mloc  = theta
    u_phi_Mloc = u_phi
    u_s_Mloc   = u_s
    omega_Mloc = omega
!     
!       do n_r=1,n_points_r(level_index)
!          E_temp_radial(n_r) = 0.0_pfdp
! !          E_kin_radial(n_r) = 0.0_pfdp
! !          do n_m=1, (nMstop-nMstart+1)
!          do n_m=nMstart,nMstop
!             m = idx2m(n_m+nMstart-1)
! 
! !               E_temp_radial(n_r)   = E_temp_radial(n_r)+cc2real(theta(n_m,n_r), m)
!               E_temp_radial(n_r)   = E_temp_radial(n_r)+ cc2real(Temp_Mloc(n_m,n_r), m)
! !               E_kin_radial(n_r)    = E_kin_radial(n_r) + cc2real(u_phi_Mloc(n_m,n_r), m) &
! !                                    &                   + cc2real(u_s_Mloc(n_m,n_r), m)
!          end do
!       end do
!       
!     call reduce_radial_on_rank(E_temp_radial, 0)
! !     call reduce_radial_on_rank(E_kin_radial, 0)
! 
!     call mpi_barrier(Comm_Pizza, ierror)
! 
!        if ( Key_Pizza == 0 ) then
! ! 
!             E_temp = rInt_R(E_temp_radial, r, rscheme)
! !             E_kin  = rInt_R(E_kin_radial, r, rscheme)
! 
!     
!          write(file_name, '(A,I0,A,A)') 'e_temp.test'
!          open(unit=n_temp_file_2D, file=file_name, position='append')
!          write(n_temp_file_2D, '(1P, 3es16.8, 3es16.8)') pf%state%t0+pf%state%dt, E_temp 
!          close(n_temp_file_2D)
! 
! !          write(file_name, '(A,I0,A,A)') 'e_kin.test'
! !          open(unit=n_temp_file_2D, file=file_name, position='append')
! !          write(n_temp_file_2D, '(1P, 3es16.8, 3es16.8)') pf%state%t0+pf%state%dt, E_kin 
! !          close(n_temp_file_2D)
! 
! 
! ! 
! !     print '("time_series: step: ",i6.3," iter: ",i4.3," level: ",i2.2," res: ",es14.7)', &
! !      &    pf%state%step+1, pf%state%iter,level_index, pf%levels(level_index)%residual
! !     
!     endif

    
    if (pf%state%step.EQ.0) then
    frame_counter=1
    endif 

    step=pf%state%step+1
    mod_step=mod(step,n_frame_step)

    if ( mod_step .EQ. 0 )  then
! 
    !        write(frame_name, '(A,I0,A,I0,A,A)') 'frame_temp_',frame_counter,'TimeProc_',Color_Pizza,'.test'
!     Temp_Mloc=theta
!        write(frame_name, '(A,I0,A,I0,A,A)') 'frame_temp_',frame_counter,'.test'
         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_temp_',frame_counter,'TimeProc_',Color_Pizza,'.test'
       call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, Temp_Mloc)   
         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_us_',frame_counter,'TimeProc_',Color_Pizza,'.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, u_s_Mloc)
         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_up_',frame_counter,'TimeProc_',Color_Pizza,'.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, u_phi_Mloc)
         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_om_',frame_counter,'TimeProc_',Color_Pizza,'.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, omega_Mloc)       
       frame_counter= frame_counter +1
       
    endif
 
!     Temp_Mloc  = theta
!     u_phi_Mloc = u_phi
!     u_s_Mloc   = u_s
!     omega_Mloc = omega
    
    call get_time_series(pf%state%t0+pf%state%dt, u_s_Mloc, u_phi_Mloc, omega_Mloc, Temp_Mloc, &
              &    us2_r, up2_r, enstrophy_r, flux_r)
   
!     call write_outputs(pf%state%t0+pf%state%dt, n_time_step, l_log, l_rst, l_frame, &
!               &             l_vphi_bal_write, l_stop_time,  us_Mloc, up_Mloc,  &
!               &             om_Mloc, temp_Mloc, dtemp_Mloc, xi_Mloc, dxi_Mloc, &
!               &             dpsidt, dTdt, dxidt)

!     do n_m=nMstart,nMstop
! !       do n_r=1,n_points_r(level_index)
!             m = idx2m(n_m)
!              if ( m == 3 ) then
!                  print*,"pf%state%t0+pf%state%dt=",pf%state%t0+pf%state%dt,"u_phi_Mloc(nm=",n_m,",nr=",n_points_r(level_index),")=",u_phi_Mloc(n_m,n_points_r(level_index))
!              endif
! !          end do
!       end do
    
    deallocate(Temp_Mloc)
    deallocate(u_phi_Mloc)
    deallocate(u_s_Mloc)
    deallocate(omega_Mloc)
    
!     print '("time_series: step: ",i6.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
!      &    pf%state%step+1, pf%state%iter,level_index, maxerr,pf%levels(level_index)%residual
!     call flush(6)
  end subroutine time_series
end module hooks
