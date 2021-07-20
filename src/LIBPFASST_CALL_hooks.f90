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

  
  
 subroutine time_series(pf, level_index)

    use namelists, only: imex_stat, n_modes_m,n_points_phi,n_points_r,n_frame_step, l_chem
    use communications, only: reduce_radial_on_rank
    use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
    use integration, only: rInt_R
    use radial_functions, only: r, rscheme
    use truncation, only: idx2m
    use parallel_mod
    use blocking, only: nMstart, nMstop!,nRstart, nRstop
    use output_frames, only: write_snapshot_mloc
    use outputs, only: get_time_series
    use radial_der, only: get_dr
    
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    complex(pfdp),      pointer :: theta(:,:),omega(:,:),xi(:,:)
    complex(pfdp),      pointer :: psi(:,:),uphi0(:,:)
    complex(pfdp),      pointer :: u_phi(:,:),u_s(:,:)
    complex(pfdp), allocatable  :: Temp_Mloc(:,:)
    complex(pfdp), allocatable  :: u_phi_Mloc(:,:)
    complex(pfdp), allocatable  :: u_s_Mloc(:,:)
    complex(pfdp), allocatable  :: omega_Mloc(:,:)
    complex(pfdp), allocatable  :: xi_Mloc(:,:)
    complex(pfdp), allocatable  :: dxi_Mloc(:,:)

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
    allocate(xi_Mloc(nMstart:nMstop,n_points_r(level_index)))
    allocate(dxi_Mloc(nMstart:nMstop,n_points_r(level_index)))


    omega => get_array2d(pf%levels(level_index)%qend,1)
    u_phi => get_array2d(pf%levels(level_index)%qend,2)
    u_s   => get_array2d(pf%levels(level_index)%qend,3)
    theta => get_array2d(pf%levels(level_index)%qend,4)
    if ( l_chem ) then
       xi => get_array2d(pf%levels(level_index)%qend,5)
    endif
    
    Temp_Mloc  = theta
    u_phi_Mloc = u_phi
    u_s_Mloc   = u_s
    omega_Mloc = omega
    if ( l_chem ) then
       xi_Mloc = xi
    endif    
    
    if (pf%state%step.EQ.0) then
       frame_counter = 1
    endif 
    
    step     = pf%state%step + 1
    mod_step = mod(step,n_frame_step)
    
    if ( mod_step .EQ. 0 )  then

    
         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_temp_', frame_counter, 'TimeProc_', Color_Pizza, '.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, Temp_Mloc)   

         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_us_', frame_counter, 'TimeProc_', Color_Pizza, '.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, u_s_Mloc)

         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_up_', frame_counter, 'TimeProc_', Color_Pizza, '.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, u_phi_Mloc)

         write(frame_name, '(A,I0,A,I0,A,A)') 'frame_om_', frame_counter, 'TimeProc_', Color_Pizza, '.test'
         call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, omega_Mloc)   

         if ( l_chem ) then
            write(frame_name, '(A,I0,A,I0,A,A)') 'frame_xi_', frame_counter, 'TimeProc_', Color_Pizza, '.test'
            call write_snapshot_mloc(frame_name, pf%state%t0+pf%state%dt, xi_Mloc)
         end if

         frame_counter = frame_counter + 1

    endif
    
    call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, n_points_r(level_index), rscheme)

    call get_time_series(pf%state%t0+pf%state%dt, u_s_Mloc, u_phi_Mloc, omega_Mloc, Temp_Mloc, &
              &    us2_r, up2_r, enstrophy_r, flux_r,xi_Mloc=xi_Mloc,dxi_Mloc=dxi_Mloc)
   
    
    deallocate(Temp_Mloc)
    deallocate(u_phi_Mloc)
    deallocate(u_s_Mloc)
    deallocate(omega_Mloc)
    deallocate(xi_Mloc)
    deallocate(dxi_Mloc)
!     
  end subroutine time_series
end module hooks
