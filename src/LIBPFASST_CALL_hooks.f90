!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pf_mod_dtype
!   use pf_mod_zndarray
  use pizza_zdnarray
  use pf_mod_utils
  implicit none
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
    use namelists, only: imex_stat, n_modes_m,n_points_phi,n_points_r
    use communications, only: reduce_radial_on_rank
    use useful, only: round_off, cc2real, cc22real, getMSD2, abortRun
    use integration, only: rInt_R
    use radial_functions, only: r, rscheme
    use truncation, only: idx2m
    use parallel_mod
    use blocking, only: nMstart, nMstop!,nRstart, nRstop

    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    complex(pfdp), pointer :: y_end(:,:)!,y_ex(:,:)
!     real(pfdp) :: maxerr
!     type(pf_zndarray_t), target :: y_exact      !<  the initial condition
    real(pfdp) :: E_temp_radial(n_points_r(level_index))
    real(pfdp) :: E_temp
    integer :: n_r,m,n_m
    character(len=144) :: file_name
    integer :: n_temp_file_2D
    integer ::  ierror
    
    y_end => get_array2d(pf%levels(level_index)%qend)
!     print*,"0 Key_Pizza =",Key_Pizza ,"pf%state%step+1=",pf%state%step+1, "pf%state%iter=",pf%state%iter

      do n_r=1,n_points_r(level_index)
         E_temp_radial(n_r) = 0.0_pfdp
!          do n_m=1,n_modes_m(level_index)
!          do n_m=nMstart, nMstop
         do n_m=1, (nMstop-nMstart+1)
            m = idx2m(n_m+nMstart-1)

               E_temp_radial(n_r)   = E_temp_radial(n_r)+cc2real(y_end(n_m,n_r), m)
!          print*,"y_end(nm=",n_m,",nr=",n_r,")=",y_end(n_m,n_r)
         end do
      end do
!     print*,"a Key_Pizza =",Key_Pizza ,"pf%state%step+1=",pf%state%step+1, "pf%state%iter=",pf%state%iter
!       do n_r=1,n_points_r(level_index)
!           print*,"Key_Pizza =",Key_Pizza ,"E_temp_radial(nr=",n_r")=", E_temp_radial(n_r)
!       enddo
    call reduce_radial_on_rank(E_temp_radial, 0)

    call mpi_barrier(Comm_Pizza, ierror)

    print*,"b Key_Pizza =",Key_Pizza ,"pf%state%step+1=",pf%state%step+1, "pf%state%iter=",pf%state%iter," res: ",pf%levels(level_index)%residual
       if ( Key_Pizza == 0 ) then
! 
            E_temp=rInt_R(E_temp_radial, r, rscheme)

    
         write(file_name, '(A,I0,A,A)') 'e_temp.test'
         open(unit=n_temp_file_2D, file=file_name, position='append')
! ! !          write(n_temp_file_2D, '(1P, es20.12, 3es16.8)') pf%state%step+1, E_temp 
! ! !          write(n_temp_file_2D, '(1P, i6.3, 3es16.8)') pf%state%step+1, E_temp 
         write(n_temp_file_2D, '(1P, 3es16.8, 3es16.8)') pf%state%t0+pf%state%dt, E_temp 
         close(n_temp_file_2D)
! 
!     print '("time_series: step: ",i6.3," iter: ",i4.3," level: ",i2.2," res: ",es14.7)', &
!      &    pf%state%step+1, pf%state%iter,level_index, pf%levels(level_index)%residual
!     
    endif
     
!     print '("time_series: step: ",i6.3," iter: ",i4.3," level: ",i2.2," error: ",es14.7," res: ",es14.7)', &
!      &    pf%state%step+1, pf%state%iter,level_index, maxerr,pf%levels(level_index)%residual
!     call flush(6)
!     call zndarray_destroy(y_exact)    
  end subroutine time_series
end module hooks
