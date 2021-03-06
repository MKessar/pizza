module update_temp_coll

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: one, zero, four, ci
   use namelists, only: kbott, ktopt, tadvz_fac, TdiffFac, BuoFac, l_buo_imp,nnodes
   use radial_functions, only: rscheme, or1, or2, dtcond, tcond, beta, &
       &                       rgrav
   use horizontal, only: hdif_T, bott_Mloc, topt_Mloc
   use blocking, only: nMstart, nMstop
   use truncation, only: n_r_max, idx2m, m2idx
   use radial_der, only: get_ddr, get_dr
   use fields, only: work_Mloc
   use algebra, only: prepare_full_mat, solve_full_mat
   use useful, only: abortRun
   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray

   implicit none
   
   private

   logical,  allocatable :: lTmat(:,:)
   real(cp), allocatable :: tMat(:, :, :,:)
   integer,  allocatable :: tPivot(:, :,:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: tMat_fac(:,:,:)
#endif
   complex(cp), allocatable :: rhs(:)

   public :: update_temp_co, initialize_temp_coll, finalize_temp_coll, &
   &         get_temp_rhs_imp_coll, finish_exp_temp_coll

contains

   subroutine initialize_temp_coll

      allocate( lTmat(nMstart:nMstop,(nnodes-1)) )
      lTmat(:,:)=.false.
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*SIZEOF_LOGICAL

      allocate( tMat(n_r_max, n_r_max, nMstart:nMstop,(nnodes-1) ))
      allocate( tPivot(n_r_max, nMstart:nMstop,(nnodes-1)) )
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*n_r_max* &
      &                 SIZEOF_DEF_REAL+n_r_max*(nMstop-nMstart+1)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate( tMat_fac(n_r_max, nMstart:nMstop,(nnodes-1) ))
      bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*  &
      &                 SIZEOF_DEF_REAL
#endif
      allocate( rhs(n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_temp_coll
!------------------------------------------------------------------------------
   subroutine finalize_temp_coll

      deallocate( rhs )
#ifdef WITH_PRECOND_S
      deallocate( tMat_fac )
#endif
      deallocate( lTmat, tMat, tPivot )

   end subroutine finalize_temp_coll
!------------------------------------------------------------------------------
   subroutine update_temp_co(temp_Mloc, dtemp_Mloc, buo_Mloc, dTdt, &
              &               lMat, l_log_next,tscheme,dtq,work_Mloc_pfasst,int_mat)

      !-- Input variables
      class(type_tscheme), intent(in),optional :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_log_next
      real(cp),intent(in),optional    :: dtq
      integer,intent(in),optional     :: int_mat
      
      !-- Output variables
      complex(cp),       intent(out) :: temp_Mloc(nMstart:nMstop, n_r_max)
      complex(cp),       intent(out) :: dtemp_Mloc(nMstart:nMstop, n_r_max)
      type(type_tarray), intent(inout) :: dTdt
      complex(cp),       intent(inout) :: buo_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),optional :: work_Mloc_pfasst(nMstart:nMstop,n_r_max)
      
      !-- Local variables
      integer :: n_r, n_m, n_r_out, m
      integer :: i_mat

      if (present(int_mat)) then
         i_mat=int_mat
      else
         i_mat=1
      endif
      if ( lMat ) then
      lTMat(:,i_mat)=.false.
      endif
      
     if (present(dtq)) then
      
         work_Mloc=work_Mloc_pfasst
      
     else
         !-- Calculation of the implicit part
         call get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc,          &
              &                     dTdt%old(:,:,tscheme%istage),   &
              &                     dTdt%impl(:,:,tscheme%istage),  &
              &                     tscheme%l_imp_calc_rhs(tscheme%istage))

         !-- Now assemble the right hand side and store it in work_Mloc
         call tscheme%set_imex_rhs(work_Mloc, dTdt, nMstart, nMstop, n_r_max)
     endif
      do n_m=nMstart, nMstop

         m = idx2m(n_m)
      if (present(dtq)) then
         if ( .not. lTmat(n_m,i_mat) ) then
!          print*,"get matrice i_mat=",i_mat,"int_mat=",int_mat
#ifdef WITH_PRECOND_S
            call get_tempMat( m, tMat(:,:,n_m,i_mat), tPivot(:,n_m,i_mat), &
                 &           tMat_fac(:,n_m,i_mat),dtq=dtq)
#else
            call get_tempMat(m, tMat(:,:,n_m,i_mat), tPivot(:,n_m,i_mat),dtq=dtq)
#endif
            lTmat(n_m,i_mat)=.true.
         end if    
      else
        if ( .not. lTmat(n_m,i_mat) ) then
#ifdef WITH_PRECOND_S
           call get_tempMat( m, tMat(:,:,n_m,i_mat), tPivot(:,n_m,i_mat), &
                &           tMat_fac(:,n_m,i_mat),tscheme)
#else
           call get_tempMat( m, tMat(:,:,n_m,i_mat), tPivot(:,n_m,i_mat),tscheme)
#endif
           lTmat(n_m,i_mat)=.true.
        end if
      endif
         !-- Inhomogeneous B.Cs (if not zero)
         rhs(1)      =topt_Mloc(n_m)
         rhs(n_r_max)=bott_Mloc(n_m)
         do n_r=2,n_r_max-1
            rhs(n_r)=work_Mloc(n_m,n_r)
         end do

#ifdef WITH_PRECOND_S
         do n_r=1,n_r_max
            rhs(n_r) = tMat_fac(n_r,n_m,i_mat)*rhs(n_r)
         end do
#endif

         call solve_full_mat(tMat(:,:,n_m,i_mat), n_r_max, n_r_max, tPivot(:, n_m,i_mat), &
              &              rhs(:))

         do n_r_out=1,rscheme%n_max
            temp_Mloc(n_m, n_r_out)=rhs(n_r_out)
         end do

      end do

      !-- set cheb modes > rscheme%n_max to zero (dealiazing)
      do n_r_out=rscheme%n_max+1,n_r_max
         do n_m=nMstart,nMstop
            temp_Mloc(n_m,n_r_out)=zero
         end do
      end do

      !-- Bring temperature back to physical space
      call rscheme%costf1(temp_Mloc, nMstart, nMstop, n_r_max)

      !-- Assemble buoyancy in case this is treated implicitly
     if ( l_buo_imp ) then
        if (present(tscheme)) then 
     
           call tscheme%assemble_implicit_buo(buo_Mloc, temp_Mloc, dTdt,      &
                &                             BuoFac, rgrav, nMstart, nMstop, &
                &                             n_r_max, .true.)
        end if
     endif
      !-- Roll the arrays before filling again the first block
     if (present(tscheme)) then 
         call tscheme%rotate_imex(dTdt, nMstart, nMstop, n_r_max)
     endif
      !-- In case log is needed on the next iteration, recalculate dT/dr
      if ( l_log_next ) then
         call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, n_r_max, rscheme)
      end if

   end subroutine update_temp_co
!------------------------------------------------------------------------------
   subroutine finish_exp_temp_coll(temp_Mloc, us_Mloc, dVsT_Mloc, buo_Mloc, &
              &                    dtemp_exp_last)

      !-- Input variables
      complex(cp), intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dVsT_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      complex(cp), intent(inout) :: buo_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(inout) :: dtemp_exp_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m

      if ( .not. l_buo_imp ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               if ( m /= 0 ) then
                  buo_Mloc(n_m,n_r)=-rgrav(n_r)*or1(n_r) &
                  &                  *BuoFac*ci*real(m,cp)*temp_Mloc(n_m,n_r)
               end if
            end do
         end do
      end if

      !-- Finish calculation of advection
      call get_dr( dVsT_Mloc, work_Mloc, nMstart, nMstop, n_r_max, &
           &       rscheme, nocopy=.true. )

      !-- Finish calculation of the explicit part for current time step
      do n_r=1,n_r_max
         do n_m=nMstart, nMstop
            dtemp_exp_last(n_m,n_r)=dtemp_exp_last(n_m,n_r)         & 
            &                       -or1(n_r)*work_Mloc(n_m,n_r)    &
            &                       -us_Mloc(n_m,n_r)*(dtcond(n_r)- &
            &                       tadvz_fac*beta(n_r)*tcond(n_r))
         end do
      end do

   end subroutine finish_exp_temp_coll
!------------------------------------------------------------------------------
   subroutine get_temp_rhs_imp_coll(temp_Mloc, dtemp_Mloc, temp_last, &
              &                     dtemp_imp_Mloc_last, l_calc_lin_rhs)

      !-- Input variables
      complex(cp), intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      logical,     intent(in) :: l_calc_lin_rhs

      !-- Output variable
      complex(cp), intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: temp_last(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) :: dtemp_imp_Mloc_last(nMstart:nMstop,n_r_max)

      !-- Local variables
      integer :: n_r, n_m, m
      real(cp) :: dm2

      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            temp_last(n_m,n_r)=temp_Mloc(n_m,n_r)
         end do
      end do

      if ( l_calc_lin_rhs ) then

         call get_ddr(temp_Mloc, dtemp_Mloc, work_Mloc, nMstart, nMstop, &
              &       n_r_max, rscheme)
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m = idx2m(n_m)
               dm2 = real(m,cp)*real(m,cp)
               dtemp_imp_Mloc_last(n_m,n_r)=TdiffFac*hdif_T(n_m)* (        &
               &                                       work_Mloc(n_m,n_r)  &
               &                         +or1(n_r)*   dtemp_Mloc(n_m,n_r)  &
               &                     -dm2*or2(n_r)*    temp_Mloc(n_m,n_r) )
            end do
         end do

      end if

   end subroutine get_temp_rhs_imp_coll
!------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_tempMat( m, tMat, tPivot, tMat_fac,tscheme,dtq)
#else
   subroutine get_tempMat( m, tMat, tPivot,tscheme,dtq)
#endif

      !-- Input variables
      class(type_tscheme), intent(in),optional :: tscheme        ! time step
      integer,             intent(in) :: m
      real(cp), intent(in),optional :: dtq
!       real(cp),  optional :: testdtq

      !-- Output variables
      real(cp), intent(out) :: tMat(n_r_max,n_r_max)
      integer,  intent(out) :: tPivot(n_r_max)
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: tMat_fac(n_r_max)
#endif

      !-- Local variables
      integer :: nR_out, nR, info, n_m
      real(cp) :: dm2

      dm2 = real(m,cp)*real(m,cp)
      n_m = m2idx(m)

      !----- Boundary coditions:
      do nR_out=1,rscheme%n_max
         if ( ktopt == 1 ) then
            tMat(1,nR_out)=rscheme%rnorm*rscheme%rMat(1,nR_out)
         else
            tMat(1,nR_out)=rscheme%rnorm*rscheme%drMat(1,nR_out)
         end if
         if ( kbott == 1 ) then
            tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%rMat(n_r_max,nR_out)
         else
            tMat(n_r_max,nR_out)=rscheme%rnorm*rscheme%drMat(n_r_max,nR_out)
         end if
      end do

      if ( rscheme%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme%n_max+1,n_r_max
            tMat(1,nR_out)      =0.0_cp
            tMat(n_r_max,nR_out)=0.0_cp
         end do
      end if
!             if (present(dtq)) then
!             print*,"present(dtq)"
!             else
!             print*,"not present(dtq)"
!             endif
      !----- Other points:
      if (present(dtq)) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               tMat(nR,nR_out)= rscheme%rnorm * (                          &
               &                                 rscheme%rMat(nR,nR_out) - &
               &                              dtq*TdiffFac*hdif_T(n_m)*(   &
               &                               rscheme%d2rMat(nR,nR_out) + &
               &          or1(nR)*              rscheme%drMat(nR,nR_out) - &
               &      dm2*or2(nR)*               rscheme%rMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               tMat(nR,nR_out)= rscheme%rnorm * (                          &
               &                                 rscheme%rMat(nR,nR_out) - &
               &              tscheme%wimp_lin(1)*TdiffFac*hdif_T(n_m)*(   &
               &                               rscheme%d2rMat(nR,nR_out) + &
               &          or1(nR)*              rscheme%drMat(nR,nR_out) - &
               &      dm2*or2(nR)*               rscheme%rMat(nR,nR_out) ) )
            end do
         end do
      endif
      
      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         tMat(nR,1)      =rscheme%boundary_fac*tMat(nR,1)
         tMat(nR,n_r_max)=rscheme%boundary_fac*tMat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         tMat_fac(nR)=one/maxval(abs(tMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nR=1,n_r_max
         tMat(nR,:) = tMat(nR,:)*tMat_fac(nR)
      end do
#endif

      !----- LU decomposition:
      call prepare_full_mat(tMat,n_r_max,n_r_max,tPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix tMat!')
      end if

   end subroutine get_tempMat
!------------------------------------------------------------------------------
end module update_temp_coll
