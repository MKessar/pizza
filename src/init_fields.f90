module init_fields
   !
   ! This module sets the initial fields either by reading in a checkpoint file
   ! or by setting a given perturbation in temperature and/or velocity specified
   ! in the input namelist
   !

   use constants, only: zero, one, two, three, ci, pi, half
   use blocking, only: nRstart, nRstop
   use communications, only: transp_r2m, r2m_fields
   use radial_functions, only: r, rscheme, or1, or2, beta, dbeta
   use namelists, only: l_start_file, dtMax, init_t, amp_t, init_u, amp_u, &
       &                radratio, r_cmb, r_icb, l_cheb_coll, l_non_rot,    &
       &                l_reset_t, l_chem, l_heat, amp_xi, init_xi
   use outputs, only: n_log_file
   use parallel_mod!, only: rank
   use blocking, only: nMstart, nMstop, nM_per_rank
   use truncation, only: m_max, n_r_max, minc, m2idx, idx2m, n_phi_max
   use useful, only: logWrite, abortRun, gausslike_compact_center, &
       &             gausslike_compact_middle, gausslike_compact_edge
   use radial_der, only: get_dr
   use fourier, only: fft
   use checkpoints, only: read_checkpoint
   use time_schemes, only: type_tscheme
   use fields
   use fieldsLast
   use precision_mod
   use output_frames, only: write_snapshot_mloc

   implicit none

   private

   public :: get_start_fields

contains

   subroutine get_start_fields(time, tscheme)

      !-- Output variables
      real(cp),            intent(out) :: time
      class(type_tscheme), intent(inout) :: tscheme

      !-- Local variables
      integer :: m, n_r, n_m
      logical :: lMat
      real(cp) :: h2
      character(len=76) :: message

   
      if ( l_start_file ) then
         call read_checkpoint(us_Mloc, up_Mloc, temp_Mloc, xi_Mloc, dpsidt, &
              &               dTdt, dxidt, time, tscheme)

         if ( l_reset_t ) time = 0.0_cp

         !-- If integration method is used, since u_s is stored, one needs to
         !-- reconstruct psi(m/=0)
         if ( .not. l_cheb_coll ) then
            do n_r=2,n_r_max
               if ( l_non_rot ) then
                  h2 = one
               else
                  h2 = r_cmb*r_cmb-r(n_r)*r(n_r)
               end if
               do n_m=nMstart,nMstop
                  m = idx2m(n_m)
                  if ( m > 0 ) then
                     psi_Mloc(n_m, n_r) = -ci*r(n_r)/real(m,cp)/h2 * us_Mloc(n_m, n_r)
                  else
                     psi_Mloc(n_m, n_r) = 0.0_cp
                  end if
               end do
            end do
            !-- Boundary point (since h2 is singular there)
            do n_m=nMstart, nMstop
               psi_Mloc(n_m,1)=zero
            end do
         end if

      else

         if ( l_heat ) temp_Mloc(:,:)=zero
         if ( l_chem ) xi_Mloc(:,:)  =zero
         us_Mloc(:,:)  =zero
         up_Mloc(:,:)  =zero
         psi_Mloc(:,:) =zero
         call dpsidt%set_initial_values()
         if ( l_heat ) call dTdt%set_initial_values()
         if ( l_chem ) call dxidt%set_initial_values()

         time=0.0_cp
         tscheme%dt(:)=dtMax

         if (Key_Pizza == 0) write(message,'(''! Using dtMax time step:'',ES16.6)') dtMax
         call logWrite(message, n_log_file)

      end if


      !-- Initialize the weights of the time scheme
      call tscheme%set_weights(lMat)

      if ( init_t /= 0 .and. l_heat ) then 
         call initProp(temp_Mloc, init_t, amp_t)
    
         call write_snapshot_mloc("frame_temp_0.test", 0.0_cp, temp_Mloc)
         
      endif
      
      if ( init_xi /= 0 .and. l_chem ) call initProp(xi_Mloc, init_xi, amp_xi)

      if ( init_u /= 0 ) then
                    call initU(us_Mloc, up_Mloc)
      call write_snapshot_mloc("frame_us_0.test", 0.0_cp, us_Mloc)
      call write_snapshot_mloc("frame_up_0.test", 0.0_cp, up_Mloc)
      
      end if 



      !-- Reconstruct missing fields, dtemp_Mloc, om_Mloc
      if ( l_heat ) call get_dr(temp_Mloc, dtemp_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      if ( l_chem ) call get_dr(xi_Mloc, dxi_Mloc, nMstart, nMstop, &
                         &      n_r_max, rscheme)
      call get_dr(up_Mloc, work_Mloc, nMstart, nMstop, n_r_max, rscheme)
      do n_r=1,n_r_max
         do n_m=nMstart,nMstop
            m = idx2m(n_m)
            om_Mloc(n_m,n_r)=work_Mloc(n_m,n_r)+or1(n_r)*up_Mloc(n_m,n_r)- &
            &                     or1(n_r)*ci*real(m,cp)*us_Mloc(n_m,n_r)
         end do
      end do
      call write_snapshot_mloc("frame_om_0.test", 0.0_cp, om_Mloc)


      !-- When not using Collocation also store temp_hat and psi_hat
      !-- This saves 2 DCTs per iteration
      if ( .not. l_cheb_coll ) then
         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               if ( l_heat ) temp_hat_Mloc(n_m,n_r)=temp_Mloc(n_m,n_r)
               if ( l_chem ) xi_hat_Mloc(n_m,n_r)=xi_Mloc(n_m,n_r)
               psi_hat_Mloc(n_m,n_r) =psi_Mloc(n_m,n_r)
            end do
         end do
         if ( l_heat ) call rscheme%costf1(temp_hat_Mloc, nMstart, nMstop, &
                            &              n_r_max)
         if ( l_chem ) call rscheme%costf1(xi_hat_Mloc, nMstart, nMstop, &
                            &              n_r_max)
         call rscheme%costf1(psi_hat_Mloc, nMstart, nMstop, n_r_max)
      end if

   end subroutine get_start_fields
!----------------------------------------------------------------------------------
   subroutine initProp(prop_Mloc, init_prop, amp_prop)

      !-- Output variables
      complex(cp), intent(inout) :: prop_Mloc(nMstart:nMstop, n_r_max)
      real(cp),    intent(in) :: amp_prop
      integer,     intent(in) :: init_prop

      !-- Local variables
      integer :: m_pertu, n_r, idx, n_phi, n_m, ir
      real(cp) :: x, c_r, rdm, c1, c2, rc, L, sigma_r
      real(cp) :: t1(n_r_max), gasp(n_r_max)
      real(cp) :: phi, phi0
      real(cp) :: phi_func(n_phi_max)

      !-- Radial dependence of perturbation in t1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         t1(n_r)=sin(pi*(r(n_r)-r_icb))
      end do

      if ( init_prop > 0 ) then ! Initialize a peculiar m mode
         
         m_pertu = init_prop

         if ( mod(m_pertu,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for temperature initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(*,*) '! Degree of mode for temperature initialisation'
            write(*,*) '! > m_max  !',m_pertu
            call abortRun('Stop run in init')
         end if

         idx = m2idx(m_pertu)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=t1(n_r)*amp_prop
               prop_Mloc(idx,n_r)=prop_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

      else if ( init_prop == -1 ) then ! bubble = Gaussian in r and phi

         sigma_r = 0.1_cp/sqrt(two)
         rc = half*(r_cmb+r_icb)
         !-- Find the closest point to the middle radius
         ir=minloc(abs(r-rc),1)
         !-- Overwrite rc to be on the grid
         rc = r(ir)
         c1 = half*(r(1)-rc)
         c2 = half*(rc-r(n_r_max))
         L = half * sigma_r

         !-- Piecewise-definition of a compact-support Gaussian-like profile
         !-- From Eq. (4.7) from Gaspari et al. (1999)
         !-- This ensures that the BCs will always be fulfilled
         do n_r=1,n_r_max
            if ( n_r == ir ) then ! Middle point
               gasp(n_r)=one
            else ! Middle and Edge parts
               if ( r(n_r) < rc-c2 ) then
                  gasp(n_r) = gausslike_compact_edge(rc-r(n_r),c2,L)
               else if ( r(n_r) >= rc-c2 .and. r(n_r) < rc ) then
                  gasp(n_r) = gausslike_compact_middle(rc-r(n_r),c2,L)
               else if ( r(n_r) > rc .and. r(n_r) <= rc+c1 ) then
                  gasp(n_r) = gausslike_compact_middle(r(n_r)-rc,c1,L)
               else if ( r(n_r) > rc+c1 ) then
                  gasp(n_r) = gausslike_compact_edge(r(n_r)-rc,c1,L)
               end if
            end if
         end do

         !-- Normalisation of the two branches by the middle value
         do n_r=1,n_r_max
            if ( n_r < ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c1,L)
            else if ( n_r > ir ) then
               gasp(n_r) = gasp(n_r)/gausslike_compact_center(c2,L)
            end if
         end do

         !-- Now finally define the bubble
         phi0 = pi/minc
         do n_r=nRstart,nRstop
            c_r = amp_prop*gasp(n_r)
            do n_phi=1,n_phi_max
               phi = (n_phi-1)*two*pi/minc/(n_phi_max)
               phi_func(n_phi)=c_r*exp(-(phi-phi0)**2/0.2_cp**2)
            end do

            !-- temp_Rloc is used as a work r-distributed array here
            call fft(phi_func, temp_Rloc(:,n_r))
         end do

         !-- MPI transpose is needed here
         call transp_r2m(r2m_fields, temp_Rloc, work_Mloc)

         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               prop_Mloc(n_m,n_r) = prop_Mloc(n_m,n_r)+work_Mloc(n_m,n_r)
            end do
         end do

      else ! random noise

         do n_r=1,n_r_max
            t1(n_r)=sin(pi*(r(n_r)-r_icb))
         end do

         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m_pertu = idx2m(n_m)
               if ( m_pertu > 0 ) then
                  call random_number(rdm)
                  prop_Mloc(n_m, n_r) = amp_prop*rdm*m_pertu**(-1.5_cp)*t1(n_r)
               end if
            end do
         end do

      end if

   end subroutine initProp
!----------------------------------------------------------------------------------
   subroutine initU(us_Mloc, up_Mloc)

      !-- Output variables
      complex(cp), intent(inout) :: us_Mloc(nMstart:nMstop, n_r_max)
      complex(cp), intent(inout) :: up_Mloc(nMstart:nMstop, n_r_max)

      !-- Local variables
      integer :: m_pertu, n_r, idx, m, n_m,n_phi
      real(cp) :: c_r
      real(cp) :: u1(n_r_max)
      real(cp) :: phi, phi0
      real(cp) :: phi_func(n_phi_max)
      complex(cp) :: dpsi_Mloc(nMstart:nMstop,n_r_max)

      !-- Radial dependence of perturbation in t1:
      do n_r=1,n_r_max
         u1(n_r)=sin(pi*(r(n_r)-r_icb))
      end do


      if ( init_u > 0 ) then
         
         m_pertu = init_u

         if ( mod(m_pertu,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for velocity initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m_pertu
            call abortRun('Stop run in init')
         end if
         if ( m_pertu > m_max ) then
            write(*,*) '! Degree of mode for velocity initialisation'
            write(*,*) '! > m_max  !',m_pertu
            call abortRun('Stop run in init')
         end if

         idx = m2idx(m_pertu)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=u1(n_r)*amp_u
               us_Mloc(idx,n_r)=us_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

         !-- Get the corresponding vorticity
         do n_r=1,n_r_max
            do n_m=nMstart, nMstop
               m = idx2m(n_m)
               om_Mloc(n_m,n_r)=-ci*m*or1(n_r)*us_Mloc(n_m,n_r)
            end do
         end do
         
      elseif (init_u < 0 ) then
         
         
         phi0 = pi/minc
         do n_r=nRstart,nRstop
            do n_phi=1,n_phi_max
               phi = (n_phi-1)*two*pi/minc/(n_phi_max)
               phi_func(n_phi)=amp_u*(((r(n_r)-r_icb)*(r(n_r)-r_cmb))**2.0_cp)*cos((3.0_cp*phi)+pi/2.0_cp)
            end do

            !-- temp_Rloc is used as a work r-distributed array here
            call fft(phi_func, temp_Rloc(:,n_r))
         end do

         !-- MPI transpose is needed here
         call transp_r2m(r2m_fields, temp_Rloc, work_Mloc)
         call get_dr(work_Mloc, dpsi_Mloc, nMstart, nMstop, n_r_max, rscheme)

         do n_r=1,n_r_max
            do n_m=nMstart,nMstop
               m=idx2m(n_m)
               us_Mloc(n_m,n_r)=us_Mloc(n_m,n_r)+ci*m*or1(n_r)*work_Mloc(n_m,n_r)
               up_Mloc(n_m,n_r)=up_Mloc(n_m,n_r)-dpsi_Mloc(n_m,n_r)-beta(n_r)*work_Mloc(n_m,n_r)
            end do
         end do

         !--stop
      else ! initialize an axisymmetric vphi

         idx = m2idx(0)
         if ( idx >= nMstart .and. idx <= nMstop ) then
            do n_r=1,n_r_max
               c_r=amp_u*sin(pi*(r(n_r)-r_icb))
               up_Mloc(idx,n_r)=up_Mloc(idx,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do
         end if

      end if

   end subroutine initU
!----------------------------------------------------------------------------------
end module init_fields
