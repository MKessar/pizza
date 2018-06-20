module galerkin

   use precision_mod
   use constants, only: one, two
   use namelists, only: radratio
   use band_matrix, only: type_bandmat_real

   implicit none

   private

   interface galerkin2cheb
      module procedure galerkin2cheb_real
      module procedure galerkin2cheb_complex
   end interface galerkin2cheb

   public :: galerkin2cheb, get_galerkin_stencil, destroy_galerkin_stencil

contains

   subroutine get_galerkin_stencil(sten, n_r_max, i_bound_type)

      !-- Input variables
      integer, intent(in) :: i_bound_type
      integer, intent(in) :: n_r_max

      !-- Output variables
      type(type_bandmat_real), intent(out) :: sten

      !-- Local variables
      real(cp) :: eta, dm
      integer :: kl, ku, i, m

      eta = radratio

      select case (i_bound_type)

         case(1) ! u(ri)=u(ro)=0 ! Dirichlet on both boundaries

            kl = 2
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               sten%dat(1,i)=-one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = one
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

         case(2) ! u'(ri)=u'(ro)=0 ! Neumann on both boundaries

            kl = 2
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               m = i-1
               sten%dat(1,i)=real((m+2)*(m+2),cp)/real(2*m*m+4*m+4,cp)
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i)=-real(m*m,cp)/real(2*m*m+4*m+4,cp)
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

         case(3) ! u(ro)=u'(ri)=0 ! Dirichlet outer/ Neumann inner

            kl = 2
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               m = i-1
               sten%dat(1,i)=one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max ) then
                  sten%dat(2,i) = -real(4*m+4,cp)/real(2*m*m+6*m+5,cp)
                  if ( i == n_r_max-1 ) sten%dat(2,i)=two*sten%dat(2,i)
               end if
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -real(2*m*m+2*m+1,cp)/real(2*m*m+6*m+5,cp)
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

         case(4) ! u'(ro)=u(ri)=0 ! Neumann outer/ Dirichlet inner

            kl = 2
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               m = i-1
               sten%dat(1,i)=one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max ) then
                  sten%dat(2,i) = real(4*m+4,cp)/real(2*m*m+6*m+5,cp)
                  if ( i == n_r_max-1 ) sten%dat(2,i)=two*sten%dat(2,i)
               end if
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -real(2*m*m+2*m+1,cp)/real(2*m*m+6*m+5,cp)
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
            end do

         case(5) ! u(ri)=u'(ri)=u(ro)=u'(ro)=0
            
            kl = 4
            ku = 0
            call sten%initialize(kl, ku, n_r_max)
            do i=1,n_r_max
               m = i-1
               sten%dat(1,i) = one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -two*real(m+2,cp)/real(m+3,cp)
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
               if ( i < n_r_max-3 ) then
                  sten%dat(5,i) = real(m+1,cp)/real(m+3,cp)
                  if ( i == n_r_max-4 ) sten%dat(5,i)=two*sten%dat(5,i)
               end if
            end do

         case(6) ! u(ri)=u''(ri)+u'(ri)=u(ro)=u''(ro)-u(ro)=0
            
            kl = 4
            ku = 0
            call sten%initialize(kl, ku, n_r_max)

            do i=1,n_r_max
               dm = real(i-1,cp)
               sten%dat(1,i)=one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max ) then
                  sten%dat(2,i) = -24.0_cp*(eta-1.0_cp)*(eta+1.0_cp)*(dm+1.0_cp)/ &
                  &               (12.0_cp*eta**2*dm**2+60.0_cp*eta**2*dm+        &
                  &                75.0_cp*eta**2+16.0_cp*eta*dm**4+160.0_cp      &
                  &               *eta*dm**3+584.0_cp*eta*dm**2+920.0_cp*eta*dm   &
                  &               +534.0_cp*eta+12.0_cp*dm**2+60.0_cp*dm+75.0_cp)
                  if ( i == n_r_max-1 ) sten%dat(2,i)=two*sten%dat(2,i)
               end if
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -2.0_cp*(dm+2.0_cp)*(12.0_cp*eta**2*dm**2+      &
                  &               48.0_cp*eta**2*dm+63.0_cp*eta**2+16.0_cp*eta*   &
                  &               dm**4+128.0_cp*eta*dm**3+424.0_cp*eta*dm**2+    &
                  &               672.0_cp*eta*dm+414.0_cp*eta+12.0_cp*dm**2+     &
                  &               48.0_cp*dm+63.0_cp)/((dm+3.0_cp)*(12.0_cp*      &
                  &               eta**2*dm**2+60.0_cp*eta**2*dm+75.0_cp*eta**2+  &
                  &               16.0_cp*eta*dm**4+160.0_cp*eta*dm**3+584.0_cp   &
                  &               *eta*dm**2+920.0_cp*eta*dm+534.0_cp*eta+12.0_cp &
                  &               *dm**2+60.0_cp*dm+75.0_cp))
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
               if ( i < n_r_max-2 ) then
                  sten%dat(4,i) = 24.0_cp*(eta-1.0_cp)*(eta+1.0_cp)*(dm+1.0_cp)/  &
                  &               (12.0_cp*eta**2*dm**2+60.0_cp*eta**2*dm+75.0_cp*&
                  &               eta**2+16.0_cp*eta*dm**4+160.0_cp*eta*dm**3+    &
                  &               584.0_cp*eta*dm**2+920.0_cp*eta*dm+534.0_cp*eta+&
                  &               12.0_cp*dm**2+60.0_cp*dm+75.0_cp)
                  if ( i == n_r_max-3 ) sten%dat(4,i)=two*sten%dat(4,i)
               end if
               if ( i < n_r_max-3 ) then
                  sten%dat(5,i) = (dm+1.0_cp)*(12.0_cp*eta**2*dm**2+36.0_cp*      &
                  &               eta**2*dm+27.0_cp*eta**2+16.0_cp*eta*dm**4+     &
                  &               96.0_cp*eta*dm**3+200.0_cp*eta*dm**2+168.0_cp   &
                  &               *eta*dm+54.0_cp*eta+12.0_cp*dm**2+36.0_cp*dm    &
                  &               +27.0_cp)/((dm+3.0_cp)*(12.0_cp*eta**2*dm**2    &
                  &               +60.0_cp*eta**2*dm+75.0_cp*eta**2+16.0_cp*eta   &
                  &               *dm**4+160.0_cp*eta*dm**3+584.0_cp*eta*dm**2    &
                  &               +920.0_cp*eta*dm+534.0_cp*eta+12.0_cp*dm**2+    &
                  &               60.0_cp*dm+75.0_cp))
                  if ( i == n_r_max-4 ) sten%dat(5,i)=two*sten%dat(5,i)
               end if
            end do

         case(7) ! u(ri)==u(ro)=u'(ri)=u''''(ro)=0

            kl = 4
            ku = 0
            call sten%initialize(kl, ku, n_r_max)

            do i=1,n_r_max
               dm = real(i-1,cp)
               sten%dat(1,i)=one
               if ( i == 1 .or. i==n_r_max ) sten%dat(1,i)=two*sten%dat(1,i)
               if ( i < n_r_max ) then
                  sten%dat(2,i) = 6.0_cp*(2.0_cp*dm**4+16.0_cp*dm**3+             &
                  &               50.0_cp*dm**2+72.0_cp*dm+35.0_cp)/((dm+4.0_cp)* &
                  &               (2.0_cp*dm**4+20.0_cp*dm**3+80.0_cp*dm**2+      &
                  &               150.0_cp*dm+105.0_cp))
                  if ( i == n_r_max-1 ) sten%dat(2,i)=two*sten%dat(2,i)
               end if
               if ( i < n_r_max-1 ) then
                  sten%dat(3,i) = -2.0_cp*(dm+2.0_cp)*(2.0_cp*dm**4+16.0_cp*dm**3+&
                  &               64.0_cp*dm**2+128.0_cp*dm+105.0_cp)/((dm+4.0_cp)&
                  &               *(2.0_cp*dm**4+20.0_cp*dm**3+80.0_cp*dm**2+     &
                  &               150.0_cp*dm+105.0_cp))
                  if ( i == n_r_max-2 ) sten%dat(3,i)=two*sten%dat(3,i)
               end if
               if ( i < n_r_max-2 ) then
                  sten%dat(4,i) = -6.0_cp*(2.0_cp*dm**4+16.0_cp*dm**3+50.0_cp*    &
                  &               dm**2+72.0_cp*dm+35.0_cp)/((dm+4.0_cp)*(2.0_cp* &
                  &               dm**4+20.0_cp*dm**3+80.0_cp*dm**2+150.0_cp*dm+  &
                  &               105.0_cp))
                  if ( i == n_r_max-3 ) sten%dat(4,i)=two*sten%dat(4,i)
               end if
               if ( i < n_r_max-3 ) then
                  sten%dat(5,i) = dm*(2.0_cp*dm**4+12.0_cp*dm**3+32.0_cp*dm**2+   &
                  &               42.0_cp*dm+17.0_cp)/((dm+4.0_cp)*(2.0_cp*dm**4+ &
                  &               20.0_cp*dm**3+80.0_cp*dm**2+150.0_cp*dm+105.0_cp))
                  if ( i == n_r_max-4 ) sten%dat(5,i)=two*sten%dat(5,i)
               end if
            end do

      end select

   end subroutine get_galerkin_stencil
!-------------------------------------------------------------------------------
   subroutine galerkin2cheb_real(sten, arr)

      !-- Input variable
      type(type_bandmat_real), intent(in) :: sten

      !-- Output variable
      real(cp), intent(inout) :: arr(sten%nlines)

      call sten%mat_vec_mul(arr)

   end subroutine galerkin2cheb_real
!-------------------------------------------------------------------------------
   subroutine galerkin2cheb_complex(sten, arr)

      !-- Input variable
      type(type_bandmat_real), intent(in) :: sten

      !-- Output variable
      complex(cp), intent(inout) :: arr(sten%nlines)

      call sten%mat_vec_mul(arr)

   end subroutine galerkin2cheb_complex
!-------------------------------------------------------------------------------
   subroutine destroy_galerkin_stencil(sten)
      !
      ! Memory deallocation of the stencils
      !
      type(type_bandmat_real), intent(in) :: sten

      call sten%finalize()

   end subroutine destroy_galerkin_stencil
!-------------------------------------------------------------------------------
end module galerkin
