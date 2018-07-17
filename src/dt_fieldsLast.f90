module fieldsLast
   !
   ! This module contains time-derivatives array of the previous time-step
   ! They are needed in the time-stepping scheme.
   !
 
   use precision_mod
   use constants, only: zero
   use mem_alloc, only: bytes_allocated
   use time_array

   implicit none

   private

   complex(cp), allocatable, target, public :: dt_fields_container_Rloc(:,:,:)
   complex(cp), allocatable, target, public :: dt_fields_container_Mloc(:,:,:)
   complex(cp), pointer, public :: dVsT_Mloc(:,:)
   complex(cp), pointer, public :: dVsOm_Mloc(:,:)
   complex(cp), allocatable, public :: buo_Mloc(:,:)
   complex(cp), allocatable, public :: dpsidt_Rloc(:,:)
   complex(cp), allocatable, public :: dtempdt_Rloc(:,:)
   complex(cp), pointer, public  :: dVsT_Rloc(:,:)
   complex(cp), pointer, public  :: dVsOm_Rloc(:,:)
   type(type_tarray), public :: dpsidt
   type(type_tarray), public :: dTdt

   public :: initialize_fieldsLast, finalize_fieldsLast

contains

   subroutine initialize_fieldsLast(nMstart,nMstop,n_m_max,nRstart,nRstop,n_r_max,&
              &                     norder_imp, norder_exp, norder_imp_lin)

      !-- Input variables
      integer, intent(in) :: nMstart
      integer, intent(in) :: nMstop
      integer, intent(in) :: n_m_max
      integer, intent(in) :: nRstart
      integer, intent(in) :: nRstop
      integer, intent(in) :: n_r_max
      integer, intent(in) :: norder_imp
      integer, intent(in) :: norder_exp
      integer, intent(in) :: norder_imp_lin

      call dpsidt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
           &                 norder_imp_lin)
      call dTdt%initialize(nMstart, nMstop, n_r_max, norder_imp, norder_exp, &
           &               norder_imp_lin)

      allocate( dt_fields_container_Mloc(nMstart:nMstop,n_r_max,2) )
      dVsT_Mloc(nMstart:,1:)     => dt_fields_container_Mloc(nMstart:,1:,1)
      dVsOm_Mloc(nMstart:,1:)    => dt_fields_container_Mloc(nMstart:,1:,2)
      allocate( buo_Mloc(nMStart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + 3*(nMstop-nMStart+1)*n_r_max*&
      &                 SIZEOF_DEF_COMPLEX

      buo_Mloc(:,:)  =zero
      dVsT_Mloc(:,:) =zero
      dVsOm_Mloc(:,:)=zero

      allocate( dtempdt_Rloc(n_m_max, nRstart:nRstop) )
      allocate( dpsidt_Rloc(n_m_max, nRstart:nRstop) )
      allocate( dt_fields_container_Rloc(n_m_max,nRstart:nRstop,2) )
      dVsT_Rloc(1:,nRstart:)    => dt_fields_container_Rloc(1:,nRstart:,1)
      dVsOm_Rloc(1:,nRstart:)   => dt_fields_container_Rloc(1:,nRstart:,2)
      bytes_allocated = bytes_allocated + 4*(nRstop-nRStart+1)*n_m_max* &
      &                 SIZEOF_DEF_COMPLEX

      dpsidt_Rloc(:,:) =zero
      dtempdt_Rloc(:,:)=zero
      dVsT_Rloc(:,:)   =zero
      dVsOm_Rloc(:,:)  =zero

   end subroutine initialize_fieldsLast
!-------------------------------------------------------------------------------
   subroutine finalize_fieldsLast

      call dTdt%finalize()
      call dpsidt%finalize()
      deallocate( dtempdt_Rloc, dpsidt_Rloc )
      deallocate( dt_fields_container_Mloc )
      deallocate( dt_fields_container_Rloc )
      deallocate( buo_Mloc )

   end subroutine finalize_fieldsLast
!-------------------------------------------------------------------------------
end module fieldsLast
