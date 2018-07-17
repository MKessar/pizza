module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use constants, only: zero
   use namelists, only: l_cheb_coll
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_m_max, n_r_max
   use blocking, only: nMstart, nMstop, nRstart, nRstop
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), allocatable, public :: psi_Mloc(:,:)
   complex(cp), allocatable, public :: dtemp_Mloc(:,:)
   complex(cp), allocatable, public :: dom_Mloc(:,:)
   complex(cp), allocatable, public :: work_Mloc(:,:)
   complex(cp), pointer, public :: temp_Mloc(:,:), om_Mloc(:,:)
   complex(cp), pointer, public :: us_Mloc(:,:), up_Mloc(:,:)
   complex(cp), allocatable, target, public :: fields_container_Mloc(:,:,:)

   complex(cp), pointer, public :: us_Rloc(:,:), up_Rloc(:,:)
   complex(cp), pointer, public :: om_Rloc(:,:), temp_Rloc(:,:)
   complex(cp), allocatable, target, public :: fields_container_Rloc(:,:,:)

   complex(cp), allocatable, public :: psi_hat_Mloc(:,:), temp_hat_Mloc(:,:)
 
   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      allocate( fields_container_Mloc(nMstart:nMstop, n_r_max, 4) )
      allocate( dtemp_Mloc(nMStart:nMstop,n_r_max) )
      allocate( psi_Mloc(nMStart:nMstop,n_r_max) )
      allocate( dom_Mloc(nMStart:nMstop,n_r_max) )
      allocate( work_Mloc(nMStart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated + &
      &                 8*(nMstop-nMstart+1)*n_r_max*SIZEOF_DEF_COMPLEX

      temp_Mloc(nMstart:,1:) => fields_container_MLoc(nMstart:,1:,1)
      us_Mloc(nMstart:,1:)   => fields_container_MLoc(nMstart:,1:,2)
      up_Mloc(nMstart:,1:)   => fields_container_MLoc(nMstart:,1:,3)
      om_Mloc(nMstart:,1:)   => fields_container_MLoc(nMstart:,1:,4)

      temp_Mloc(:,:) =zero
      dtemp_Mloc(:,:)=zero
      psi_Mloc(:,:)  =zero
      om_Mloc(:,:)   =zero
      dom_Mloc(:,:)  =zero
      us_Mloc(:,:)   =zero
      up_Mloc(:,:)   =zero
      work_Mloc(:,:) =zero

      if ( .not. l_cheb_coll ) then
         allocate( psi_hat_Mloc(nMstart:nMstop,n_r_max) )
         allocate( temp_hat_Mloc(nMstart:nMstop,n_r_max) )
         bytes_allocated = bytes_allocated+(nMstop-nMstart+1)*n_r_max*&
         &                 SIZEOF_DEF_COMPLEX
         psi_hat_Mloc(:,:) =zero
         temp_hat_Mloc(:,:)=zero
      else
         ! We have to allocate the arrays to make the code runs when
         ! debug flags are turned on
         allocate (psi_hat_Mloc(0,0), temp_hat_Mloc(0,0))
      end if

      allocate( fields_container_Rloc(n_m_max, nRstart:nRstop, 4) )
      bytes_allocated = bytes_allocated + &
      &                 4*n_m_max*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      temp_Rloc(1:,nRstart:) => fields_container_RLoc(1:,nRstart:,1)
      us_Rloc(1:,nRstart:)   => fields_container_RLoc(1:,nRstart:,2)
      up_Rloc(1:,nRstart:)   => fields_container_RLoc(1:,nRstart:,3)
      om_Rloc(1:,nRstart:)   => fields_container_RLoc(1:,nRstart:,4)

      temp_Rloc(:,:) =zero
      us_Rloc(:,:)   =zero
      up_Rloc(:,:)   =zero
      om_Rloc(:,:)   =zero

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( fields_container_Rloc, fields_container_Mloc )
      deallocate( work_Mloc, dom_Mloc )
      if ( .not. l_cheb_coll ) deallocate( psi_hat_Mloc, temp_hat_Mloc )
      deallocate( dtemp_Mloc, psi_Mloc )

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
