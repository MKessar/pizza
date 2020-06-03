module blocking

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, n_m_max, m2idx
   use parallel_mod!, only: NProc_Pizza, rank

   implicit none

   private

   type :: load
      integer :: nStart
      integer :: nStop
      integer :: n_per_rank
      integer :: n_points
   end type load

   type(load), public, allocatable :: radial_balance(:)
   type(load), public, allocatable :: m_balance(:)

   integer, public :: nRstart, nRstop, nR_per_rank
   integer, public :: nMstart, nMstop, nm_per_rank
   logical, public :: l_rank_has_m0

   public :: set_mpi_domains, destroy_mpi_domains

contains

   subroutine set_mpi_domains

      integer :: idx

      allocate ( radial_balance(0:NProc_Pizza-1) )
      allocate ( m_balance(0:NProc_Pizza-1) )

      call getBlocks(radial_balance, n_r_max, NProc_Pizza)
      nRstart = radial_balance(Key_Pizza)%nStart
      nRstop = radial_balance(Key_Pizza)%nStop
      nR_per_rank = radial_balance(Key_Pizza)%n_per_rank
      call getBlocks(m_balance, n_m_max, NProc_Pizza)
      nMstart = m_balance(Key_Pizza)%nStart
      nMstop = m_balance(Key_Pizza)%nStop
      nm_per_rank = m_balance(Key_Pizza)%n_per_rank

      bytes_allocated = bytes_allocated+8*NProc_Pizza*SIZEOF_INTEGER

      idx = m2idx(0)

      if ( idx>=nMstart .and. idx<=nMstop ) then
         l_rank_has_m0=.true.
      else
         l_rank_has_m0=.false.
      end if

   end subroutine set_mpi_domains
!------------------------------------------------------------------------------
   subroutine destroy_mpi_domains

      deallocate( m_balance, radial_balance )

   end subroutine destroy_mpi_domains
!------------------------------------------------------------------------------
   subroutine getBlocks(bal, n_points, NProc_Pizza)

      type(load), intent(inout) :: bal(0:)
      integer, intent(in) :: NProc_Pizza
      integer, intent(in) :: n_points

      integer :: n_points_loc, check, p
      
      n_points_loc = n_points/NProc_Pizza

      check = mod(n_points,NProc_Pizza)-1

      bal(0)%nStart = 1

      do p =0, NProc_Pizza-1
         if ( p /= 0 ) bal(p)%nStart=bal(p-1)%nStop+1
         bal(p)%n_per_rank=n_points_loc
         if ( p <= check ) then
            bal(p)%n_per_rank=n_points_loc+1
         end if
         bal(p)%nStop=bal(p)%nStart+bal(p)%n_per_rank-1
         bal(p)%n_points=n_points
      end do

   end subroutine getBlocks
!------------------------------------------------------------------------------
end module blocking
