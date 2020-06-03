module parallel_mod
   use mpimod

   implicit none

   integer :: rank, n_procs, ierr
   integer :: GLOBAL_COMM
   
   integer :: Color_Pizza, NProc_Pizza,Key_Pizza
   integer :: Comm_Pizza

   integer :: Color_LIBPFASST, NProc_LIBPFASST,Key_LIBPFASST
   integer :: Comm_LIBPFASST
   integer :: err
   
contains

   subroutine initialize_mpi

      call MPI_init(ierr)
      call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call MPI_comm_size(MPI_COMM_WORLD,n_procs, ierr)
      
!       NProc_LIBPFASST=1 
!       if (rank.EQ.0) then
!      
!       print*, "-------------------------------------------------------------"
!       print*, "Temporary/test value for NProc_LIBPFASST=",NProc_LIBPFASST
!       print*, "NProc_LIBPFASST to be chosen from input in actual simulations"
!       print*, "-------------------------------------------------------------"
!       endif

!       call mpi_split_time_parallel


   end subroutine initialize_mpi
   
!-----------------------------------------------------------------------------
   
   subroutine mpi_split_time_parallel(NProc_time)
!    use namelists, only: NProc_time
   integer :: i,j
   integer :: NProc_time
        NProc_LIBPFASST=NProc_time
        
        NProc_Pizza = n_procs/NProc_LIBPFASST ! this must be a integer 
        if (mod(n_procs,NProc_LIBPFASST).NE.0) then
            print*,"n_procs must be a multiple of NProc_LIBPFASST"
        endif
        if (rank.EQ.0) then
           print*, "-------------------------------------------------------------"
           print*, "NProc_Pizza =",NProc_Pizza
           print*, "NProc_time  =",NProc_LIBPFASST
           print*, "-------------------------------------------------------------"
        endif

        
        Key_Pizza    = mod(rank,NProc_Pizza)
        
        Color_LIBPFASST = Key_Pizza
        Key_LIBPFASST    = rank/NProc_Pizza 
        Color_Pizza = Key_LIBPFASST

        call MPI_Comm_split( MPI_COMM_WORLD, Color_Pizza, Key_Pizza, Comm_Pizza,err)



! 
!         if (rank.EQ.0) then
!            print*, "-------------------------------------------------------------"
!            print*, "Key_LIBPFASST and Color_LIBPFASST before second com split"
!            print*, "-------------------------------------------------------------"
!         endif
! 
!         do i=0,n_procs-1
!            if (rank.EQ.i) then
!               print*,"rank=", rank, "Key_LIBPFASST=",Key_LIBPFASST, "Color_LIBPFASST=",Color_LIBPFASST, "Comm_LIBPFASST=",Comm_LIBPFASST
!            endif
!            call MPI_Barrier(MPI_COMM_WORLD, ierr)
!         enddo 
        
        call MPI_Comm_split( MPI_COMM_WORLD, Color_LIBPFASST, Key_LIBPFASST, Comm_LIBPFASST,err)

        !-------------------------------------------------------------
        ! test prints, to be deleted after testing!
        !-------------------------------------------------------------
        if (rank.EQ.0) then
           print*, "-------------------------------------------------------------"
           print*, "Rank and Comm pizza"
        endif
        

        do i=0,n_procs-1
           if (rank.EQ.i) then
              print*,"rank=", rank, "Key_Pizza=",Key_Pizza, "Color_Pizza=", Color_Pizza, "Comm_Pizza=",Comm_Pizza
           endif
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
        enddo
        
        if (rank.EQ.0) then
           print*, "-------------------------------------------------------------"
        endif
        
        
        
        if (rank.EQ.0) then
           print*, "-------------------------------------------------------------"
           print*, "Rank and Comm LibPFASST"
        endif
        
        do i=0,n_procs-1
           if (rank.EQ.i) then
              print*,"rank=", rank, "Key_LIBPFASST=",Key_LIBPFASST, "Color_LIBPFASST=",Color_LIBPFASST, "Comm_LIBPFASST=",Comm_LIBPFASST
           endif
           call MPI_Barrier(MPI_COMM_WORLD, ierr)
        enddo
        
        if (rank.EQ.0) then
           print*, "-------------------------------------------------------------"
        endif
        
        
        
        
        !-------------------------------------------------------------
        ! End of test prints, to be deleted after testing!
        !-------------------------------------------------------------

   end subroutine mpi_split_time_parallel

!-----------------------------------------------------------------------------
   subroutine finalize_mpi

      call MPI_finalize(ierr)

   end subroutine finalize_mpi
!-----------------------------------------------------------------------------
end module parallel_mod
