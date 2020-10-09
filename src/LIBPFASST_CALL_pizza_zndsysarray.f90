!!  N-dimensional complex system of arrays encapsulation
!
! This file is part of LIBPFASST.
!

!> !> Module to define and encapsulation for a system of N-ndimensional complex arrays 
!!
!! When a new solution is created by a PFASST level, this encapsulation
!! uses the levels 'lev_shape' attribute to create a new multi-component array with that
!! shape.  Thus, the 'lev_shape' attributes of the PFASST levels should be
!! set appropriately.  The last component of lev_shape is the number of components in the system
!!
!! For example, before calling pf_pfasst_run we can
!! set the arr_shape of the coarsest level by doing:
!!
!!   allocate(pf%levels(1)%lev_shape(3))
!!   pf%levels(1)%lev_shape = [ nx, ny, 3 ]
!!
!! Which would imply that a 3 component system of two-dimensional solutions.
!!
!! The helper routines get_array1d, get_array2d, get_array3d, etc can be used to
!! extract pointers to a component of  encapsulated system without
!! performing any copies, or return a pointer to the whole array.
!!
module pizza_zndsysarray
  use parallel_mod
  use precision_mod
  use iso_c_binding
  use pf_mod_dtype
  use pf_mod_stop
  implicit none

  !>  Type to create and destroy the arrays
  type, extends(pf_factory_t) :: pizza_zndsysarray_factory_t
   contains
     procedure :: create_single  => pizza_zndsysarray_create_single
     procedure :: create_array  => pizza_zndsysarray_create_array
     procedure :: destroy_single => pizza_zndsysarray_destroy_single
     procedure :: destroy_array => pizza_zndsysarray_destroy_array
  end type pizza_zndsysarray_factory_t

  !>  Type to extend the abstract encap and set procedure pointers
  type, extends(pf_encap_t) :: pizza_zndsysarray_t
     integer             :: ndim    !  The spatial dimension of each component in system
     integer             :: ncomp  !  The number of components in the system
     integer             :: ndof   !  The number of variables in each component
     integer,    allocatable :: arr_shape(:)
     complex(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => pizza_zndsysarray_setval
     procedure :: copy => pizza_zndsysarray_copy
     procedure :: norm => pizza_zndsysarray_norm
     procedure :: pack => pizza_zndsysarray_pack
     procedure :: unpack => pizza_zndsysarray_unpack
     procedure :: axpy => pizza_zndsysarray_axpy
     procedure :: eprint => pizza_zndsysarray_eprint
     procedure, private  :: get_sysarray_1d, get_sysarray_2d,get_sysarray_3d,get_sysarray_4d
     generic :: get_sysarray => get_sysarray_1d ,get_sysarray_2d,get_sysarray_3d,get_sysarray_4d
     procedure, private  :: get_array_2d
!      procedure, private  :: get_array_1d  ,get_array_2d,get_array_3d,get_array_4d
!      generic :: get_array => get_array_1d ,get_array_2d,get_array_3d,get_array_4d
     generic :: get_array => get_array_2d
  end type pizza_zndsysarray_t

contains
  !>  Subroutine to allocate the array and set the size parameters
  subroutine pizza_zndsysarray_build(q, arr_shape)
    class(pf_encap_t), intent(inout) :: q
    integer,           intent(in   ) :: arr_shape(:)

    integer :: ierr

    select type (q)
    class is (pizza_zndsysarray_t)
       allocate(q%arr_shape(SIZE(arr_shape)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                      
       
       q%ndim   = SIZE(arr_shape)-1
       q%ncomp = arr_shape(q%ndim+1)
       q%ndof = product(arr_shape(1:q%ndim))
       q%arr_shape = arr_shape

!        allocate(q%flatarray(product(arr_shape)),stat=ierr)
       allocate(q%flatarray(product(arr_shape(1:q%ndim))*(q%ncomp)),stat=ierr)
!        allocate(q%flatarray(product(arr_shape(1:q%ndim))*(q%ncomp-1)+arr_shape(2)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
       
    end select
  end subroutine pizza_zndsysarray_build

  !> Subroutine to  create a single array
  subroutine pizza_zndsysarray_create_single(this, x, level_index, lev_shape)
    class(pizza_zndsysarray_factory_t), intent(inout)     :: this
    class(pf_encap_t),      intent(inout), allocatable :: x
    integer,                intent(in   )              :: level_index
    integer,                intent(in   )              :: lev_shape(:)
    integer :: i,ierr
    allocate(pizza_zndsysarray_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                             
    call pizza_zndsysarray_build(x, lev_shape)
  end subroutine pizza_zndsysarray_create_single

  !> Subroutine to create an array of arrays
  subroutine pizza_zndsysarray_create_array(this, x, n, level_index,  lev_shape)
    class(pizza_zndsysarray_factory_t), intent(inout)     :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer,                intent(in   )              :: n
    integer,                intent(in   )              ::  level_index
    integer,                intent(in   )              ::  lev_shape(:)
    integer :: i,ierr
    allocate(pizza_zndsysarray_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                                 
    do i = 1, n
       call pizza_zndsysarray_build(x(i), lev_shape)
    end do
  end subroutine pizza_zndsysarray_create_array

  !>  Subroutine to destroy array
  subroutine pizza_zndsysarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pizza_zndsysarray_t), pointer :: pizza_zndsysarray_obj

    pizza_zndsysarray_obj => cast_as_pizza_zndsysarray(encap)

    deallocate(pizza_zndsysarray_obj%arr_shape)
    deallocate(pizza_zndsysarray_obj%flatarray)

    nullify(pizza_zndsysarray_obj)

  end subroutine pizza_zndsysarray_destroy

  !> Subroutine to destroy an single array
  subroutine pizza_zndsysarray_destroy_single(this, x)
    class(pizza_zndsysarray_factory_t), intent(inout)     :: this
    class(pf_encap_t),      intent(inout), allocatable :: x

    select type (x)
    class is (pizza_zndsysarray_t)
       deallocate(x%arr_shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)
  end subroutine pizza_zndsysarray_destroy_single


  !> Subroutine to destroy an array of arrays
  subroutine pizza_zndsysarray_destroy_array(this, x)
    class(pizza_zndsysarray_factory_t), intent(inout)     :: this
    class(pf_encap_t),      intent(inout), allocatable :: x(:)
    integer                                            :: i

    select type(x)
    class is (pizza_zndsysarray_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%arr_shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
  end subroutine pizza_zndsysarray_destroy_array


  !>  The following are the base subroutines that all encapsulations must provide
  !!
  
  !> Subroutine to set array to a scalare  value.
  subroutine pizza_zndsysarray_setval(this, val, flags)
    class(pizza_zndsysarray_t), intent(inout)           :: this
    real(pfdp),     intent(in   )           :: val
    integer,        intent(in   ), optional :: flags
    this%flatarray = val
  end subroutine pizza_zndsysarray_setval

  !> Subroutine to copy an array
  subroutine pizza_zndsysarray_copy(this, src, flags)
    class(pizza_zndsysarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: src
    integer,           intent(in   ), optional :: flags
    select type(src)
    type is (pizza_zndsysarray_t)
       this%flatarray = src%flatarray
    class default
       call pf_stop(__FILE__,__LINE__,'unsupported type ')
    end select
  end subroutine pizza_zndsysarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine pizza_zndsysarray_pack(this, z, flags)
    class(pizza_zndsysarray_t), intent(in   ) :: this
    real(pfdp),     intent(  out) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: ntot
    ntot = this%ndof*this%ncomp
    
    z(1:ntot) = REAL(this%flatarray,pfdp)
    z(ntot+1:2*ntot) = AIMAG(this%flatarray)
  end subroutine pizza_zndsysarray_pack

  !> Subroutine to unpack a flatarray after receiving
  subroutine pizza_zndsysarray_unpack(this, z, flags)
    class(pizza_zndsysarray_t), intent(inout) :: this
    real(pfdp),     intent(in   ) :: z(:)
    integer,     intent(in   ), optional :: flags

    integer :: ntot    
    ntot = this%ndof*this%ncomp
    this%flatarray = z(1:ntot)
!    this%flatarray = this%flatarray + cmplx(0.0,1.0,pfdp)*z(ntot+1:2*ntot)
    this%flatarray =  cmplx(z(1:ntot),z(ntot+1:2*ntot),pfdp)        
  end subroutine pizza_zndsysarray_unpack

  !> Subroutine to define the norm of the array (here the max norm)
  function pizza_zndsysarray_norm(this, flags) result (norm)
    class(pizza_zndsysarray_t), intent(in   ) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    real(pfdp) :: norm_local
    
    norm_local = maxval(abs(this%flatarray))
    
      call MPI_ALLREDUCE(norm_local, norm, 1, MPI_DEF_REAL, &
           &          MPI_MAX, Comm_Pizza, ierr)    
  end function pizza_zndsysarray_norm

  !> Subroutine to compute y = a x + y where a is a scalar and x and y are arrays
  subroutine pizza_zndsysarray_axpy(this, a, x, flags)
    class(pizza_zndsysarray_t),    intent(inout)           :: this
    class(pf_encap_t), intent(in   )           :: x
    real(pfdp),        intent(in   )           :: a
    integer,           intent(in   ), optional :: flags
    if (a .eq. 0.0_pfdp) return
    select type(x)
    type is (pizza_zndsysarray_t)
       if (a .eq. 1.0_pfdp) then
          this%flatarray = x%flatarray + this%flatarray
       elseif (a .eq. -1.0_pfdp) then
          this%flatarray = -x%flatarray + this%flatarray
       else
          this%flatarray = a*x%flatarray + this%flatarray
       end if
    class default
       call pf_stop(__FILE__,__LINE__,'unsupported type ')
    end select
  end subroutine pizza_zndsysarray_axpy

  !>  Subroutine to print the array to the screen (mainly for debugging purposes)
  subroutine pizza_zndsysarray_eprint(this,flags)
    class(pizza_zndsysarray_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags
    !  Just print the first few values
        print *, this%flatarray(1:10)
    !print *, this%flatarray
  end subroutine pizza_zndsysarray_eprint


  function cast_as_pizza_zndsysarray(encap_polymorph) result(pizza_zndsysarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pizza_zndsysarray_t), pointer :: pizza_zndsysarray_obj
    
    select type(encap_polymorph)
    type is (pizza_zndsysarray_t)
       pizza_zndsysarray_obj => encap_polymorph
    end select
  end function cast_as_pizza_zndsysarray

  !>  Helper function to return the array part, these are called with get_array
  subroutine get_sysarray_1d(this,r,flags) 
    class(pizza_zndsysarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:)
    integer,    intent(in   ), optional :: flags
    r => this%flatarray
  end subroutine get_sysarray_1d
  subroutine get_sysarray_2d(this,r,flags) 
    class(pizza_zndsysarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2)) => this%flatarray
  end subroutine get_sysarray_2d
  subroutine get_sysarray_3d(this,r,flags) 
    class(pizza_zndsysarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3)) => this%flatarray
  end subroutine get_sysarray_3d
  subroutine get_sysarray_4d(this,r,flags) 
    class(pizza_zndsysarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3),1:this%arr_shape(4)) => this%flatarray
  end subroutine get_sysarray_4d
  
  !>  Helper function to return the array part or whole
  function get_array1d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:)
    select type (x)
    type is (pizza_zndsysarray_t)
       if (n .eq. 0) then  !  Return the whole array
          if (x%ncomp .eq. 1 .and. x%ndim .eq. 1) then
             r => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       else                !  Return the nth component
          if (x%ndim .eq. 1) then
             r => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       end if
    end select
  end function get_array1d
  

  function get_array2d(x,n,flags) result(r)
    class(pf_encap_t), target,intent(in) :: x
    integer, intent(in) :: n    
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:)

    
!     print*,"get_array2d n=",n
    select type (x)
    type is (pizza_zndsysarray_t)
!     print*,"get_array2d 1"
       if (n .eq. 0) then
!     print*,"get_array2d 2"
          if (x%ndim .eq. 1) then
             r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
          end if
       else              !  Return the nth component
!     print*,"get_array2d 3"
          if (x%ndim .eq. 2) then
             r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
          else
             call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
          end if
!     elseif  (n .LT. 4) then              !  Return the nth component
! !     print*,"get_array2d 3"
!        if (x%ndim .eq. 2) then
!           r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!        end if
!     elseif  (n .eq. 4) then              !  Return the nth component
! !     print*,"get_array2d 4"
!        if (x%ndim .eq. 2) then
!           r(1:1,1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*(n-1)+x%arr_shape(2))
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!        end if
       endif
    end select
!     print*,"get_array2d 5"
!         print*,"size r",SIZE(r)

  end function get_array2d
  
! 
!   function get_array3d(x,n,flags) result(r)
!     class(pf_encap_t), target,intent(in) :: x
!     integer, intent(in) :: n
!     integer,           intent(in   ), optional :: flags
!     complex(pfdp), pointer :: r(:,:,:)
! 
!     select type (x)
!     type is (pizza_zndsysarray_t)
!        if (n .eq. 0) then        !  Return the whole array  
!           if (x%ndim .eq. 2) then
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!           end if
!        else                     !  Return the nth component
!           if (x%ndim .eq. 3) then
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
!           end if
!        end if
!     end select
!   end function get_array3d
! 
!   function get_array4d(x,n,flags) result(r)
!     class(pf_encap_t), target,intent(in) :: x
!     integer, intent(in) :: n
!     integer,           intent(in   ), optional :: flags
!     complex(pfdp), pointer :: r(:,:,:,:)
! 
!     select type (x)
!     type is (pizza_zndsysarray_t)
!        if (n .eq. 0) then        !  Return the whole array
!           if (x%ndim .eq. 3) then
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
!           end if
!        else                       !  Return the nth component
!           if (x%ndim .eq. 4) then   
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 4. ndim=',x%ndim)
!           end if
!        end if
!     end select
!   end function get_array4d

  !>  Helper function to return the array part or whole
!   subroutine get_array_1d(x,r,n,flags) 
!     class(pizza_zndsysarray_t), target,intent(in) :: x
!     integer, intent(in) :: n
!     integer,           intent(in   ), optional :: flags
!     complex(pfdp), pointer :: r(:)
!     if (n .eq. 0) then  !  Return the whole array
!        if (x%ncomp .eq. 1 .and. x%ndim .eq. 1) then
!           r => x%flatarray
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
!           end if
!        else                !  Return the nth component
!           if (x%ndim .eq. 1) then
!              r => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
!           end if
!        end if
! 
!   end subroutine get_array_1d
!   

  subroutine get_array_2d(x,r,n,flags) 
    class(pizza_zndsysarray_t), target,intent(in) :: x
    integer, intent(in) :: n    
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:)

    if (n .eq. 0) then
       if (x%ndim .eq. 1) then
          r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray
       else
          call pf_stop(__FILE__,__LINE__,'bad dimension, must be 1. ndim=',x%ndim)
       end if
    else              !  Return the nth component
       if (x%ndim .eq. 2) then
          r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
       else
          call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
       end if
!     elseif  (n .LT. 4) then              !  Return the nth component
!        if (x%ndim .eq. 2) then
!           r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!        end if
!     elseif  (n .eq. 4) then              !  Return the nth component
!        if (x%ndim .eq. 2) then
!           r(1:1,1:x%arr_shape(2)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*(n-1)+x%arr_shape(2))
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!        end if
        
    endif
   end subroutine get_array_2d
  
! 
!   subroutine get_array_3d(x,r,n,flags) 
!     class(pizza_zndsysarray_t), target,intent(in) :: x
!     integer, intent(in) :: n
!     integer,           intent(in   ), optional :: flags
!     complex(pfdp), pointer :: r(:,:,:)
! 
!     if (n .eq. 0) then        !  Return the whole array  
!        if (x%ndim .eq. 2) then
!           r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 2. ndim=',x%ndim)
!        end if
!     else                     !  Return the nth component
!        if (x%ndim .eq. 3) then
!           r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!        else
!           call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
!        end if
!     end if
! 
!   end subroutine get_array_3d
! 
!   subroutine get_array_4d(x,r,n,flags) 
!     class(pizza_zndsysarray_t), target,intent(in) :: x
!     integer, intent(in) :: n
!     integer,           intent(in   ), optional :: flags
!     complex(pfdp), pointer :: r(:,:,:,:)
! 
!     select type (x)
!     type is (pizza_zndsysarray_t)
!        if (n .eq. 0) then        !  Return the whole array
!           if (x%ndim .eq. 3) then
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 3. ndim=',x%ndim)
!           end if
!        else                       !  Return the nth component
!           if (x%ndim .eq. 4) then   
!              r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3),1:x%arr_shape(4)) => x%flatarray(x%ndof*(n-1)+1:x%ndof*n)
!           else
!              call pf_stop(__FILE__,__LINE__,'bad dimension, must be 4. ndim=',x%ndim)
!           end if
!        end if
!     end select
!   end subroutine get_array_4d
  
end module pizza_zndsysarray
