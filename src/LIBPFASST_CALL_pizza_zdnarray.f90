!
! This file is part of LIBPFASST.
!
!> Example of using LIBPFASST.
!!
!!  This program solves the 1-d advection diffusion problem on a periodic domain

!>  The main program here just initializes mpi, calls the solver and then finalizes mpi
module pizza_zdnarray 
!   use pf_mod_mpi
  use parallel_mod
  use iso_c_binding  
  use pf_mod_dtype
  use pf_mod_stop
  use precision_mod
  implicit none
  
  !>  Factory for making zndarray
  type, extends(pf_factory_t) :: pizza_zndarray_factory_t
   contains
     procedure :: create_single => pizza_zndarray_create_single
     procedure :: create_array => pizza_zndarray_create_array
     procedure :: destroy_single => pizza_zndarray_destroy_single
     procedure :: destroy_array => pizza_zndarray_destroy_array
  end type pizza_zndarray_factory_t
  
  !>  Complex N-dimensional array type,  extends the abstract encap type
  type, extends(pf_encap_t) :: pizza_zndarray_t
     integer :: ndim
     integer,    allocatable :: arr_shape(:)     
     complex(pfdp), allocatable :: flatarray(:)
   contains
     procedure :: setval => pizza_zndarray_setval
     procedure :: copy => pizza_zndarray_copy
     procedure :: norm => pizza_zndarray_norm
     procedure :: pack => pizza_zndarray_pack
     procedure :: unpack => pizza_zndarray_unpack
     procedure :: axpy => pizza_zndarray_axpy
     procedure :: eprint => pizza_zndarray_eprint
     procedure, private  :: get_array_1d  ,get_array_2d,get_array_3d,get_array_4d
     generic :: get_array => get_array_1d ,get_array_2d,get_array_3d,get_array_4d
  end type pizza_zndarray_t
  
contains
  
  function cast_as_pizza_zndarray(encap_polymorph) result(pizza_zndarray_obj)
    class(pf_encap_t), intent(in), target :: encap_polymorph
    type(pizza_zndarray_t), pointer :: pizza_zndarray_obj
    
    select type(encap_polymorph)
    type is (pizza_zndarray_t)
       pizza_zndarray_obj => encap_polymorph
    end select
  end function cast_as_pizza_zndarray
  
  !> Allocates complex ndarray
  subroutine pizza_zndarray_build(q, shape_in)
    class(pf_encap_t), intent(inout) :: q
    !type(pf_zndarray_t), intent(inout) :: q    
    integer,           intent(in   ) :: shape_in(:)
    
    integer :: ierr
    select type (q)
    class is (pizza_zndarray_t)
       allocate(q%arr_shape(SIZE(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                         
       allocate(q%flatarray(product(shape_in)),stat=ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                         
       q%ndim   = SIZE(shape_in)
       q%arr_shape = shape_in
       q%flatarray = cmplx(66666666.0, -66666666.0,pfdp)
    end select
  end subroutine pizza_zndarray_build
  
  !> Wrapper routine for allocation of a single zndarray type array
  subroutine pizza_zndarray_create_single(this, x, level_index,  lev_shape)
    class(pizza_zndarray_factory_t), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x
    integer, intent(in) :: level_index
    integer, intent(in) :: lev_shape(:)

    integer :: ierr
    allocate(pizza_zndarray_t::x,stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)                      
    call pizza_zndarray_build(x, lev_shape)

  end subroutine pizza_zndarray_create_single

  !> Wrapper routine for looped allocation of many zndarray type arrays
  subroutine pizza_zndarray_create_array(this, x, n, level_index,  lev_shape)
    class(pizza_zndarray_factory_t), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x(:)
    integer, intent(in) :: n
    integer, intent(in) :: level_index
    integer, intent(in) :: lev_shape(:)
    integer :: i,ierr

    allocate(pizza_zndarray_t::x(n),stat=ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'allocate fail, error=',ierr)    
    
    do i = 1, n
!       allocate(pf_zndarray_t::x(i),stat=ierr)
       call pizza_zndarray_build(x(i), lev_shape)
    end do

  end subroutine pizza_zndarray_create_array

  !>  Subroutine to destroy array
  subroutine pizza_zndarray_destroy(encap)
    class(pf_encap_t), intent(inout) :: encap
    type(pizza_zndarray_t), pointer :: pizza_zndarray_obj
    
    pizza_zndarray_obj => cast_as_pizza_zndarray(encap)
    deallocate(pizza_zndarray_obj%arr_shape)
    deallocate(pizza_zndarray_obj%flatarray)
    nullify(pizza_zndarray_obj)

  end subroutine pizza_zndarray_destroy

  !> Subroutine to destroy an single array
  subroutine pizza_zndarray_destroy_single(this, x)
    class(pizza_zndarray_factory_t), intent(inout) :: this
    class(pf_encap_t), intent(inout), allocatable :: x

    select type (x)
    class is (pizza_zndarray_t)
       deallocate(x%arr_shape)
       deallocate(x%flatarray)
    end select
    deallocate(x)

  end subroutine pizza_zndarray_destroy_single

  !> Wrapper routine for looped allocation of many zndarray type arrays
  subroutine pizza_zndarray_destroy_array(this, x)
    class(pizza_zndarray_factory_t), intent(inout)  :: this
    class(pf_encap_t), intent(inout),allocatable :: x(:)
    integer :: i

    select type(x)
    class is (pizza_zndarray_t)
       do i = 1,SIZE(x)
          deallocate(x(i)%arr_shape)
          deallocate(x(i)%flatarray)
       end do
    end select
    deallocate(x)
    
  end subroutine pizza_zndarray_destroy_array

  !>  The following are the base subroutines that all encapsulations must provide
  !!

  !> Subroutine to set array to a scalar  value
  subroutine pizza_zndarray_setval(this, val, flags)
    class(pizza_zndarray_t), intent(inout) :: this
    real(pfdp), intent(in) :: val
    integer, intent(in), optional :: flags
    complex(pfdp) :: zval

    zval = cmplx(val, 0.0, pfdp)
    this%flatarray = zval
  end subroutine pizza_zndarray_setval

  !> Subroutine to copy an array
  subroutine pizza_zndarray_copy(this, src, flags)
    class(pizza_zndarray_t), intent(inout) :: this
    class(pf_encap_t), intent(in) :: src
    integer, intent(in), optional :: flags
    class(pizza_zndarray_t), pointer :: pizza_zndarray_src

    pizza_zndarray_src => cast_as_pizza_zndarray(src)

    this%flatarray =  pizza_zndarray_src%flatarray
  end subroutine pizza_zndarray_copy

  !> Subroutine to pack an array into a flat array for sending
  subroutine pizza_zndarray_pack(this, z,flags)
    class(pizza_zndarray_t), intent(in) :: this
    real(pfdp), intent(out) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: i

    do i = 1,product(this%arr_shape)
       z(2*i-1) = REAL(this%flatarray(i))
       z(2*i)    = AIMAG(this%flatarray(i))
    end do
  end subroutine pizza_zndarray_pack

  !> Subroutine to unpack to a flatarray after receiving
  subroutine pizza_zndarray_unpack(this, z,flags)
    class(pizza_zndarray_t), intent(inout) :: this
    real(pfdp), intent(in) :: z(:)
    integer,     intent(in   ), optional :: flags
    integer :: i

    do i = 1,product(this%arr_shape)
       this%flatarray(i) = cmplx(z(2*i-1), z(2*i), pfdp)
    enddo
  end subroutine pizza_zndarray_unpack

  ! Compute norm of array

  function pizza_zndarray_norm(this,flags) result (norm)
    class(pizza_zndarray_t), intent(in) :: this
    integer,     intent(in   ), optional :: flags
    real(pfdp) :: norm
    real(pfdp) :: norm_local

    norm_local = maxval(abs(this%flatarray))
    
      call MPI_ALLREDUCE(norm_local, norm, 1, MPI_DEF_REAL, &
           &          MPI_MAX, Comm_Pizza, ierr)
  end function pizza_zndarray_norm
  ! Compute y = a x + y where a is a scalar and x and y are arrays
  subroutine pizza_zndarray_axpy(this, a, x, flags)
    class(pizza_zndarray_t), intent(inout) :: this
    class(pf_encap_t), intent(in) :: x
    real(pfdp), intent(in) :: a
    integer, intent(in), optional :: flags
    class(pizza_zndarray_t), pointer :: pizza_zndarray_obj

    if (a .eq. 0.0_pfdp) return
    
    pizza_zndarray_obj => cast_as_pizza_zndarray(x)
    if (a .eq. 1.0_pfdp) then
       this%flatarray = pizza_zndarray_obj%flatarray + this%flatarray
    elseif (a .eq. -1.0_pfdp) then
       this%flatarray = -pizza_zndarray_obj%flatarray + this%flatarray
    else
       this%flatarray = a*pizza_zndarray_obj%flatarray + this%flatarray
    end if
    
  end subroutine pizza_zndarray_axpy

  subroutine pizza_zndarray_eprint(this,flags)
    class(pizza_zndarray_t), intent(inout) :: this
    integer,           intent(in   ), optional :: flags

    print*, this%flatarray(1)
  end subroutine pizza_zndarray_eprint


  !>  Helper function to return the array part, these are called with get_array
  subroutine get_array_1d(this,r,flags) 
    class(pizza_zndarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:)
    integer,    intent(in   ), optional :: flags
    r => this%flatarray
  end subroutine get_array_1d
  subroutine get_array_2d(this,r,flags) 
    class(pizza_zndarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2)) => this%flatarray
  end subroutine get_array_2d
  subroutine get_array_3d(this,r,flags) 
    class(pizza_zndarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3)) => this%flatarray
  end subroutine get_array_3d
  subroutine get_array_4d(this,r,flags) 
    class(pizza_zndarray_t), target, intent(in) :: this
    complex(pfdp), pointer, intent(inout) :: r(:,:,:,:)
    integer,    intent(in   ), optional :: flags
    r(1:this%arr_shape(1),1:this%arr_shape(2),1:this%arr_shape(3),1:this%arr_shape(4)) => this%flatarray
  end subroutine get_array_4d

  !>  Helper function to return the array part
  function get_array1d(x,flags) result(r)
    class(pf_encap_t), target, intent(in) :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:)
    select type (x)
    type is (pizza_zndarray_t)
       r => x%flatarray
    end select
  end function get_array1d
  

  function get_array2d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:)

    select type (x)
    type is (pizza_zndarray_t)
       r(1:x%arr_shape(1),1:x%arr_shape(2)) => x%flatarray
    end select
  end function get_array2d
  

  function get_array3d(x,flags) result(r)
    class(pf_encap_t), intent(in),target :: x
    integer,           intent(in   ), optional :: flags
    complex(pfdp), pointer :: r(:,:,:)

    select type (x)
    type is (pizza_zndarray_t)
       r(1:x%arr_shape(1),1:x%arr_shape(2),1:x%arr_shape(3)) => x%flatarray
    end select
  end function get_array3d

end module pizza_zdnarray
