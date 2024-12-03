! # 1 "radiation/radiation_matrix.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_matrix.f90"
! radiation_matrix.f90 - spartacus matrix operations
!
! (c) copyright 2014- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! author:  robin hogan
! email:   r.j.hogan@ecmwf.int
!
! modifications
!   2018-10-15  r. hogan  added fast_expm_exchange_[23]
!
! this module provides the neccessary mathematical functions for the
! spartacus radiation scheme: matrix multiplication, matrix solvers
! and matrix exponentiation, but (a) multiple matrices are operated on
! at once with array access indended to facilitate vectorization, and
! (b) optimization for 2x2 and 3x3 matrices.  there is probably
! considerable scope for further optimization. note that this module
! is not used by the mcica solver.

module radiation_matrix

  use parkind1, only : jprb

  implicit none
  public

  ! codes to describe sparseness pattern, where the shortwave
  ! pattern is of the form:
  ! (x x x)
  ! (x x x)
  ! (0 0 x)
  ! where each element may itself be a square matrix.  
  integer, parameter :: imatrixpatterndense     = 0
  integer, parameter :: imatrixpatternshortwave = 1

  public  :: mat_x_vec, singlemat_x_vec, mat_x_mat, &
       &     singlemat_x_mat, mat_x_singlemat, &
       &     identity_minus_mat_x_mat, solve_vec, solve_mat, expm, &
       &     fast_expm_exchange_2, fast_expm_exchange_3, &
       &     sparse_x_dense

  private :: solve_vec_2, solve_vec_3, solve_mat_2, &
       &     solve_mat_3, lu_factorization, lu_substitution, solve_mat_n, &
       &     diag_mat_right_divide_3

  interface fast_expm_exchange
    module procedure fast_expm_exchange_2, fast_expm_exchange_3
  end interface fast_expm_exchange

contains

  ! --- matrix-vector multiplication ---

  !---------------------------------------------------------------------
  ! treat a as n m-by-m square matrices (with the n dimension varying
  ! fastest) and b as n m-element vectors, and perform matrix-vector
  ! multiplications on first iend pairs
  function mat_x_vec(n,iend,m,a,b,do_top_left_only_in)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                   :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:) :: a
    real(jprb), intent(in), dimension(:,:)   :: b
    logical,    intent(in), optional         :: do_top_left_only_in
    real(jprb),             dimension(iend,m):: mat_x_vec

    integer :: j1, j2
    logical :: do_top_left_only

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',0,hook_handle)

    if (present(do_top_left_only_in)) then
      do_top_left_only = do_top_left_only_in
    else
      do_top_left_only = .false.
    end if

    ! array-wise assignment
    mat_x_vec = 0.0_jprb

    if (do_top_left_only) then
      mat_x_vec(1:iend,1) = a(1:iend,1,1)*b(1:iend,1)
    else
      do j1 = 1,m
        do j2 = 1,m
          mat_x_vec(1:iend,j1) = mat_x_vec(1:iend,j1) &
               &               + a(1:iend,j1,j2)*b(1:iend,j2)
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_vec',1,hook_handle)

  end function mat_x_vec


  !---------------------------------------------------------------------
  ! treat a as an m-by-m square matrix and b as n m-element vectors
  ! (with the n dimension varying fastest), and perform matrix-vector
  ! multiplications on first iend pairs
  function singlemat_x_vec(n,iend,m,a,b)

!    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                    :: n, m, iend
    real(jprb), intent(in), dimension(m,m)    :: a
    real(jprb), intent(in), dimension(:,:)    :: b
    real(jprb),             dimension(iend,m) :: singlemat_x_vec

    integer    :: j1, j2
!    real(jphook) :: hook_handle

!    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',0,hook_handle)

    ! array-wise assignment
    singlemat_x_vec = 0.0_jprb

    do j1 = 1,m
      do j2 = 1,m
        singlemat_x_vec(1:iend,j1) = singlemat_x_vec(1:iend,j1) &
             &                    + a(j1,j2)*b(1:iend,j2)
      end do
    end do

!    if (lhook) call dr_hook('radiation_matrix:single_mat_x_vec',1,hook_handle)

  end function singlemat_x_vec


  ! --- square matrix-matrix multiplication ---

  !---------------------------------------------------------------------
  ! treat a and b each as n m-by-m square matrices (with the n
  ! dimension varying fastest) and perform matrix multiplications on
  ! all n matrix pairs
  function mat_x_mat(n,iend,m,a,b,i_matrix_pattern)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    integer,    intent(in), optional            :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:)    :: a, b

    real(jprb),             dimension(iend,m,m) :: mat_x_mat
    integer    :: j1, j2, j3
    integer    :: mblock, m2block
    integer    :: i_actual_matrix_pattern
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = imatrixpatterndense
    end if

    ! array-wise assignment
    mat_x_mat = 0.0_jprb

    if (i_actual_matrix_pattern == imatrixpatternshortwave) then
      ! matrix has a sparsity pattern
      !     (c d e)
      ! a = (f g h)
      !     (0 0 i)
      mblock = m/3
      m2block = 2*mblock 
      ! do the top-left (c, d, f, g)
      do j2 = 1,m2block
        do j1 = 1,m2block
          do j3 = 1,m2block
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + a(1:iend,j1,j3)*b(1:iend,j3,j2)
          end do
        end do
      end do
      do j2 = m2block+1,m
        ! do the top-right (e & h)
        do j1 = 1,m2block
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + a(1:iend,j1,j3)*b(1:iend,j3,j2)
          end do
        end do
        ! do the bottom-right (i)
        do j1 = m2block+1,m
          do j3 = m2block+1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + a(1:iend,j1,j3)*b(1:iend,j3,j2)
          end do
        end do
      end do
    else
      ! ordinary dense matrix
      do j2 = 1,m
        do j1 = 1,m
          do j3 = 1,m
            mat_x_mat(1:iend,j1,j2) = mat_x_mat(1:iend,j1,j2) &
                 &                  + a(1:iend,j1,j3)*b(1:iend,j3,j2)
          end do
        end do
      end do
    end if

    if (lhook) call dr_hook('radiation_matrix:mat_x_mat',1,hook_handle)

  end function mat_x_mat


  !---------------------------------------------------------------------
  ! treat a as an m-by-m matrix and b as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function singlemat_x_mat(n,iend,m,a,b)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(m,m)      :: a
    real(jprb), intent(in), dimension(:,:,:)    :: b
    real(jprb),             dimension(iend,m,m) :: singlemat_x_mat

    integer    :: j1, j2, j3
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',0,hook_handle)

    ! array-wise assignment
    singlemat_x_mat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          singlemat_x_mat(1:iend,j1,j2) = singlemat_x_mat(1:iend,j1,j2) &
               &                        + a(j1,j3)*b(1:iend,j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:singlemat_x_mat',1,hook_handle)

  end function singlemat_x_mat


  !---------------------------------------------------------------------
  ! treat b as an m-by-m matrix and a as n m-by-m square matrices
  ! (with the n dimension varying fastest) and perform matrix
  ! multiplications on the first iend matrix pairs
  function mat_x_singlemat(n,iend,m,a,b)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                      :: n, m, iend
    real(jprb), intent(in), dimension(:,:,:)    :: a
    real(jprb), intent(in), dimension(m,m)      :: b

    real(jprb),             dimension(iend,m,m) :: mat_x_singlemat
    integer    :: j1, j2, j3
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',0,hook_handle)

    ! array-wise assignment
    mat_x_singlemat = 0.0_jprb

    do j2 = 1,m
      do j1 = 1,m
        do j3 = 1,m
          mat_x_singlemat(1:iend,j1,j2) = mat_x_singlemat(1:iend,j1,j2) &
               &                        + a(1:iend,j1,j3)*b(j3,j2)
        end do
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:mat_x_singlemat',1,hook_handle)

  end function mat_x_singlemat


  !---------------------------------------------------------------------
  ! compute i-a*b where i is the identity matrix and a & b are n
  ! m-by-m square matrices
  function identity_minus_mat_x_mat(n,iend,m,a,b,i_matrix_pattern)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)                   :: n, m, iend
    integer,    intent(in), optional         :: i_matrix_pattern
    real(jprb), intent(in), dimension(:,:,:) :: a, b
    real(jprb),             dimension(iend,m,m) :: identity_minus_mat_x_mat

    integer    :: j
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',0,hook_handle)

    if (present(i_matrix_pattern)) then
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,a,b,i_matrix_pattern)
    else
      identity_minus_mat_x_mat = mat_x_mat(n,iend,m,a,b)
    end if

    identity_minus_mat_x_mat = - identity_minus_mat_x_mat
    do j = 1,m
      identity_minus_mat_x_mat(1:iend,j,j) &
           &     = 1.0_jprb + identity_minus_mat_x_mat(1:iend,j,j)
    end do

    if (lhook) call dr_hook('radiation_matrix:identity_mat_x_mat',1,hook_handle)

  end function identity_minus_mat_x_mat


  
  !---------------------------------------------------------------------
  ! replacement for matmul in the case that the first matrix is sparse
  function sparse_x_dense(sparse, dense)

    real(jprb), intent(in) :: sparse(:,:), dense(:,:)
    real(jprb) :: sparse_x_dense(size(sparse,1),size(dense,2))

    integer :: j1, j2, j3 ! loop indices
    integer :: n1, n2, n3 ! array sizes

    n1 = size(sparse,1)
    n2 = size(sparse,2)
    n3 = size(dense,2)
    
    sparse_x_dense = 0.0_jprb
    do j2 = 1,n2
      do j1 = 1,n1
        if (sparse(j1,j2) /= 0.0_jprb) then
          sparse_x_dense(j1,:) = sparse_x_dense(j1,:) + sparse(j1,j2)*dense(j2,:)
        end if
      end do
    end do
    
  end function sparse_x_dense
  

  ! --- repeatedly square a matrix ---

  !---------------------------------------------------------------------
  ! square m-by-m matrix "a" nrepeat times. a will be corrupted by
  ! this function.
  function repeated_square(m,a,nrepeat,i_matrix_pattern)
    integer,    intent(in)           :: m, nrepeat
    real(jprb), intent(inout)        :: a(m,m)
    integer,    intent(in), optional :: i_matrix_pattern
    real(jprb)                       :: repeated_square(m,m)

    integer :: j1, j2, j3, j4
    integer :: mblock, m2block
    integer :: i_actual_matrix_pattern

    if (present(i_matrix_pattern)) then
      i_actual_matrix_pattern = i_matrix_pattern
    else
      i_actual_matrix_pattern = imatrixpatterndense
    end if

    if (i_actual_matrix_pattern == imatrixpatternshortwave) then
      ! matrix has a sparsity pattern
      !     (c d e)
      ! a = (f g h)
      !     (0 0 i)
      mblock = m/3
      m2block = 2*mblock
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        ! do the top-left (c, d, f & g)
        do j2 = 1,m2block
          do j1 = 1,m2block
            do j3 = 1,m2block
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + a(j1,j3)*a(j3,j2)
            end do
          end do
        end do
        do j2 = m2block+1, m
          ! do the top-right (e & h)
          do j1 = 1,m2block
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + a(j1,j3)*a(j3,j2)
            end do
          end do
          ! do the bottom-right (i)
          do j1 = m2block+1, m
            do j3 = m2block+1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + a(j1,j3)*a(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          a = repeated_square
        end if
      end do
    else
      ! ordinary dense matrix
      do j4 = 1,nrepeat
        repeated_square = 0.0_jprb
        do j2 = 1,m
          do j1 = 1,m
            do j3 = 1,m
              repeated_square(j1,j2) = repeated_square(j1,j2) &
                   &                 + a(j1,j3)*a(j3,j2)
            end do
          end do
        end do
        if (j4 < nrepeat) then
          a = repeated_square
        end if
      end do
    end if

  end function repeated_square


  ! --- solve linear equations ---

  !---------------------------------------------------------------------
  ! solve ax=b to obtain x.  version optimized for 2x2 matrices using
  ! cramer's method: "a" contains n 2x2 matrices and "b" contains n
  ! 2-element vectors; returns a^-1 b.
  pure subroutine solve_vec_2(n,iend,a,b,x)

    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  a(1:iend,1,1)*a(1:iend,2,2) &
         &                - a(1:iend,1,2)*a(1:iend,2,1))

    x(1:iend,1) = inv_det*(a(1:iend,2,2)*b(1:iend,1)-a(1:iend,1,2)*b(1:iend,2))
    x(1:iend,2) = inv_det*(a(1:iend,1,1)*b(1:iend,2)-a(1:iend,2,1)*b(1:iend,1))

  end subroutine solve_vec_2


  !---------------------------------------------------------------------
  ! solve ax=b to obtain x, i.e. the matrix right-hand-side version of
  ! solve_vec_2, with a, x and b all containing n 2x2 matrices;
  ! returns a^-1 b using cramer's method.
  pure subroutine solve_mat_2(n,iend,a,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(in)  :: b(:,:,:)
    real(jprb), intent(out) :: x(:,:,:)

    real(jprb) :: inv_det(iend)

    inv_det = 1.0_jprb / (  a(1:iend,1,1)*a(1:iend,2,2) &
         &                - a(1:iend,1,2)*a(1:iend,2,1))

    x(1:iend,1,1) = inv_det*( a(1:iend,2,2)*b(1:iend,1,1) &
         &                   -a(1:iend,1,2)*b(1:iend,2,1))
    x(1:iend,2,1) = inv_det*( a(1:iend,1,1)*b(1:iend,2,1) &
         &                   -a(1:iend,2,1)*b(1:iend,1,1))
    x(1:iend,1,2) = inv_det*( a(1:iend,2,2)*b(1:iend,1,2) &
         &                   -a(1:iend,1,2)*b(1:iend,2,2))
    x(1:iend,2,2) = inv_det*( a(1:iend,1,1)*b(1:iend,2,2) &
         &                   -a(1:iend,2,1)*b(1:iend,1,2))

  end subroutine solve_mat_2


  !---------------------------------------------------------------------
  ! solve ax=b optimized for 3x3 matrices, using lu
  ! factorization and substitution without pivoting.
  pure subroutine solve_vec_3(n,iend,a,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(in)  :: b(:,:)
    real(jprb), intent(out) :: x(:,:)

    real(jprb), dimension(iend) :: l21, l31, l32
    real(jprb), dimension(iend) :: u22, u23, u33
    real(jprb), dimension(iend) :: y2, y3

    ! some compilers unfortunately don't support assocate
    !    associate (u11 => a(:,1,1), u12 => a(:,1,2), u13 => a(1,3), &
    !         y1 => b(:,1), x1 => solve_vec3(:,1), &
    !         x2 => solve_vec3(:,2), x3 => solve_vec3(:,3))

    ! lu decomposition:
    !     ( 1        )   (u11 u12 u13)
    ! a = (l21  1    ) * (    u22 u23)
    !     (l31 l32  1)   (        u33)
    l21 = a(1:iend,2,1) / a(1:iend,1,1)
    l31 = a(1:iend,3,1) / a(1:iend,1,1)
    u22 = a(1:iend,2,2) - l21*a(1:iend,1,2)
    u23 = a(1:iend,2,3) - l21*a(1:iend,1,3)
    l32 =(a(1:iend,3,2) - l31*a(1:iend,1,2)) / u22
    u33 = a(1:iend,3,3) - l31*a(1:iend,1,3) - l32*u23

    ! solve ly = b by forward substitution
    y2 = b(1:iend,2) - l21*b(1:iend,1)
    y3 = b(1:iend,3) - l31*b(1:iend,1) - l32*y2

    ! solve ux = y by back substitution
    x(1:iend,3) = y3/u33
    x(1:iend,2) = (y2 - u23*x(1:iend,3)) / u22
    x(1:iend,1) = (b(1:iend,1) - a(1:iend,1,2)*x(1:iend,2) &
         &         - a(1:iend,1,3)*x(1:iend,3)) / a(1:iend,1,1)
    !    end associate

  end subroutine solve_vec_3


  !---------------------------------------------------------------------
  ! solve ax=b optimized for 3x3 matrices, using lu factorization and
  ! substitution with no pivoting.
  pure subroutine solve_mat_3(n,iend,a,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(in)  :: b(:,:,:)
    real(jprb), intent(out) :: x(:,:,:)

    real(jprb), dimension(iend) :: l21, l31, l32
    real(jprb), dimension(iend) :: u22, u23, u33
    real(jprb), dimension(iend) :: y2, y3

    integer :: j

    !    associate (u11 => a(:,1,1), u12 => a(:,1,2), u13 => a(1,3))
    ! lu decomposition:
    !     ( 1        )   (u11 u12 u13)
    ! a = (l21  1    ) * (    u22 u23)
    !     (l31 l32  1)   (        u33)
    l21 = a(1:iend,2,1) / a(1:iend,1,1)
    l31 = a(1:iend,3,1) / a(1:iend,1,1)
    u22 = a(1:iend,2,2) - l21*a(1:iend,1,2)
    u23 = a(1:iend,2,3) - l21*a(1:iend,1,3)
    l32 =(a(1:iend,3,2) - l31*a(1:iend,1,2)) / u22
    u33 = a(1:iend,3,3) - l31*a(1:iend,1,3) - l32*u23

    do j = 1,3
      ! solve ly = b(:,:,j) by forward substitution
      ! y1 = b(:,1,j)
      y2 = b(1:iend,2,j) - l21*b(1:iend,1,j)
      y3 = b(1:iend,3,j) - l31*b(1:iend,1,j) - l32*y2
      ! solve ux(:,:,j) = y by back substitution
      x(1:iend,3,j) = y3 / u33
      x(1:iend,2,j) = (y2 - u23*x(1:iend,3,j)) / u22
      x(1:iend,1,j) = (b(1:iend,1,j) - a(1:iend,1,2)*x(1:iend,2,j) &
           &          - a(1:iend,1,3)*x(1:iend,3,j)) / a(1:iend,1,1)
    end do

  end subroutine solve_mat_3


  !---------------------------------------------------------------------
  ! return x = b a^-1 = (a^-t b)^t optimized for 3x3 matrices, where b
  ! is a diagonal matrix, using lu factorization and substitution with
  ! no pivoting.
  pure subroutine diag_mat_right_divide_3(n,iend,a,b,x)
    integer,    intent(in)  :: n, iend
    real(jprb), intent(in)  :: a(iend,3,3)
    real(jprb), intent(in)  :: b(iend,3)
    real(jprb), intent(out) :: x(n,3,3)

    real(jprb), dimension(iend) :: l21, l31, l32
    real(jprb), dimension(iend) :: u22, u23, u33
    real(jprb), dimension(iend) :: y2, y3

    !    associate (u11 => a(:,1,1), u12 => a(:,1,2), u13 => a(1,3))
    ! lu decomposition of the *transpose* of a:
    !       ( 1        )   (u11 u12 u13)
    ! a^t = (l21  1    ) * (    u22 u23)
    !       (l31 l32  1)   (        u33)
    l21 = a(1:iend,1,2) / a(1:iend,1,1)
    l31 = a(1:iend,1,3) / a(1:iend,1,1)
    u22 = a(1:iend,2,2) - l21*a(1:iend,2,1)
    u23 = a(1:iend,3,2) - l21*a(1:iend,3,1)
    l32 =(a(1:iend,2,3) - l31*a(1:iend,2,1)) / u22
    u33 = a(1:iend,3,3) - l31*a(1:iend,3,1) - l32*u23

    ! solve x(1,:) = a^-t ( b(1) )
    !                     (  0   )
    !                     (  0   )
    ! solve ly = b(:,:,j) by forward substitution
    ! y1 = b(:,1)
    y2 = - l21*b(1:iend,1)
    y3 = - l31*b(1:iend,1) - l32*y2
    ! solve ux(:,:,j) = y by back substitution
    x(1:iend,1,3) = y3 / u33
    x(1:iend,1,2) = (y2 - u23*x(1:iend,1,3)) / u22
    x(1:iend,1,1) = (b(1:iend,1) - a(1:iend,2,1)*x(1:iend,1,2) &
         &          - a(1:iend,3,1)*x(1:iend,1,3)) / a(1:iend,1,1)

    ! solve x(2,:) = a^-t (  0   )
    !                     ( b(2) )
    !                     (  0   )
    ! solve ly = b(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = b(1:iend,2)
    y3 = - l32*b(1:iend,2)
    ! solve ux(:,:,j) = y by back substitution
    x(1:iend,2,3) = y3 / u33
    x(1:iend,2,2) = (b(1:iend,2) - u23*x(1:iend,2,3)) / u22
    x(1:iend,2,1) = (-a(1:iend,2,1)*x(1:iend,2,2) &
         &           -a(1:iend,3,1)*x(1:iend,2,3)) / a(1:iend,1,1)

    ! solve x(3,:) = a^-t (  0   )
    !                     (  0   )
    !                     ( b(3) )
    ! solve ly = b(:,:,j) by forward substitution
    ! y1 = 0
    ! y2 = 0
    ! y3 = b(1:iend,3)
    ! solve ux(:,:,j) = y by back substitution
    x(1:iend,3,3) = b(1:iend,3) / u33
    x(1:iend,3,2) = -u23*x(1:iend,3,3) / u22
    x(1:iend,3,1) = (-a(1:iend,2,1)*x(1:iend,3,2) &
         &          - a(1:iend,3,1)*x(1:iend,3,3)) / a(1:iend,1,1)

  end subroutine diag_mat_right_divide_3


  !---------------------------------------------------------------------
  ! treat a as n m-by-m matrices and return the lu factorization of a
  ! compressed into a single matrice (with l below the diagonal and u
  ! on and above the diagonal; the diagonal elements of l are 1). no
  ! pivoting is performed.
  pure subroutine lu_factorization(n, iend, m, a, lu)
    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(out) :: lu(iend,m,m)

    real(jprb) :: s(iend)
    integer    :: j1, j2, j3

    ! this routine is adapted from an in-place one, so we first copy
    ! the input into the output.
    lu(1:iend,1:m,1:m) = a(1:iend,1:m,1:m)

    do j2 = 1, m
      do j1 = 1, j2-1
        s = lu(1:iend,j1,j2)
        do j3 = 1, j1-1
          s = s - lu(1:iend,j1,j3) * lu(1:iend,j3,j2)
        end do
        lu(1:iend,j1,j2) = s
      end do
      do j1 = j2, m
        s = lu(1:iend,j1,j2)
        do j3 = 1, j2-1
          s = s - lu(1:iend,j1,j3) * lu(1:iend,j3,j2)
        end do
        lu(1:iend,j1,j2) = s
      end do
      if (j2 /= m) then
        s = 1.0_jprb / lu(1:iend,j2,j2)
        do j1 = j2+1, m
          lu(1:iend,j1,j2) = lu(1:iend,j1,j2) * s
        end do
      end if
    end do

  end subroutine lu_factorization


  !---------------------------------------------------------------------
  ! treat lu as an lu-factorization of an original matrix a, and
  ! return x where ax=b. lu consists of n m-by-m matrices and b as n
  ! m-element vectors.
  pure subroutine lu_substitution(n,iend,m,lu,b,x)
    ! check: dimensions should be ":"?
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: lu(iend,m,m)
    real(jprb), intent(in) :: b(:,:)
    real(jprb), intent(out):: x(iend,m)

    integer :: j1, j2

    x(1:iend,1:m) = b(1:iend,1:m)

    ! first solve ly=b
    do j2 = 2, m
      do j1 = 1, j2-1
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*lu(1:iend,j2,j1)
      end do
    end do
    ! now solve ux=y
    do j2 = m, 1, -1
      do j1 = j2+1, m
        x(1:iend,j2) = x(1:iend,j2) - x(1:iend,j1)*lu(1:iend,j2,j1)
      end do
      x(1:iend,j2) = x(1:iend,j2) / lu(1:iend,j2,j2)
    end do

  end subroutine lu_substitution


  !---------------------------------------------------------------------
  ! return matrix x where ax=b. lu, a, x, b all consist of n m-by-m
  ! matrices.
  pure subroutine solve_mat_n(n,iend,m,a,b,x)
    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: a(:,:,:)
    real(jprb), intent(in) :: b(:,:,:)
    real(jprb), intent(out):: x(iend,m,m)

    real(jprb) :: lu(iend,m,m)

    integer :: j

    call lu_factorization(n,iend,m,a,lu)

    do j = 1, m
      call lu_substitution(n,iend,m,lu,b(1:,1:m,j),x(1:iend,1:m,j))
!      call lu_substitution(n,iend,m,lu,b(1:n,1:m,j),x(1:iend,1:m,j))
    end do

  end subroutine solve_mat_n


  !---------------------------------------------------------------------
  ! solve ax=b, where a consists of n m-by-m matrices and x and b
  ! consist of n m-element vectors. for m=2 or m=3, this function
  ! calls optimized versions, otherwise it uses general lu
  ! decomposition without pivoting.
  function solve_vec(n,iend,m,a,b)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in) :: n, m, iend
    real(jprb), intent(in) :: a(:,:,:)
    real(jprb), intent(in) :: b(:,:)

    real(jprb)             :: solve_vec(iend,m)
    real(jprb)             :: lu(iend,m,m)
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_vec',0,hook_handle)

    if (m == 2) then
      call solve_vec_2(n,iend,a,b,solve_vec)
    elseif (m == 3) then
      call solve_vec_3(n,iend,a,b,solve_vec)
    else
      call lu_factorization(n,iend,m,a,lu)
      call lu_substitution(n,iend,m,lu,b,solve_vec)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_vec',1,hook_handle)

  end function solve_vec


  !---------------------------------------------------------------------
  ! solve ax=b, where a, x and b consist of n m-by-m matrices. for m=2
  ! or m=3, this function calls optimized versions, otherwise it uses
  ! general lu decomposition without pivoting.
  function solve_mat(n,iend,m,a,b)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)  :: n, m, iend
    real(jprb), intent(in)  :: a(:,:,:)
    real(jprb), intent(in)  :: b(:,:,:)

    real(jprb)              :: solve_mat(iend,m,m)
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:solve_mat',0,hook_handle)

    if (m == 2) then
      call solve_mat_2(n,iend,a,b,solve_mat)
    elseif (m == 3) then
      call solve_mat_3(n,iend,a,b,solve_mat)
    else
      call solve_mat_n(n,iend,m,a,b,solve_mat)
    end if

    if (lhook) call dr_hook('radiation_matrix:solve_mat',1,hook_handle)

  end function solve_mat


  ! --- matrix exponentiation ---

  !---------------------------------------------------------------------
  ! perform matrix exponential of n m-by-m matrices stored in a (where
  ! index n varies fastest) using the higham scaling and squaring
  ! method. the result is placed in a. this routine is intended for
  ! speed so is accurate only to single precision.  for simplicity and
  ! to aid vectorization, the pade approximant of order 7 is used for
  ! all input matrices, perhaps leading to a few too many
  ! multiplications for matrices with a small norm.
  subroutine expm(n,iend,m,a,i_matrix_pattern)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,    intent(in)      :: n, m, iend
    real(jprb), intent(inout)   :: a(n,m,m)
    integer,    intent(in)      :: i_matrix_pattern

    real(jprb), parameter :: theta(3) = (/4.258730016922831e-01_jprb, &
         &                                1.880152677804762e+00_jprb, &
         &                                3.925724783138660e+00_jprb/) 
    real(jprb), parameter :: c(8) = (/17297280.0_jprb, 8648640.0_jprb, &
         &                1995840.0_jprb, 277200.0_jprb, 25200.0_jprb, &
         &                1512.0_jprb, 56.0_jprb, 1.0_jprb/)

    real(jprb), dimension(iend,m,m) :: a2, a4, a6
    real(jprb), dimension(iend,m,m) :: u, v

    real(jprb) :: norma(iend), sum_column(iend)

    integer    :: j1, j2, j3
    real(jprb) :: frac(iend)
    integer    :: expo(iend)
    real(jprb) :: scaling(iend)

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:expm',0,hook_handle)

    norma = 0.0_jprb

    ! compute the 1-norms of a
    do j3 = 1,m
      sum_column(:) = 0.0_jprb
      do j2 = 1,m
        do j1 = 1,iend
          sum_column(j1) = sum_column(j1) + abs(a(j1,j2,j3))
        end do
      end do
      do j1 = 1,iend
        if (sum_column(j1) > norma(j1)) then
          norma(j1) = sum_column(j1)
        end if
      end do
    end do

    frac = fraction(norma/theta(3))
    expo = exponent(norma/theta(3))
    where (frac == 0.5_jprb)
      expo = expo - 1
    end where

    where (expo < 0)
      expo = 0
    end where

    ! scale the input matrices by a power of 2
    scaling = 2.0_jprb**(-expo)
    do j3 = 1,m
      do j2 = 1,m
        a(1:iend,j2,j3) = a(1:iend,j2,j3) * scaling
      end do
    end do
    ! pade approximant of degree 7
    a2 = mat_x_mat(n,iend,m,a, a, i_matrix_pattern)
    a4 = mat_x_mat(n,iend,m,a2,a2,i_matrix_pattern)
    a6 = mat_x_mat(n,iend,m,a2,a4,i_matrix_pattern)

    v = c(8)*a6 + c(6)*a4 + c(4)*a2
    do j3 = 1,m
      v(:,j3,j3) = v(:,j3,j3) + c(2)
    end do
    u = mat_x_mat(n,iend,m,a,v,i_matrix_pattern)
    v = c(7)*a6 + c(5)*a4 + c(3)*a2
    ! add a multiple of the identity matrix
    do j3 = 1,m
      v(:,j3,j3) = v(:,j3,j3) + c(1)
    end do

    v = v-u
    u = 2.0_jprb*u
    a(1:iend,1:m,1:m) = solve_mat(n,iend,m,v,u)

    ! add the identity matrix
    do j3 = 1,m
      a(1:iend,j3,j3) = a(1:iend,j3,j3) + 1.0_jprb
    end do

    ! loop through the matrices
    do j1 = 1,iend
      if (expo(j1) > 0) then
        ! square matrix j1 expo(j1) times          
        a(j1,:,:) = repeated_square(m,a(j1,:,:),expo(j1),i_matrix_pattern)
      end if
    end do

    if (lhook) call dr_hook('radiation_matrix:expm',1,hook_handle)

  end subroutine expm


  !---------------------------------------------------------------------
  ! return the matrix exponential of n 2x2 matrices representing
  ! conservative exchange between spartacus regions, where the
  ! matrices have the structure
  !   (-a   b)
  !   ( a  -b)
  ! and a and b are assumed to be positive or zero.  the solution uses
  ! putzer's algorithm - see the appendix of hogan et al. (gmd 2018)
  subroutine fast_expm_exchange_2(n,iend,a,b,r)

    use ecradhook, only : lhook, dr_hook, jphook

    integer,                      intent(in)  :: n, iend
    real(jprb), dimension(n),     intent(in)  :: a, b
    real(jprb), dimension(n,2,2), intent(out) :: r

    real(jprb), dimension(iend) :: factor

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',0,hook_handle)

    ! security to ensure that if a==b==0 then the identity matrix is returned
    factor = (1.0_jprb - exp(-(a(1:iend)+b(1:iend))))/max(1.0e-12_jprb,a(1:iend)+b(1:iend))

    r(1:iend,1,1) = 1.0_jprb - factor*a(1:iend)
    r(1:iend,2,1) = factor*a(1:iend)
    r(1:iend,1,2) = factor*b(1:iend)
    r(1:iend,2,2) = 1.0_jprb - factor*b(1:iend)

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_2',1,hook_handle)

  end subroutine fast_expm_exchange_2


  !---------------------------------------------------------------------
  ! return the matrix exponential of n 3x3 matrices representing
  ! conservative exchange between spartacus regions, where the
  ! matrices have the structure
  !   (-a   b   0)
  !   ( a -b-c  d)
  !   ( 0   c  -d)
  ! and a-d are assumed to be positive or zero.  the solution uses the
  ! diagonalization method and is a slight generalization of the
  ! solution provided in the appendix of hogan et al. (gmd 2018),
  ! which assumed c==d.
  subroutine fast_expm_exchange_3(n,iend,a,b,c,d,r)

    use ecradhook, only : lhook, dr_hook, jphook

    real(jprb), parameter :: my_epsilon = 1.0e-12_jprb

    integer,                      intent(in)  :: n, iend
    real(jprb), dimension(n),     intent(in)  :: a, b, c, d
    real(jprb), dimension(n,3,3), intent(out) :: r

    ! eigenvectors
    real(jprb), dimension(iend,3,3) :: v

    ! non-zero eigenvalues
    real(jprb), dimension(iend) :: lambda1, lambda2

    ! diagonal matrix of the exponential of the eigenvalues
    real(jprb), dimension(iend,3) :: diag

    ! result of diag right-divided by v
    real(jprb), dimension(iend,3,3) :: diag_rdivide_v

    ! intermediate arrays
    real(jprb), dimension(iend) :: tmp1, tmp2

    integer :: j1, j2

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',0,hook_handle)

    ! eigenvalues lambda1 and lambda2
    tmp1 = 0.5_jprb * (a(1:iend)+b(1:iend)+c(1:iend)+d(1:iend))
    tmp2 = sqrt(max(0.0_jprb, tmp1*tmp1 - (a(1:iend)*c(1:iend) &
         &                    + a(1:iend)*d(1:iend) + b(1:iend)*d(1:iend))))
    ! the eigenvalues must not be the same or the lu decomposition
    ! fails; this can occur occasionally in single precision, which we
    ! avoid by limiting the minimum value of tmp2
    tmp2 = max(tmp2, epsilon(1.0_jprb) * tmp1)
    lambda1 = -tmp1 + tmp2
    lambda2 = -tmp1 - tmp2

    ! eigenvectors, with securities such that if a--d are all zero
    ! then v is non-singular and the identity matrix is returned in r;
    ! note that lambdax is typically negative so we need a
    ! sign-preserving security
    v(1:iend,1,1) = max(my_epsilon, b(1:iend)) &
         &  / sign(max(my_epsilon, abs(a(1:iend) + lambda1)), a(1:iend) + lambda1)
    v(1:iend,1,2) = b(1:iend) &
         &  / sign(max(my_epsilon, abs(a(1:iend) + lambda2)), a(1:iend) + lambda2)
    v(1:iend,1,3) = b(1:iend) / max(my_epsilon, a(1:iend))
    v(1:iend,2,:) = 1.0_jprb
    v(1:iend,3,1) = c(1:iend) &
         &  / sign(max(my_epsilon, abs(d(1:iend) + lambda1)), d(1:iend) + lambda1)
    v(1:iend,3,2) = c(1:iend) &
         &  / sign(max(my_epsilon, abs(d(1:iend) + lambda2)), d(1:iend) + lambda2)
    v(1:iend,3,3) = max(my_epsilon, c(1:iend)) / max(my_epsilon, d(1:iend))
    
    diag(:,1) = exp(lambda1)
    diag(:,2) = exp(lambda2)
    diag(:,3) = 1.0_jprb

    ! compute diag_rdivide_v = diag * v^-1
    call diag_mat_right_divide_3(iend,iend,v,diag,diag_rdivide_v)

    ! compute v * diag_rdivide_v
    do j1 = 1,3
      do j2 = 1,3
        r(1:iend,j2,j1) = v(1:iend,j2,1)*diag_rdivide_v(1:iend,1,j1) &
             &          + v(1:iend,j2,2)*diag_rdivide_v(1:iend,2,j1) &
             &          + v(1:iend,j2,3)*diag_rdivide_v(1:iend,3,j1)
      end do
    end do

    if (lhook) call dr_hook('radiation_matrix:fast_expm_exchange_3',1,hook_handle)

  end subroutine fast_expm_exchange_3

!  generic :: fast_expm_exchange => fast_expm_exchange_2, fast_expm_exchange_3

 
end module radiation_matrix
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

