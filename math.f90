module math
  
  integer, parameter :: r4k = KIND(1.0)                  !single precision
  integer, parameter :: r8k = SELECTED_REAL_KIND(2*PRECISION(1.0_r4k))
  integer, parameter :: rk=r8k !this is the real kinddd we are using 

  INTEGER, PARAMETER :: isoluku  = SELECTED_INT_KIND(12)
  
  !pi = 4.0_rk*atan(1.0_rk)
  real(kind = rk), parameter :: pi = 3.14159265358979_rk 
  real(kind = rk), parameter :: kB = 3.16681050847798e-06_rk
  real(kind = rk), parameter :: meV2Hartree = 1.0_rk/27.211385_rk/1000
  real(kind = rk), parameter :: nm2bohr = 10.0_rk/0.529177_rk
  INTEGER, parameter              :: dimensio = 2
  real(kind = rk), save      :: dielectric_constant
  real(kind = rk), save      :: omega2

  real(kind = rk), save :: lambdaMS,ZZMS

contains
  
  FUNCTION gaussinen() RESULT( Gauss_taul )
    !--------------------------------------------------------
    ! Tama generoi normaalisti jakautuneita satunnaislukuja
    !
    
    !Box-Muller algoritmi
    
    
    REAL(kind = rk), DIMENSION(dimensio)      :: Gauss_taul
    REAL(kind = rk), DIMENSION(2)      :: ksi
    
    CALL RANDOM_NUMBER( ksi )
    
    gauss_taul(1) = SQRT( -2.0_rk*LOG(ksi(1)) )*COS( 2.0_rk*pi*ksi(2) )
    gauss_taul(2) = SQRT( -2.0_rk*LOG(ksi(1)) )*SIN( 2.0_rk*pi*ksi(2) )
    
    if (dimensio>2) then
       !Toinen kerta toden sanoo
       CALL RANDOM_NUMBER( ksi )
       
       gauss_taul(dimensio) = SQRT( -2.0_rk*LOG(ksi(1)) )*COS( 2.0_rk*pi*ksi(2) )
    end if
    
  END FUNCTION gaussinen

  function linspace(x,y,dim) result(r)
    implicit none
    integer :: dim, i
    real(kind = rk) :: x,y, dr
    real(kind = rk), dimension(dim) :: r

    
    dr = (y-x)/(dim-1)
    
    do i = 1,dim
       r(i) = x + (i-1)*dr
    end do
    
  end function linspace

  
  FUNCTION erfc(x)
    !Numerical Recipes:
    !Returns the complementary error function erfc(x) with
    !fractional error everywhere less than 1.2 function erfcc(x)
    !Numerical Recipes:
    IMPLICIT NONE
    REAL(kind=rk) :: x,z,t,ans,erfc
    z=ABS(x);
    t=1.0_rk/(1.0_rk+z/2);
    ans=t*EXP(-z*z-1.26551223_rk+t*(1.00002368_rk+t*(0.37409196_rk+t*(0.09678418_rk+&
        t*(-0.18628806_rk+t*(0.27886807_rk+t*(-1.13520398_rk+t*(1.48851587_rk+&
        t*(-0.82215223_rk+t*0.17087277_rk)))))))));
    IF (x>=0.0_rk) THEN
       erfc = ans
    ELSE
       erfc = 2.0_rk - ans
    END IF
  END FUNCTION erfc



  FUNCTION erff(t)
    IMPLICIT NONE
    REAL(kind=rk) :: t,erff

    erff = 1.0_rk - erfc(t)

  END FUNCTION erff
  
  subroutine LUdecomp(A,dim)
    implicit none
    integer, intent(in) :: dim
    real(kind = rk), dimension(dim,dim),intent(inout) :: A
    real(kind = rk), dimension(dim,dim) :: L,U
    integer :: i,j,k

    L = 0.0_rk
    U = 0.0_rk

    do i=1,dim
  
       ! non-zero diagonal assumed
  
       L(i,i) = 1.0_rk
  
       do j=i,dim
          U(i,j) = A(i,j)
          do k=1,i-1
             U(i,j) = U(i,j) - L(i,k)*U(k,j)
          end do
       end do
  
       do j=i+1,dim
          L(j,i) = A(j,i)
          do k=1,i-1
             L(j,i) = L(j,i) - L(j,k)*U(k,i)
          end do
          if (abs(U(i,i))<1.0e-100_rk) then
             L(j,i) = L(j,i)*sign(1.0e100_rk,U(i,i))
          else
             L(j,i) = L(j,i)/U(i,i)
          end if
       end do
    end do
    
    do i=1,dim
       do j=i,dim
          A(j,i) = L(j,i)
          A(i,j) = U(i,j)
       end do
    end do

  end subroutine LUdecomp


  
  FUNCTION det(A,M) RESULT( detvalue )
    IMPLICIT NONE
    !EXTERNAL :: SGETRF
    REAL(KIND=RK), DIMENSION(M,M), INTENT(out) :: A
    REAL(KIND=RK)                             :: detvalue
    INTEGER, INTENT(in)              :: M
    !INTEGER                          :: INFO,i,j
    integer :: i,j
    !INTEGER, DIMENSION(M)            :: IPIV

    !DO i = 1,M
    !  IPIV(i) = i
    !END DO
    !INFO = 0

    ! LU decomposition
    !CALL SGETRF(M,M,A,M,IPIV,INFO)
    
    call Pivoting(A,M,j)
    call LUdecomp(A,M)
    
    ! value of the determinant
    !j = 0
    detvalue = 1.0_rk
    DO i = 1,M
       !IF (IPIV(i) .NE. i) THEN
       !  j = j + 1
       !END IF
       
       detvalue = detvalue*A(i,i)
       
    END DO
    
    detvalue = ((-1.0)**j)*detvalue
    
  END FUNCTION det

  subroutine Pivoting(A,M,j)
    implicit none
    integer, intent(in) :: M
    real(kind=rk), dimension(M,M), intent(inout) :: A
    real(kind=rk), dimension(1,M) :: row
    real(kind=rk), dimension(M,1) :: column
    integer, intent(out) :: j
    integer :: i,k,mm,nn,ind
    real(kind=rk) :: largest
    j=0
    mm=1
    nn=1
    
    
    do ind=1,M-1
       largest=0.0_rk 
       do i=ind,M
          do k=ind,M
             if (abs(A(i,k))>largest) then
                largest=abs(A(i,k))
                mm=i
                nn=k
             end if
          end do
       end do
       
       if (mm .ne. ind) then
          row(1,:)=A(mm,:)
          A(mm,:)=A(ind,:)
          A(ind,:)=row(1,:)
          j=j+1
       end if
       
       if (nn .ne. ind) then
          column(:,1)=A(:,nn)
          A(:,nn)=A(:,ind)
          A(:,ind)=column(:,1)
          j=j+1
       end if
    end do
    
  end subroutine Pivoting


  function index_check(ind,limit_down,limit_up) result(inbounds)
    implicit none
    integer :: ind, limit_down, limit_up
    logical :: inbounds
    
    if (ind >=1 .and. ind >= limit_down .and. ind <= limit_up) then
       inbounds = .true.
    else
       inbounds = .false.
    end if
    
  end function index_check
  
  
end module math
