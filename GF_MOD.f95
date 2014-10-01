!compile with 
! f2py --fcompiler=gfortran -llapack --opt=-O -c -m FMod F_GF_MOD.f95 
! install lapack64
! install blast64

!Just the single integral (other performed by contour integration)
module shared_data
implicit none
  save
  ! Mathematical parameters
  complex(8), parameter :: im = (0.0d0,1.0d0)
  real(8), parameter :: pi = acos(-1.0d0)
  ! Material parameters
  real(8), parameter :: t = -1.0d0
  
  
end module


! Fortran inverions Routine
subroutine inv(A, Ainv, n)
implicit none
    ! Input arguments
    integer, intent(in) :: n
    complex(8), dimension(n,n), intent(in) :: A
    complex(8), dimension(n,n), intent(out) :: Ainv
    ! Dummy arguments
    complex(8), dimension(n) :: work  ! work array for LAPACK
    integer, dimension(n) :: ipiv   ! pivot indices
    integer :: info
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    ! ZGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call ZGETRF(n, n, Ainv, n, ipiv, info)

    ! This is overkill, you either need to ignore this, write it to a file, or write it to the error output.
    if (info /= 0) then
      stop 'Matrix is numerically singular!'
    end if

    ! ZGETRI computes the inverse of a matrix using the LU factorization
    ! computed by ZGETRF.
    call ZGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
      stop 'Matrix inversion failed!'
    end if

end subroutine inv



! Bulk Graphene Green's Function
function gf_bulk_kz(m,n,s_lat,E,kZ)
use shared_data
implicit none
  complex(8) :: GF_bulk_kz
  
  complex(8), intent(in) :: E
  real(8), intent(in) :: kZ
  integer, intent(in) :: m,n
  integer, intent(in) :: s_lat

  complex(8) :: f 
  complex(8) :: q
  complex(8) :: Const, Den
  integer :: sig

  q = acos( (E**2 - t**2 - 4.0d0*t**2 *cos(kZ)**2)/(4.0d0*t**2 *cos(kZ) ) )

  if (aimag(q) < 0.0d0) q = -q

  Const = im/(4*pi*t**2)
  Den = cos(kZ)*sin(q)

  if (s_lat == 0) then
    sig = sign(1,m+n)
    gf_bulk_kz = Const*E*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/ Den 
    
  else if (s_lat > 0) then 
    sig = sign(1,m+n)
    f = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(sig*im*q) )
    gf_bulk_kz = Const*f*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den  
    
  else 
    sig = sign(1,m+n-1)
    f = t*( 1.0d0 + 2.0d0*cos(kZ)*exp(-sig*im*q) )
    gf_bulk_kz = Const*f*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den 
    
  end if
        
end function gf_bulk_kz 


! Bulk Graphene Green's Function with strain
function gf_bulk_kz_strain( t1, t2, m, n , s, E, kZ )
use shared_data
implicit none
  complex(8) :: gf_bulk_kz_strain
  
  complex(8), intent(in) :: E
  real(8), intent(in) :: kZ
  real(8), intent(in) :: t1, t2
  integer, intent(in) :: m,n
  integer, intent(in) :: s

  complex(8) :: f 
  complex(8) :: q
  complex(8) :: Const, Den
  integer :: sig

  q = acos( (E**2 - (t1*t1) - 4.0d0*t2*t2*cos(kZ)**2)/(4.0d0*t1*t2*cos(kZ) ) )

  if (aimag(q) < 0.0d0) q = -q

  Const = im/(4*pi*t1*t2)
  Den = cos(kZ)*sin(q)

  if (s == 0) then
    sig = sign(1,m+n)
    gf_bulk_kz_strain = Const*E*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/ Den 
    
  else if (s > 0) then 
    sig = sign(1,m+n)
    f = ( t1 + 2.0d0*t2*cos(kZ)*exp(sig*im*q) )
    gf_bulk_kz_strain = Const*f*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den  
    
  else 
    sig = sign(1,m+n-1)
    f = ( t1 + 2.0d0*t2*cos(kZ)*exp(-sig*im*q) )
    gf_bulk_kz_strain = Const*f*exp( im*(sig*q*(m+n) + kZ*(m-n) ) )/Den 
    
    
  end if
end function gf_bulk_kz_strain 
        

    
