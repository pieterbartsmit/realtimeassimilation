module modPar
  !
  integer(4)          ,parameter  :: rkind = 8          ! Default real
  integer(4)          ,parameter  :: ikind = 4         ! Default integer
  real(rkind)         ,parameter  :: g = 9.81          ! Gravitational constant
  complex(kind=rkind) ,parameter  :: i1 = cmplx(0.,1.)  ! Imaginary unit
  real(kind=rkind)    ,parameter  :: PI = acos(-1._4)   ! Pi
  real(kind=rkind)    ,parameter  :: PI2 = 2.*Pi        ! 2*Pi
  real(kind=rkind)    ,parameter  :: kd_deep = 10.       ! Deep water limit
  real(kind=rKind)                :: drylim = 0.1

  character(len=100) :: msg

  ! Numerical Parameters
  integer(kind=iKind)             :: num_stepSizeControlMethod = 1
  integer(kind=iKind)             :: num_minLevel = 4
  integer(kind=iKind)             :: num_maxStep = 100000 !Number   of steps
  real(kind=rKind)                :: num_frac    = 0.025   !Fraction of wavelength
  real(kind=rKind)                :: num_frac_k  = 0.001   !Fraction of wavelength
  real(kind=rKind)                :: num_dist_rel   = .05
  real(kind=rKind)                :: num_dist_ratio = 1.05
  
  real(kind=rKind)                :: num_stepsize_maxdegPerStep  = 0.5    !Fraction of wavelengthnum_stepsize_fracBatGrid
  real(kind=rKind)                :: num_stepsize_fracBatGrid = 1
  real(kind=rKind)                :: num_stepsize_fracK = 0.001
  real(kind=rKind)                :: num_stepsize_accuracy = 0.001

    integer, parameter :: inumpar   = 23
    integer, parameter :: iDistf   = 23 
    integer, parameter :: idepf   = 22  
    integer, parameter :: irotf   = 21
    integer, parameter :: iyf     = 20
    integer, parameter :: ixf     = 19
    integer, parameter :: idnormal= 18
    integer, parameter :: inormal = 17
    integer, parameter :: iRot    = 16
    integer, parameter :: iFric   = 15
    integer, parameter :: iDist   = 14
    integer, parameter :: iTime   = 13
    integer, parameter :: ic0     = 12
    integer, parameter :: icg0    = 11
    integer, parameter :: id0     = 10
    integer, parameter :: ith0    =  9
    integer, parameter :: iy0     =  8
    integer, parameter :: ix0     =  7
    integer, parameter :: ic1     =  6
    integer, parameter :: icg1    =  5
    integer, parameter :: id1     =  4
    integer, parameter :: ith1    =  3
    integer, parameter :: iy1     =  2
    integer, parameter :: ix1     =  1  
  

  real(kind=rKind),allocatable    :: xlim(:)
  real(kind=rKind),allocatable    :: ylim(:)
  !

  real(kind=rKind),allocatable    :: Dep(:,:)
  real(kind=rKind),allocatable    :: x(:)
  real(kind=rKind),allocatable    :: y(:)
  integer(kind=rKind),allocatable :: bndlist(:)
  integer(kind=rKind),allocatable :: nbnd
  
  real(kind=rKind)                :: dx_bat
  real(kind=rKind)                :: dy_bat  

 ! real(kind=rKind),allocatable    :: cont_xy(:,:)
 ! integer(kind=iKind)             :: cont_n
  real(kind=rKind)                :: cont_d
  logical                         :: cont_set = .false.

  integer(kind=iKind)             :: mx
  integer(kind=iKind)             :: my

 ! integer(kind=iKind)             :: stepsizelimit( 4 ) = 0

  logical          :: spherical = .true.
  real(kind=rKind) :: Rearth    = 6371008.8

!   2.E7/PI 
  real(kind=rKind) :: Lat0
  real(kind=rKind) :: Lon0    
  
  
end module modPar
