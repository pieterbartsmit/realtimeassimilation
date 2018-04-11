module modTools

  use modPar
  use IO
  implicit none
  
contains

    subroutine interpDep(Xi,Yi,Zi,np)
! 
!*****************************************************************************
!
!

!-----------------------------------------------------------------------------
!                                   INTERFACE
!-----------------------------------------------------------------------------

!   --------- DESCRIPTION  --------

!   Function which finds the bi-linear interpolation Zi at Xi,Yi from the function Z defined at X,Y.
!   Note that both grids (X,Y and Xi,Yi) must be defined in the same coordinate system. For values outside the input
!   grid the closest defined values are taken (Nearest interpolation).
!   

!   --------- DEPENDENCIES --------
    use modPar 
    
    implicit none

!   --------- VARIABLES --------
!
!   ARGUMENT     
    real(kind=rKind),intent(in) ,dimension(np) :: Xi       ! x - coordinate to interpolate Z
    real(kind=rKind),intent(in) ,dimension(np) :: Yi       ! y - coordinate to interpolate Z
    real(kind=rKind),intent(out),dimension(np) :: Zi       ! Interpolated values at (Xi,Yi)
    
    integer(kind=iKind),intent(in)                  :: np       ! Number of points in x-dir
        
!   LOCAL
    integer(kind=iKind)                             :: ip      ! loop index
    integer(kind=iKind),dimension(2)                :: ix       ! index of points surrounding Xi(ixi)
    integer(kind=iKind),dimension(2)                :: iix       ! index of points surrounding Xi(ixi)    
    integer(kind=iKind),dimension(2)                :: iy       ! index of points surrounding Yi(iyi)
    integer(kind=iKind)                             :: is       ! index of points surrounding Yi(iyi)    
    real(kind=rKind)                                :: dx       ! meshsize input grid
    real(kind=rKind)                                :: dy       ! meshsize input grid
    real(kind=rKind)                                :: dxb      ! relative distance
    real(kind=rKind)                                :: dyb      ! relative distance
    real(kind=rKind)                                :: dif      ! relative distance    
    
!    
!-----------------------------------------------------------------------------
!                                   SOURCE
!-----------------------------------------------------------------------------
!
! Interpolation by searching, appears to give better results than gridinterp, but why? Roundoff?
!
    !

    dx = (X(mx) - X(1))/(mx-1)
    dy = (Y(my) - Y(1))/(my-1)
    !
    do ip = 1, np
       !
       ! Find interpolation values
       !

       if ( Yi(ip) >= Y(my) ) then
          !
          iy(1) = my
          iy(2) = my
          !
       elseif (Yi(ip) <= Y(1) ) then
          !
          iy(1) = 1
          iy(2) = 1
          !
       else
          !
          dif = ( Yi(ip) - Y(1) )/dy
          iy(1) = floor( dif  ) + 1
          iy(2) =  iy(1)+1             
         ! searchloop: do is = 1,my-1
             !
             !if ( Yi(ip) >= Y(is) .and. Yi(ip) < Y(is+1)  ) then
                !
              !  iy(1) = is
             !   iy(2) = is+1
            !    exit searchloop
                !
           !  endif
             !
          !enddo searchloop
          !
          !call binsearch( Y , my,yi(ip) , iy )
          !
       end if
       
       if ( Xi(ip) > X(mx) ) then
          !
          ix(1) = mx
          ix(2) = mx
          !
       elseif (Xi(ip) < X(1) ) then
          !
          ix(1) = 1
          ix(2) = 1
          !
       else
          !
          dif = ( Xi(ip) - X(1) )/dx 
          ix(1) = floor( dif ) + 1
          ix(2) = ix(1)+1          
          !
       end if
 
       !

       if (  iy(2) > iy(1) ) then
          !
          dyb = ( Yi(ip)-Y(iy(1)) ) /dy
          !
       else
          !
          dyb = 0.
          !
       endif

       if ( ix(2) > ix(1) ) then
          !
          dxb = ( Xi(ip)-X(ix(1)) ) /dx
          !
       else
          !
          dxb = 0.
          !
       endif

       !
       !bi-linear interpolation from surrounding values
       
       Zi(ip) = dep(ix(1),iy(1))*(1.-dyb)*(1.-dxb) + &
                dep(ix(2),iy(1))*(1.-dyb)*(   dxb) + &
                dep(ix(2),iy(2))*(   dyb)*(   dxb) + &
                dep(ix(1),iy(2))*(   dyb)*(1.-dxb)   

        !
    enddo
    !
  end subroutine interpDep

  subroutine binsearch( x , mx,xi , lim )
    !
    !   --------- DEPENDENCIES --------
    use modPar, only: rkind,ikind
    
    integer(kind=iKind),intent(in)    :: mx
    integer(kind=iKind),intent(inout) :: lim(2)
    real(kind=rKind),intent(in)       :: x(mx)
    real(kind=rKind),intent(in)       :: xi
    integer(kind=iKind)               :: inbetween
    integer(kind=iKind)               :: ifailsafe    

    lim(1) = 1
    lim(2) = mx
    ifailsafe = 1
    binsrch: do
       !
       ifailsafe = ifailsafe+1
       inbetween = ( lim(2) + lim(1) ) / 2
       if ( xi >= x(lim(1)) .and. xi < x( inbetween ) ) then
          !
          lim(2) = inbetween
          !
       else
          !
          lim(1) = inbetween
          !
       endif
       !

       if ( lim(2) - lim(1) <= 1 ) then
          !
          exit binsrch
          !
       elseif ( ifailsafe > mx ) then
          !
          exit binsrch
          !
       endif
       !
    enddo binsrch

    lim(2) = lim(1) + 1
  end subroutine binsearch
  
  



  function inBndlist( iBnd )
    !
    ! check if Ibnd is in the bndlist
    !
    logical :: inBndlist
    integer(kind=iKind), intent(in) :: iBnd
    integer(kind=iKind) :: jBnd    
    
    inBndlist = .false.
    if (allocated( bndlist )) then
       !
       bndloop: do jBnd = 1 , nbnd
          !
          if (bndlist( jBnd ) == iBnd ) then
             inBndlist = .true.
             exit bndloop
          endif
          !
       enddo bndloop
       !
    endif
    
    
  end function inBndlist

  
  
  function inDomain( x,y)
    !
    integer(kind=ikind) :: inDomain

    real(kind=rKind) :: x
    real(kind=rKind) :: y


    !
    indomain = 0
    if ( x <= xlim(2) .and. x >= xlim(1) .and. y <= ylim(2) .and. y >= ylim(1) ) then
       !
       indomain = -1
    else
       !
       if ( x < xlim(1) ) then
          ! W
          indomain = 1
          !
       elseif ( x > xlim(2) ) then
          ! E
          indomain = 2
          !
       elseif ( y > ylim(2) ) then
          ! N
          indomain = 4
          !
       elseif ( y < ylim(1) ) then
          ! S
          indomain = 3
          !
       endif
       !
    endif
    !
  end function inDomain

  subroutine lineCrossing( linescross, xo,yo, xa1,xa2,ya1,ya2,xb1,xb2,yb1,yb2)
    !
    ! Do two line segments cross?
    !
    logical         ,intent(out) :: linescross
    real(kind=rKind),intent(out) :: xo
    real(kind=rKind),intent(out) :: yo    
    
    real(kind=rKind),intent(in) :: xa1
    real(kind=rKind),intent(in) :: xa2
    real(kind=rKind),intent(in) :: ya1
    real(kind=rKind),intent(in) :: ya2
    real(kind=rKind),intent(in) :: xb1
    real(kind=rKind),intent(in) :: xb2
    real(kind=rKind),intent(in) :: yb1
    real(kind=rKind),intent(in) :: yb2    

    real(kind=rKind)            :: mat(2,2)
    real(kind=rKind)            :: inv(2,2)
    real(kind=rKind)            :: Det    
    real(kind=rKind)            :: rhs(2)
    real(kind=rKind)            :: vec(2)       

    real(kind=rKind)            :: dxa,dxb,dya,dyb

    ! ex1 = xa1 + alpha * dxa = xb1 + beta * dxb 
    ! ex2 = ya1 + alpha * dya = yb1 + beta * dyb
    !
    ! or
    !
    ! 
    
    dxa = xa2-xa1
    dya = ya2-ya1
    dxb = xb2-xb1
    dyb = yb2-yb1

    mat(1,1) =   dxa
    mat(1,2) = - dxb
    mat(2,1) =   dya
    mat(2,2) = - dyb

    !call log_real('mat(1,1)', mat(1,1),1)
    !call log_real('mat(1,2)', mat(1,2),1)
    !call log_real('mat(2,1)', mat(2,1),1)
    !call log_real('mat(2,2)', mat(2,2),1)
    
    rhs(1)   = xb1 - xa1
    rhs(2)   = yb1 - ya1

    Det      = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
        
    if ( abs(Det) < 1.e-14 ) then
       !
       linescross = .false.
       xo = 0.
       yo = 0.
       return
       !
    endif

    call matinv(mat,inv)
    vec = matmul( inv , RHS )

    !
    linescross = .false.
    !
    xo = xa1 + dxa * vec(1)
    yo = ya1 + dya * vec(1)

    if ( (vec(1) >= 0.) .and. (vec(1) <= 1.) .and. (vec(2) >= 0.) .and. (vec(2) <= 1.) ) then
       !
       ! The lines are crossing; what is the mutual angle?
       !
       linescross = .true.                   
       !
    endif
    
   ! call log_real('xo=',xo,1)
   ! call log_real('yo=',yo,1)    
    !
  end subroutine lineCrossing

  subroutine matinv( A , B)
    !
    ! 2 by 2 matrix inversion of matrix A
    !
    real(kind=rKind), intent(in )  :: A(2,2)
    real(kind=rKind), intent(out)  :: B(2,2)
    real(kind=rKind)               :: D
    !
    D = 1./ ( A(1,1)*A(2,2) - A(2,1)*A(1,2) )
    B(1,1) = A(2,2)
    B(2,2) = A(1,1)
    B(1,2) = -A(1,2)
    B(2,1) = -A(2,1)
    B      = B*D
    !
  end subroutine matinv

  elemental function disperK( omega , dep )
    !
    ! Checked this algorithm on tue april 10, 2018
    !
    real(kind=rKind) :: disperK
    real(kind=rKind),intent(in) :: omega
    real(kind=rKind),intent(in) :: dep

    real(kind=rKind) :: k_dw
    real(kind=rKind) :: kd
    real(kind=rKind) :: k_guess
    real(kind=rKind) :: F
    real(kind=rKind) :: dF
    real(kind=rKind) :: th

    integer(kind=iKind) :: iter
    
    k_dw = omega**2 / g
    
    kd   = k_dw * dep
    
    if ( kd > 1 ) then
       !
       k_guess = k_dw
       !
    else
       !
       k_guess =  omega/ sqrt( g * dep )
       !
    endif
    !

    !
    disperK = k_guess
    ITERLOOP: do iter = 1,3
       !
       kd =  disperK * dep
       th = tanh( kd )
       F  = k_dw - disperK*th

       if ( abs(F) < k_dw * 1.d-6 ) then
          !
          exit ITERLOOP
          !
       endif
       
       dF = - th - kd + kd * th**2

       
       disperK = disperK - F / dF
       !
    enddo ITERLOOP
    !  
    !
  end function disperK

  elemental real(kind=rKind) function an_dcdh(k,dep,c)
    !
    
    implicit none
    !
    real(kind=rKind),intent(in)     :: k
    real(kind=rKind),intent(in)     :: c          
    real(kind=rKind),intent(in)     :: dep
    real(kind=rKind)                :: kd
    !

    kd = k*dep
    an_dcdh = .5 * g / c * ( 1 - tanh(kd)**2 )
    return
    if (kd<0.001)     then
       !
       an_dcdh = 0.5_rkind * sqrt(g)/sqrt(dep)
       !
    elseif (kd>10) then
       !
       an_dcdh = 0.0_rKind
       !
    else
       !            
       an_dcdh = .5 * g / c * ( 1 - tanh(kd)**2 )
       !
    endif
         !
  end function an_dcdh

  elemental real(kind=rKind) function an_dkdh(k,dep)
    !

    implicit none
    !
    real(kind=rKind),intent(in)     :: k          
    real(kind=rKind),intent(in)     :: dep
    real(kind=rKind)                :: kd
    real(kind=rKind)                :: ch2
    real(kind=rKind)                :: th    
    !
    kd = k*dep
    ch2 = cosh(kd)**2
    th = tanh(kd)
    !
    if (kd<0.001)     then
       !
       an_dkdh = - k / 2. / dep
       !
    elseif (kd>10) then
       !
       an_dkdh = 0.0_rKind
       !
    else
       !            
       an_dkdh = - g * k**2 / ch2 / ( g * th + g * k * dep / ch2 )
       !
    endif
    !
  end function an_dkdh


  subroutine dDepdS( x, y, dDdx , dDdy, dDdS,dx, dy )
    !
    real( kind = rKind ), intent(in) :: x
    real( kind = rKind ), intent(in) :: y
    real( kind = rKind ), intent(in) :: dx
    real( kind = rKind ), intent(in) :: dy
    real( kind = rKind ), intent(out) :: dDdx
    real( kind = rKind ), intent(out) :: dDdy          
    real( kind = rKind ), intent(out) :: dDds    

    real( kind = rKind )             :: xp(4)
    real( kind = rKind )             :: yp(4)
    real( kind = rKind )             :: di(4)



    xp = [ x - 0.5*dx , x + 0.5 * dx,  x          , x            ]
    yp = [ y          , y           ,  y - 0.5*dy , y + 0.5 * dy ]

    call  interpDep(xp,yp,di,4)

    dDdx = ( di(2) - di(1) )/ dx
    dDdy = ( di(4) - di(3) )/ dy

    dDdS = sqrt( dDdx**2 + dDdy**2 )

  end subroutine dDepdS
        
  real(kind=rKind) function parderK(x,y,freq,dx,dy)
    !
    real( kind = rKind ), intent(in) :: x
    real( kind = rKind ), intent(in) :: y
    real( kind = rKind ), intent(in) :: freq
    real( kind = rKind ), intent(in) :: dx
    real( kind = rKind ), intent(in) :: dy

    real( kind = rKind )             :: xp(3)
    real( kind = rKind )             :: yp(3)
    real( kind = rKind )             :: ds,di(3),omega(3)
    real( kind = rKind )             :: k(3)

    omega = freq
    xp = [ x - 0.5*dx , x, x + 0.5 * dx ]
    yp = [ y - 0.5*dy , y, y + 0.5 * dy ]

    ds = sqrt(dx**2 + dy**2)
    call  interpDep(xp,yp,di,3)
    k = disperK( omega , di )


    parderK = ( k(3) - k(1) ) / ds / k(2)

  end function parderK

  elemental real(kind=rKind) function bisectAngle( angle1 , angle2)
    !
    ! This function returns the angle inbetween two given angles
    !
    real(kind = rKind ), intent(in) :: angle1
    real(kind = rKind ), intent(in) :: angle2

    real( kind = rKind )            :: vec1(2)
    real( kind = rKind )            :: vec2(2)
    real( kind = rKind )            :: vec(2)

    vec1 = [cos(angle1),sin(angle1) ]
    vec2 = [cos(angle2),sin(angle2) ]
    vec  = vec1 / 2 + vec2 / 2

    bisectAngle = atan2( vec(2) , vec(1) )

  end function bisectAngle
end module modTools
