module modRayIntegration
  !
  use modPar
  use IO
  !
  implicit none

  !
contains
  !
  subroutine rk4step( yin,yout,freq,ds, dry  )
    !
    ! nth-order Runge Kutta method
    !
    !   
    ! OUTPUT
    !-------------------------------------------------------------------------------------------------------------------
    ! 
    !
    real(kind=rKind), intent(out) :: yout(3)  !solution vector: [x,y,theta]

    !   
    ! INPUT
    !-------------------------------------------------------------------------------------------------------------------
    !     
    real(kind=rKind)   , intent(in) :: yin(3)  !x-velocity component
    real(kind=rKind)   , intent(in) :: freq
    real(kind=rKind)   , intent(in) :: ds !stepsize
    logical            , intent(out) :: dry
    !   
    ! LOCAL
    !-------------------------------------------------------------------------------------------------------------------
    ! 
    !

    !
    ! Private variables
    !
    integer(kind=iKind)       :: ik
    integer(kind=iKind)       :: nRkStage !Number of RKStages
    real(kind = rKind)        :: stagecoef(4)
    real(kind = rKind)        :: sumcoef(4)
    real(kind = rKind)        :: k(3,4)
    integer(kind=iKind)       :: ind(4)
    real(kind=rKind)          :: yk(3)
    !   
    ! IMPLEMENTATION
    !-------------------------------------------------------------------------------------------------------------------    
    !


    
    !
    !------------------------------------------------------------------
    !                          IMPLEMENTATION
    !------------------------------------------------------------------
    !
    nRkStage = 4
    stagecoef = (/ 0. , 0.5, 0.5, 1. /)
    sumcoef   = (/ 1. , 2. , 2. , 1. /)/6.
    ind       = (/ 0,1,1,1 /)

    !
    yout = yin
    k    = 0.
    do ik = 1,nRkStage

       yk  = yin + ds * stagecoef(ik) * k(:,ik - ind(ik) )
       !
       ! evaluate difgl
       !
       call difvglspherical( k(: , ik) ,yk, freq, ds, dry)

       if ( dry ) return
       
       yout = yout + sumcoef(ik) * ds * k(: , ik)
       !
    enddo

    !
  end subroutine rk4step  
  !

  subroutine updateDistTimeFric( dist, time, intFric,vecold , vecnew,freq, dsi)
    !
    use modTools
    use IO
    
    real(kind=rKind)   , intent(in)  :: vecold(3)
    real(kind=rKind)   , intent(in)  :: vecnew(3)
    real(kind=rKind)   , intent(inout) :: dist,time,intFric
    real(kind=rKind)   , intent(in) :: freq
    real(kind=rKind)  , intent(in)   :: dsi
    
    real(kind=rKind)                 :: di(2)
    real(kind=rKind)                 :: xi(2)
    real(kind=rKind)                 :: yi(2)
    real(kind=rKind)                 :: k(2)    
    real(kind=rKind)                 :: ds
    real(kind=rKind)                 :: dt    
    real(kind=rKind)                 :: cg(2)
    real(kind=rKind)                 :: kd(2)
    real(kind=rKind)                 :: jonFric(2)
    real(kind=rKind)                 :: Jac(2)

    xi = (/ vecold(1), vecnew(1) /)
    yi = (/ vecold(2), vecnew(2) /)
    
    call gridInterp2(X,Y,mx,my,Dep,xi(1),yi(1),di(1),1,1) ! Using gridInterp for points is really akward...
    call gridInterp2(X,Y,mx,my,Dep,xi(2),yi(2),di(2),1,1)

    !
    ! Use the "actual" stepsize, rather than the update computational stepsize...
    !
    ds   = dsi !sqrt( ( xi(2) - xi(1) )**2 + ( yi(2) - yi(1) )**2 )
    ds = abs(ds)
    dist = dist + ds

    k(1)  = disperK( freq , di(1) )
    k(2)  = disperK( freq , di(2) )

    kd = k * di
    if ( maxval( kd ) > 10. ) then
       !
       cg = freq/k * .5
       jonFric = 0.
       !
    else
       !
       
       cg    = freq / k * ( .5 + k * di / sinh( 2. * k * di ) )
       Jac   = 1./ (cg / k)
       jonFric = - 0.038 * ( freq   / g / sinh( kd ) ) **2 !* Jac
       !
    endif
    !
    dt = ds * .5 * ( 1/cg(1) + 1/cg(2) )
    time = time + dt
    intFric = intFric + dt * (jonFric(1) + jonFric(2)) *.5
    !
  end subroutine updateDistTimeFric

subroutine difvglspherical( out, in , freq , ds , dry )
    !
    use modTools
    use IO
    
    real(kind=rKind)   , intent(in)  :: in( 3)
    real(kind=rKind)   , intent(out) :: out(3)
    real(kind=rKind)   , intent(in)  :: freq
    real(kind=rKind)   , intent(in)  :: ds
    logical            , intent(out) :: dry

    ! LOCAL
    real(kind=rKind)                 :: sth,cth, dcdx,dcdy,dkdh,dhdx,dhdy
    real(kind=rKind)                 :: w(  5 )   
    real(kind=rKind)                 :: c(  5 )        
    real(kind=rKind)                 :: di( 5 )
    real(kind=rKind)                 :: xi( 5 )
    real(kind=rKind)                 :: yi( 5 )
    real(kind=rKind)                 :: k(  5 )
    real(kind=rKind)                 :: dstep, lat, lon

    dry = .false.

    dstep = 0.5*min( dx_bat , dy_bat )! 1. ! abs(ds/10.)
   ! call log_real('ja',dstep,0)
    xi = [in(1),in(1)+dstep,in(1)-dstep,in(1)      ,in(1)      ]
    yi = [in(2),in(2)      ,in(2)      ,in(2)+dstep,in(2)-dstep]

    call interpDep(xi,yi,di,5)

    if ( di(1) <= drylim ) then
       !
       dry = .true.
       return
       !
    endif

    
    w  = freq
    k  = disperK( w , di )
    c = freq / k 
    dcdx = 0.
    dhdx = 0.    
    if ( di(2) >= drylim .and. di(3) >= drylim) then
       !
       dcdx  = (c(2) - c(3)) / 2. / dstep
       dhdx  = (di(2) - di(3)) / 2. / dstep
       !
    else
       !
       if ( di(2) <= drylim .and. di(3) >=drylim ) then
          !
          dcdx  = (c(1) - c(3)) / dstep
          dhdx  = (di(1) - di(3)) / dstep
          !
       elseif ( di(3) <= drylim .and. di(2) >= drylim ) then
          !
          dcdx  = (c(2) - c(1)) / dstep
          dhdx  = (di(2) - di(1)) / dstep
          !
       else
          !
          dcdx = 0.
          dhdx = 0.
          !
       endif
       !
    endif
    dcdy = 0.
    dhdy = 0.        
    if ( di(4) >= drylim .and. di(5) >= drylim) then
       !
       dcdy  = (c(4) - c(5)) / 2. / dstep
       dhdy  = (di(4) - di(5)) / 2. / dstep       
       !
    else
       !
       if ( di(4) <= drylim .and. di(5) >=drylim ) then
          !
          dcdy  = (c(1) - c(5)) / dstep
          dhdy  = (di(1) - di(5)) / dstep          
          !
       elseif ( di(5) <= drylim .and. di(4) >= drylim ) then
          !
          dcdy  = (c(4) - c(1)) / dstep
          dhdy  = (di(4) - di(1)) / dstep 
          !
       else
          !
          dcdy = 0.
          dhdy = 0.
          !
       endif
       !
    endif    
    !
    dkdh =  an_dkdh(k(1),di(1))
    !
    if (spherical) then
       !
       dcdx = - freq / (k(1)**2) * dkdh * dhdx
       dcdy = - freq / (k(1)**2) * dkdh * dhdy
       !

       
       cth = cos(in(3))
       sth = sin(in(3))
       !
       lon = Lon0 + in(1)
       lat = Lat0 + in(2)       
       
       out(1) = cth / Rearth / cos( lat *pi/180. )
       out(2) = sth / Rearth


       !write( msg , '("angles ", f15.8, f15.8, f7.3, f7.3, f7.3  )' ) Rearth ,lat ,cth,sth,cos(lat)     
       !call wlog( msg , 0 )

       out(3) = 1./c(1)/Rearth *  ( dcdx*sth/cos(lat*pi/180.) - dcdy*cth) - cth * tan(lat*pi/180)/Rearth
       !
    else
       !
       dkdh =  an_dkdh(k(1),di(1))
       dcdx = - freq / (k(1)**2) * dkdh * dhdx
       dcdy = - freq / (k(1)**2) * dkdh * dhdy
       !
       cth = cos(in(3))
       sth = sin(in(3))
       !
       out(1) = cth
       out(2) = sth
       out(3) = 1./c(1) *  ( dcdx*sth - dcdy*cth )       
       !
    endif     
    !
  end subroutine difvglspherical
  !
end module modRayIntegration
