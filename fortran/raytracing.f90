module raytracing

  use iso_c_binding    
  use modPar
  use IO

  
  implicit none

  
  
contains 
  !
  subroutine setnum( frac , maxstep ) bind(c)
    !
    ! This subroutine controls the numerical integration parameters
    !
    real( kind = rkind )   , intent(in) :: frac
    integer( kind = ikind ), intent(in) :: maxstep 

    num_maxstep = maxstep
    num_stepsize_accuracy    = frac

    write( msg , '("\n")')
    write( msg , '("Numerical Parameters:")' )
    write( msg , '("---------------------")' )
    write( msg , '("Accuracy set to: ", f2.0)' ) frac
    write( msg , '("Relative distance: ", f2.0)' ) num_dist_rel
    write( msg , '("Distance ratio: ", f2.0)' ) num_dist_ratio
    call wlog( msg, -1)
    !
    
    !
  end subroutine setnum
  !

  subroutine setbndlist( inbndlist , n ) bind(c)
    !
    ! This subroutine controls the integration domain
    !
    integer(kind=iKind),intent(in)   :: inbndlist(n)
    integer(kind=iKind),intent(in)   :: n
    !
    if ( allocated(bndlist) ) then
       !
       deallocate(bndlist)
       !
    endif
    allocate( bndlist(n) )
    bndlist = inbndlist

    
    nbnd = n
    !
  end subroutine setbndlist
  
  !
  subroutine setdom( inxlim , inylim ) bind(c)
    !
    ! This subroutine controls the integration domain
    !
    real(kind=rKind),intent(in)   :: inxlim(2)
    real(kind=rKind),intent(in)   :: inylim(2)

    call wlog( msg , 0 )
    
    if (allocated(xlim)) then
       !
       deallocate(xlim,ylim)
       !
    endif
    allocate( xlim(2), ylim(2) )

    lon0 = inxlim(1)
    lat0 = inylim(1)
    
    xlim = inxlim - lon0
    ylim = inylim - lat0

    !
  end subroutine setDom
  !

  subroutine calc( matrix, xp,yp,np, freq, nfreq,angles,nang, nsubrays ) bind(c)
    !
    ! Input variables
    !-----------------
    integer( kind = iKind),intent(in)    :: nsubrays
    integer( kind = iKind),intent(in)    :: nang
    integer( kind = iKind)               :: isec
    real(    kind = rKind),intent(in)    :: angles(nang)
    real(    kind = rKind),intent(in)    :: xp(np)
    real(    kind = rKind),intent(in)    :: yp(np)
    integer( kind = iKind),intent(in)    :: np    
    real(    kind = rKind),intent(in)    :: freq(nfreq)
    integer( kind = iKind),intent(in)    :: nfreq
    real(    kind = rKind),intent(inout)    :: matrix(nang,nang,np,nfreq)
    real(    kind = rKind)               :: ang(nang)    
    !
    ! Local variables
    !-----------------
    integer( kind = iKind)               :: nout(nang-1)
    integer( kind = iKind)               :: ip, ifreq, iang

    !
    ! Set angles to -pi < .. < pi
    !
    ang = atan2( sin(angles) , cos(angles) )
    !
    ! Simulation
    !
    write( msg , '("Number of angles ", i04  )' ) nang
    call wlog( msg , 0 )      
    write( msg , '("Number of sub rays ", i04  )' ) nsubrays
    call wlog( msg , 0 )      
    do ifreq = 1 , nfreq
       !
       write( msg , '("Freq ", i04, " out of " ,i04  )' ) ifreq , nfreq
       call wlog( msg , 0 )       
       do ip = 1 , np
          !
          write( msg , '("Point ", i04, " out of " ,i04  )' ) ip , np
          call wlog( msg , 0 )          
          call fullcircle( matrix( 1 , 1 , ip,ifreq) ,nout,xp(ip),yp(ip), &
               freq(ifreq),ang,nang,nsubrays)
          !
          !
       enddo
       !
    enddo
    !
  end subroutine calc
  
  subroutine fullcircle( matrix,nout,xp,yp,freq,angles,nang,maxrays) bind(c)

    use modTools
    !
    !
    ! Input variables
    !-----------------
    integer( kind = iKind),intent(in)    :: maxrays
    integer( kind = iKind),intent(inout) :: nout(nang-1) 
    integer( kind = iKind)               :: nang
    integer( kind = iKind)               :: isec
    integer( kind = iKind)               :: isubray    
    integer( kind = iKind)               :: isubrays    
    integer( kind = iKind)               :: thread_num
    integer( kind = iKind)               :: thread_id
    real(    kind = rKind),intent(in)    :: angles(nang)
    real(    kind = rKind),intent(in)    :: xp
    real(    kind = rKind),intent(in)    :: yp
    real(    kind = rKind),intent(in)    :: freq
    integer(    kind = iKind)            :: iang(2)    
    real(    kind = rKind)               :: dang
    real(    kind = rKind)               :: dang2    
    real(    kind = rKind)               :: fac(2)    
    real(    kind = rKind)               :: inpar(4)
    real(    kind = rKind)               :: out(3)
    real(    kind = rKind)               :: outpar(inumpar)    
    real(    kind = rKind),intent(inout) :: matrix( nang , nang )
    real(    kind = rKind)               :: ccg
    real(    kind = rKind)               :: angdif
    real(    kind = rKind)               :: xpl
    real(    kind = rKind)               :: ypl
    real(    kind = rKind)               :: direction    
    logical                              :: tmp
    integer( kind = iKind )              :: jang
    
    integer( kind = iKind )              :: iDirection
    logical                              :: found
    
    !
    dang = angles(2) - angles(1)
    xpl = xp - lon0
    ypl = yp - lat0

    if ( .not. allocated(bndlist) ) then
       !
       nbnd = 4
       allocate(bndlist(nbnd) )
       bndlist = [1,2,3,4]       
       !
    end if
    !
    !
    !$OMP PARALLEL PRIVATE( isec,isubray,inpar,out,outpar,ccg,iang,fac,jang,dang, iDirection,dang2,direction,found  )
    !
    !
    !$OMP DO
    !
    do isec = 1, nang
       !
       do isubray = 1, maxrays
       do iDirection =  1,2          
          !
          inpar = [ xpl , ypl , angles(isec) - dang/2. + (isubray-1.)*dang/maxrays , freq ]

          if (iDirection == 1) then
             !Backward phase
             direction = -1._rKind         
          else
             !Forward phase
             direction = 1._rKind

          endif
          call tracing( outpar(1), out(1),inpar(1),1,1,1,direction)          
          !
          ! Reason of termination (Rot) == 10 + [ 1..4 ] if terminating on
          ! W=1,E=2,S=3,N=4 boundaries
          !
          if (  inBndlist( int(outpar(iRot)) - 10 ) ) then
             !
             
             ccg = outpar(ic1)*outpar(icg1) / ( outpar(ic0) *  outpar(icg0) )

             ! Limit to -pi <  ... < pi
             outpar(3) = atan2( sin(outpar(3)),cos(outpar(3)))
             !
             if (outpar(3) < angles(1)) then
                !
                outpar(3) = outpar(3) + PI2
                !
             endif
                
             if (outpar(3) > angles(nang)) then
                !
                outpar(3) = outpar(3) - PI2
                !
             endif

             found = .false.
             searchloop: do jang = 1, nang
                !
                if ( outpar(3) >= angles(jang) ) then
                   !
                   if (jang==nang) then
                      iang(1) = nang
                      iang(2) = 1
                      dang2    = angles(1) + pi2 - angles(nang)
                      fac(1) = 1. - ( outpar(3) - angles(nang)    )/dang2
                      fac(2) = 1. - ( angles(1) + pi2 - outpar(3) )/dang2
                      found = .true.
                      exit searchloop
                   endif

                   if (outpar(3) < angles(jang+1)) then
                      !
                      dang2 = angles(jang+1) - angles(jang)
                      iang(1) = jang
                      iang(2) = jang+1
                      fac(1) = 1. - ( outpar(3     ) - angles(jang) )/dang2
                      fac(2) = 1. - ( angles(jang+1) - outpar(3)    )/dang2
                      found = .true.
                      exit searchloop
                      !
                   endif
                   !
                endif
                !
             enddo searchloop

             if ( found ) then
                !
                matrix( iang(2) , isec ) = matrix( iang(2) , isec ) + ccg * fac(2)/ maxrays
                matrix( iang(1) , isec ) = matrix( iang(1) , isec ) + ccg * fac(1)/ maxrays
                !
             else
                !
                ! This branch should never be entered - it is though; not sure why yet
                !
                write( msg , '("Point ", f8.3, f8.3, f8.3  )' ) outpar(3) , angles(1)  , angles(2) 
                call wlog( msg , 0 )     
                !
             endif
             !
          else
             !
             ccg = 0.
             !
          endif          
          !
       enddo
       enddo
       !
    enddo
    !
    !$OMP END DO
    !$OMP END PARALLEL
    !   
  end subroutine fullcircle  
 
  
  !
  subroutine setbat( di , xi, yi,nxi,nyi ) bind(c)
    !
    ! This subroutine controls the bathymetry
    !
    real(kind=rKind), intent(in) :: di(nxi,nyi)
    real(kind=rKind), intent(in) :: yi(nyi)
    real(kind=rKind), intent(in) :: xi(nxi)
    integer(kind=iKind)          :: nxi
    integer(kind=iKind)          :: nyi    
    
    if ( allocated(dep) ) then
       !
       deallocate(dep,x,y)
       !
    endif
    !
    allocate( dep( nxi,nyi),x(nxi), y(nyi) )
    mx = nxi
    my = nyi
    
    dep = di
    x   = xi - lon0
    y   = yi - lat0

    dx_bat = x(2) - x(1)
    dy_bat = y(2) - y(1)
    !
  end subroutine setBat
  !

  subroutine setcont( d ) bind(c)
    !
    real(kind=rKind)   ,intent(in) :: d
    !
    cont_d = d
    cont_set = .true.
    
  end subroutine setCont

  
  subroutine tracing( outpar, out,inpar,nRay, &
       maxout,outstep,direction) bind(c)
    !
    use modTools
    use modRayIntegration
    
    !
    real(kind=rKind)   ,intent(in)  :: inpar(  4, nray )
    real(kind=rKind)   ,intent(in)  :: direction 
    real(kind=rKind)   ,intent(out) :: outpar( inumpar, nray )
    ! 1 ,2 ,3  ,4 ,5  ,6 ,7 ,8 ,9  ,10,11 ,12,13  ,14  ,15  ,16 ,17    ,18
    ! x1,y1,th1,d1,cg1,c1,x0,y0,th0,d0,cg0,c0,time,dist,fric,rot,normal,dnormal
    !    
    integer(kind=iKind),intent(in) :: nray

    real(kind=rKind)   ,intent(out) :: out( 3,maxout,nray)
    integer(kind=iKind)   ,intent(in)  :: maxout
    integer(kind=iKind)   ,intent(in)  :: outstep      

    real(kind=rKind),allocatable :: vecin(:,:)
    real(kind=rKind),allocatable :: vecout(:,:)
    real(kind=rKind)             :: vecin2(3)    
    real(kind=rKind)             :: vecout2(3)    
    real(kind=rKind)             :: tmp(1),tmpK
    real(kind=rKind)             :: ds,dDdx,dDdy,dDdS,xc,yc

    logical                      :: linescross, isdry
    logical                      :: dry(3)

    integer(kind=iKind)         :: istep,iRay,iout, jstep


    real(kind=rKind)                :: dcoarse(3)
    real(kind=rKind)                :: dsol(4)    
    real(kind=rKind)                :: norm(3)
    real(kind=rKind)                :: acc(4),dist,distf
    character(len=100)              :: message
    !
    ! INIT
    !    
    allocate(vecin (3,nray)); vecin  = inpar(1:3,:)
    allocate(vecout(3,nray)); vecout = inpar(1:3,:)

    
    !
    ! Do for each ray
    !  
    do iRay  = 1 , nRay
       !
       outpar( : , iray )  = 0.
       call interpDep([inpar(1,iray)],[inpar(2,iray)],tmp(1),1)
       tmpK = disperK( inpar(4,iray), tmp(1) )

       if (tmp(1) < 0.) then
          !
          write(  message , '("WARNING: depth negative at starting point")' )
          write(*,*) message
          call wlog( message , 0 )
          cycle
          !
       endif
       
       outpar( ix0 , iray )  = inpar(1,iray)
       outpar( iy0 , iray )  = inpar(2,iray)
       outpar( ith0, iray )  = inpar(3,iray)
       outpar( id0 , iray )  = tmp(1)
       outpar( ic0 , iray )  = inpar( 4,iray) / tmpK
       outpar( icg0 , iray ) = outpar( ic0 , iray )  * &
             ( .5 +tmpK * tmp(1) / sinh( 2. * tmpK * tmp(1) ) )
       
       iout = 1
       if (maxout > 1) then
          !
          out( 1, iout,iray ) = inpar(1,iray)
          out( 2, iout,iray ) = inpar(2,iray)
          out( 3, iout,iray ) = inpar(3,iray)
          iout = 2
          !
       endif
       !
       linescross = .false.

       msg = ''

       outpar( iROT , iRay) = -1 !Default is: reached maximum number of steps...
       RAYLOOP: do istep = 1 , num_maxstep
          !
          ! Are we still in the computational domain?
          !
          if ( inDomain( vecin(1,iRay),vecin(2,iRay)  ) >= 0 ) then
             !
             outpar( iROT,iRay) = 10 + inDomain( vecin(1,iRay),vecin(2,iRay)  )
             !Termination is: hit bounding box of bathy/region
             exit RAYLOOP
             !
          endif
          !
          ! Ok, so how large will our stepsize be?
          !       
          if (  num_stepSizeControlMethod == 0 ) then
             !
             ! Estimation based on bottom grid
             !
             ds = direction * stepsize( vecin(1,iRay), vecin(2,iRay), inpar(4,iRay),vecin(3,iRay) )             
             call rk4step( vecin(1,iRay),vecout(1,iRay),inpar(4,iRay), ds  , dry(1))
             !
          elseif (  num_stepSizeControlMethod == 1 ) then
             !
             ! Adaptive Runge Kutta as described in numerical recipes
             !
             !
             !The bottom gradients are in degrees, switch to radians
             !implies the additional Jacobian Factor pi/180
             ds  = direction * min(pi/180*dx_bat* Rearth*cos(lat0*pi/180.), &
                  pi/180*dy_bat * Rearth) * num_stepsize_fracBatGrid
             !
             
             ACCURACYLOOP: do jstep = 1 , 10
                !
                ! The coarse step
                !
                vecin2 = vecin(:,iRay)
                vecout2 = vecin2
                call rk4step( vecin2(1),vecout2(1),inpar(4,iRay),  ds,dry(1)  )
                dcoarse = vecin2 - vecout2
                !
                ! 2 finer steps
                !
                call rk4step( vecin2(1),vecout(1,iRay),inpar(4,iRay), ds /2.,dry(2) )
                vecin2 = vecout(:,iRay)
                call rk4step( vecin2(1),vecout(1,iRay),inpar(4,iRay), ds /2.,dry(3)  )
                !
                ! Calculate accuracy
                !
                dsol(1:3) = vecout2( :  ) - vecout( : , iRay )


                dsol(4) = vecin(3,Iray)  - vecout2(3)
                !
                call interpDep([vecin(1,iRay)],[vecin(2,iRay)],tmp(1),1)

                !Compare accuracy with size of the bottom grid
                acc(1) = num_stepsize_accuracy * min(dx_bat* Rearth*cos(lat0*pi/180.),dy_bat * Rearth)
                acc(2) = acc(1)
                acc(3) = pi2/3600.
                acc(4) = pi2/72.
                !
                !num_stepsize_accuracy * min(  sqrt(g * tmp(1) ) /  inpar(4,iRay)  *pi2 , abs(ds) )

                if (dry(1)  .or. dry(2) .or. dry(3)  ) then
                   !
                   ds = ds / 2.
                   cycle accuracyloop
                   !
                endif
                !
                ! new stepsize
                !ds = ds * abs(minval( (abs( acc) / (abs(dsol)) )**(0.2) ))
                
                !
                ! If the prescribed accuracy hit...                !
                !
                ! The accuracy checks for (1) if the difference between solutions
                ! is smaller then 0.1 degree, and (2) that the ray did not turn more than
                ! 5 degrees over a step.
                if ( dsol(3)**2 <= acc(3)**2 .and. &
                     dsol(4)**2 <= acc(4)**2 ) then
                   !
                   ! then see if we can enlarge step?
                   !

                   !
                   ! and exit
                   !
                   exit accuracyloop
                   !
                else
                   !
                   ds = ds / 2.
                   !
                endif
                !
             enddo ACCURACYLOOP
             !
             ! And do free Richardson's extrapolation
             !
             !vecout( : , iRay ) = vecout( : , iRay ) + dsol(1:3)/15.
             !
          endif
          
          !
          isdry = .false.
          if (dry(1)  .or. dry(2) .or. dry(3)  ) then
             !
             ! So we are on land ...
             !
             outpar(iRot,iRay) = 100 !Termination is: hit land
             vecout(1,iRay) = vecin(1,iRay)
             vecout(2,iRay) = vecin(2,iRay)
             vecout(3,iRay) = vecin(3,iRay)
             isdry = .true.
             !
          else
             !
             !
             ! Did we cross the boundary?
             !
             call checkcross( xc, yc,linescross , vecin(1,iRay),vecout(1,iRay) )
             !

             if (linescross ) then
                !             
                outpar(iRot,iRay) = 1 !Termination is: hit prescribed line

                dist  = sqrt(  (xc - vecin(1,iRay))**2 +  (yc - vecin(2,iRay))**2 )
                distf = sqrt(  (vecout(1,iRay) - vecin(1,iRay))**2 + &
                     (vecout(2,iRay) - vecin(2,iRay))**2 )

                vecout(1,iRay) =  xc
                vecout(2,iRay) =  yc
                !
                ! Interpolate the angle..., hijack the already declared variables for this :)
                !             
                vecout(3,iRay) =  vecin(3,iRay) +dist * (vecout(3,iRay) - vecin(3,iray)) / distf             
                !
             endif
             
             !
          endif
          
         
          
          if (maxout > 1 .and. iout <= maxout .and. ( mod(istep,outstep )==0 .or. linescross .or. isdry) ) then
             !
             out( 1, iout,iray ) = vecout(1,iRay)
             out( 2, iout,iray ) = vecout(2,iRay)
             out( 3, iout,iray ) = vecout(3,iRay)

             iout = iout+1
             !
          endif


          call updateDistTimeFric(  outpar( iDist , iRay), outpar( iTime , iRay), outpar( iFric , iRay),  &
                                    vecin(1,iRay) ,vecout(1,iRay),inpar(4,iRay), ds)
          vecin(:,iRay) = vecout(:,iRay)    
          if (linescross) exit RAYLOOP


          !
       enddo RAYLOOP
       !
       outpar(1,iRay)   = vecin(1,iRay)
       outpar(2,iRay)   = vecin(2,iRay)
       outpar(3,iRay)   = vecin(3,iRay)
       !

       if ( linescross ) then
          !
          ! If we crossed the desired depth contour, set the depth to that contour
          !
          tmp(1) = cont_d
          !
       else
          !
          call interpDep([vecin(1,iRay)],[vecin(2,iRay)],tmp(1),1)
          !
       endif
       
       tmpK = disperK( inpar(4,iray), tmp(1) )       
       outpar( id1 , iray )  = tmp(1)
       outpar( ic1 , iray )  = inpar( 4,iray) / tmpK
       outpar( icg1 , iray ) = outpar( ic1 , iray )  * &
            ( .5 +tmpK * tmp(1) / sinh( 2. * tmpK * tmp(1) ) )
       !
       call dDepdS(  vecin(1,iRay),  vecin(2,iRay), dDdx , dDdy, dDdS, abs(dx_bat), abs(dy_bat) )
       outpar( inormal , iray )  = atan2( dDdy , dDdx )
       outpar( idnormal , iray ) = atan2( sin(outpar(3,iRay)) , cos(outpar(3,iRay)) ) - outpar( inormal , iray )
       !
       if ( outpar( idnormal , iray ) >   pi ) outpar( idnormal , iray ) = outpar( idnormal , iray ) - pi2
       if ( outpar( idnormal , iray ) < - pi ) outpar( idnormal , iray ) = outpar( idnormal , iray ) + pi2
       !
    enddo
    !
  end subroutine tracing
  !

  !
  real(kind=rKind) function stepsize( xi , yi, freq,theta)
    !
    !
    use modTools
    !
    real(kind=rKind),intent(in) :: xi
    real(kind=rKind),intent(in) :: yi
    real(kind=rKind),intent(in) :: freq,theta
    real(kind=rKind)            :: dxb
    real(kind=rKind)            :: dyb,ds

    real(kind=rKind)            :: d(1),k, dCy,dCx,dCs,dThds

    call interpDep([xi],[yi],d,1)
    k = disperK( freq , d(1) )

    dxb = x(2)-x(1)
    dyb = y(2)-y(1)
    !
    ! Stepsize is at most a fraction of the bathymetric grid
    !
    stepsize = min(dxb,dyb) * num_stepsize_fracBatGrid   ! pi2 / k * num_frac

    !
    ! Stepsize is at most a fraction of the wavenumber
    !
    !stepsize = min( stepsize ,  pi2 / k * num_stepsize_fracK )
    !ds = stepsize
    !dCx = parderK( xi,yi,freq,stepsize,0._rKind)
    !dCy = parderK( xi,yi,freq,0._rKind,stepsize)
    !dCs = sqrt( dCx ** 2 + dCy ** 2 )

    !
    ! Stepsize is at most a fraction of the wavenumber
    !
    !stepsize = min( stepsize , num_frac_k / dcs )

    !
    ! Only xxx degree per step
    !
    dthds = ( dCx*sin(theta) - dCy*cos(theta) )
    if ( abs(dthds) > 0. ) then
       !
       ds = num_stepsize_maxdegPerStep / 180. * pi2 / abs(dthds)
       !
       if ( ds < stepsize ) then
          !
          stepsize = ds
          !
       endif
       !
    endif
    !

    !call log_int( 'numTh', stepsizelimit(3),0 )
    !call log_int( 'numTo', stepsizelimit(4),0 )  
    !
  end function stepsize
  !

  !
  subroutine localcontour( xi , yi , xc , yc , dc, found)
    !
    ! Find the local contour line segment; note, this algorithms fails for saddle points
    !
    use modTools
    !    
    logical                , intent(out) :: found
    real( kind = rKind )   , intent(in ) :: xi(2)  ,yi(2)
    real( kind = rKind )   , intent(out) :: xc(2) , yc(2)
    real( kind = rKind )   , intent(in)  :: dc

    real( kind = rKind )                 :: dx , dy, ds
    real( kind = rKind )                 :: xp(4) , yp(4), dp(4), xm ,ym
    integer( kind = iKind )              :: iSides, nSides, jm, jp

    !
    ! The local points may not contain the sought for depth, in that case found
    ! will return as false
    !
    found = .false.
    dx = abs( xi(2) - xi(1) )
    dy = abs( yi(2) - yi(1) )

    ds = max( dx , dy )
    xm = sum(xi)/2.
    ym = sum(yi)/2.
    !
    
    !
    ! The Interpolation points
    xp = [ xm - 0.5*ds , xm + 0.5 * ds,  xm + 0.5 * ds , xm - 0.5 * ds ]
    yp = [ ym - 0.5*ds , ym - 0.5 * ds,  ym + 0.5 * ds , ym + 0.5 * ds ]
    !
    call  interpDep(xp,yp,dp,4)
    
    !
    ! If the contour depth isn't within the parameter, return...
    !
    if ( .not. ( minval( dp) <= dc .and. maxval( dp ) >= dc ) ) return

    !
    nsides = 0
    SIDESLOOP: do isides = 1 , 4
       !
       jm = isides
       jp = isides + 1
       if (jp > 4) jp = 1
       
       if ( maxval( [ dp(jm) , dp(jp) ] ) >= dc   .and. &
            minval( [ dp(jm) , dp(jp) ] ) <  dc ) then
          !

          nsides = nsides + 1
          dx = xp(jp) - xp(jm)
          dy = yp(jp) - yp(jm)
          !
          if ( abs( dx ) > 0. ) then
             !
             xc( nsides ) = dx * ( dc - dp(jm) )/ ( dp( jp ) - dp(jm) ) + xp(jm)
             yc( nsides ) = yp(jm)
             !
          else
             !
             yc( nsides ) = dy * ( dc - dp(jm) )/ ( dp( jp ) - dp(jm) ) + yp(jm)
             xc( nsides ) = xp(jm)             
             !
          endif
          !
          if (nsides == 2) then
             !
             found = .true.
             exit SIDESLOOP
             !
          endif
          !
       endif
       !
    enddo SIDESLOOP
    !
    !    
  end subroutine localcontour
  !
  !
  subroutine checkcross( xc , yc, linescross , pa1,pa2 )
    !
    ! Here we check if the step crosses a depth contour
    !
    use modTools
    
    logical             , intent(inout) :: linescross
    real( kind = rKind) , intent(in)    :: pa1(2)
    real( kind = rKind) , intent(in)    :: pa2(2)
    real( kind = rKind ), intent(out)   :: xc
    real( kind = rKind ), intent(out)   :: yc    
    !
    real( kind = rKind )                :: d
    real( kind = rKind )                :: dp(2)
    real( kind = rKind )                :: xp(2)
    real( kind = rKind )                :: yp(2)
    real( kind = rKind )                :: cont_x(2)
    real( kind = rKind )                :: cont_y(2)        
    
    integer( kind = iKind )             :: iseg
    logical                             :: found

    xp(1) = pa1(1)
    xp(2) = pa2(1)
    yp(1) = pa1(2)
    yp(2) = pa2(2)        
    xp(1) = pa1(1)
    xp(2) = pa2(1)
    
    call interpDep(xp,yp,dp,2)
    d = minval( dp )

    linescross = .false.
    if (  (d - cont_d ) > 0.7 * cont_d ) return

    found = .false.
    call localcontour( xp , yp , cont_x , cont_y , cont_d, found)
    !
    if ( found ) then
       !            
       
       call  lineCrossing( linescross, xc,yc, &
            xp(1)  , xp(2), &
            yp(1)  , yp(2), &
            cont_x(1), cont_x(2), &
            cont_y(1), cont_y(2) )
       
       if ( linescross ) then
          !
          return
          !
       endif
       !
    endif
    !
    ! Did we cross the boundary?
    !
    ! contloop: do iseg = 1,cont_n-1
    !    !
    !    call  lineCrossing( linescross, xc,yc, &
    !         xp(1)  , xp(2), &
    !         yp(1)  , yp(2), &
    !         cont_xy(1,iSeg), cont_xy(1,iSeg+1), &
    !         cont_xy(2,iSeg), cont_xy(2,iSeg+1) )
    !    !
    !    if ( linescross ) then
    !       !
    !        call log_real( 'ja' , cont_x(1), 1 )
    !       return
    !       !
    !    endif
    !    !
    ! enddo contloop
    ! !

    !
    if ( d < cont_d ) then
       !
       linescross = .true.
       xc = xp(2)
       yc = yp(2)
       return
       !
    endif
    !
  end subroutine checkcross
  !

  !
  subroutine outputspec( Eout ,Ebnd,fbnd,dirbnd, dtbnd,nfbnd,ndirbnd, ntbnd, data ,fdata, &
       dirdata, nout , nray,ndir,nf,np,numray,kdleak) bind(c)
    !
    use modtools, only:disperK
    
    real( kind = rKind )  , intent(out):: Eout( 8,ndir , nfbnd,np )
    real( kind = rKind )  , intent(in) :: Ebnd( ndirbnd , nfbnd ,ntbnd )
    integer( kind = iKind), intent(in) :: ntbnd
    real( kind = rKind )  , intent(in) :: dtbnd
    real( kind = rKind )  , intent(in) :: dirbnd( ndirbnd )
    integer(kind=iKind)   , intent(in) :: ndirbnd    
    real( kind = rKind )  , intent(in) :: fbnd( nfbnd )
    integer(kind=iKind)   , intent(in) :: nfbnd
    real( kind = rKind )  , intent(in) :: data( nout , nray , ndir , nf , np )
    real( kind = rKind )  , intent(in) :: fdata( nf )
    real( kind = rKind )  , intent(in) :: dirdata( ndir )
    real( kind = rKind )  , intent(in) :: kdleak(2)    
    integer(kind=iKind)   , intent(in) :: nout,nray,ndir,nf,np
    integer(kind = iKind) , intent(in) :: numray( ndir , nf , np )

    integer(kind = iKind)              :: ip , jf, idir, iRay, it(2)
    integer(kind = iKind)              :: ifloc(2),iffac(2),df,iType,iang(2),ifdata,ifout,dt
    real( kind = rKind )               :: K , KD,KD2, Dep, Jac, Fric, dTheta
    real( kind = rKind )               :: dThetaTotal,ang1,ang2,dang,fac(4),Eintp, Dep2
    real( kind = rKind )               :: cs,sn,da
!    real( kind = rKind )               :: depmat( 8 , ndirbnd , ndirbnd , nfbnd , ntbnd , np ) !The Dependency Matrix...   
    !
    dang = dirbnd(2) - dirbnd(1)
    Eout = 0.
    !
    do ip = 1 ,np
       !
       do ifout = 1 , nfbnd          
          !
          ifdata = minloc( fdata - fbnd(ifout),1 )
          !
          do idir = 1 , ndir
             !
             dThetaTotal = abs( data( 9 , 1 , idir , ifdata , ip ) &
                  - data( 9 , numray(idir ,ifdata ,  ip) , idir , ifdata , ip ) )
             !
             do iRay = 1 , numray( idir , ifdata , ip )
                !
                Dep  = data( 22 , iray , idir, ifdata , ip )
                Dep2 = data( 4 , iray , idir, ifdata , ip )                
                !
                ! Relative depth of the forward traced ray (i.e. did the
                ! ray exit the domain at deep water?)
                !
                K  = disperK( fbnd(ifout), data(22,iray,idir,ifdata,ip) )
                KD = K * data(22,iray , idir , ifdata ,ip )
                !
                ! Relative depth of the backward traced ray (i.e., did
                ! the ray come from deep water?)
                !
                K  = disperK( fbnd(ifout), data(4,iray,idir,ifdata,ip) )
                KD2 = K * data(4,iray , idir , ifdata ,ip )                
                !
                ! Did the ray originate nearshore ( irot = 1 ) or off-
                ! shore ( irot = 0 ), and did it end nearshore (irotf = 1)
                ! or offshore ( irotf = 0 )
                !
                if (       data( irot , iray,idir,ifdata,ip ) == 1 &
                     .and. data( irotf, iray,idir,ifdata,ip ) == 1 ) then
                   !
                   ! Trapped waves
                   !
                   iType = 1
                   !
                elseif (   data( irot , iray,idir,ifdata,ip ) == 1 &
                     .and. data( irotf, iray,idir,ifdata,ip ) == 0 &
                     .and. ( ( kd  > kdleak(1) ) .or. ( dep  > kdleak(2) ) ) ) then
                   !
                   ! "Free" or "Leaky" waves
                   !
                   iType = 2
                   !
                elseif (   data( irot , iray,idir,ifdata,ip ) == 0 &                     
                     .and. ( ( kd2 > kdleak(1) ) .or. ( Dep2 > kdleak(2) ) ) ) then
                   ! .and. data( irotf, iray,idir,ifdata,ip ) == 1 &
                   ! "Incident" waves from distant shores
                   !
                   itype = 4
                   !                   
                else
                   !
                   ! Can't determine waves
                   !
                   iType = 3
                   !
                   if (.not. data( irot , iray,idir,ifdata,ip ) == 1) then
                      !
                      cycle
                      !
                   endif
                   !
                endif
                !
                ! Jacobian
                !
                Jac =  data( icg1 , iray,idir,ifdata,ip ) * data( ic1 , iray,idir,ifdata,ip ) &
                   / ( data( icg0 , iray,idir,ifdata,ip ) * data( ic0 , iray,idir,ifdata,ip ) )
                !
                ! Friction
                !
                Fric = exp( data( iFric , iRay , idir , ifdata , ip ) )                
                !
                ! Integration factor
                !
                dTheta = 0.
                !
                if  (Iray > 1) then
                   !
                   dTheta =  abs( data( 9 , iray - 1 , idir , ifdata , ip ) &
                  - data( 9 , iray , idir , ifdata , ip ) ) / 2.
                   !
                endif
                !
                
                !
                if  (Iray < numray( idir ,ifdata ,  ip ) ) then
                   !
                   dTheta = dTheta + abs( data( 9 , iray + 1 , idir , ifdata , ip ) &
                  - data( 9 , iray , idir , ifdata , ip ) ) / 2.
                   !
                endif
                !
                dTheta = dTheta / dThetaTotal
                !
                !  Interpolate in directions...
                !-------------------------------

                !
                ! Calculate the angle with the normal...
                !
                ang1 = data( inormal , iRay, idir ,ifdata , ip )
                ang2 = data( ith1    , iRay, idir ,ifdata , ip )

                cs  = cos( ang1 ) * cos( ang2 ) + sin( ang1 ) * sin( ang2 )
                sn  = sin( ang1 ) * cos( ang2 ) - cos( ang1 ) * sin( ang2 )
                da  = atan2( sn , cs )
                !
                ! find surrounding points
                !
                if ( da < dirbnd(1) .or. da > dirbnd(ndirbnd) ) then
                   !
                   if ( itype .ne. 4)  cycle
                   !
                end if
                !
                if ( da > 0.5 *pi .or. da < -0.5*pi ) then
                   !
                   if (da>0.5*pi) then
                      da = pi-da
                   endif

                   if (da<-0.5*pi) then
                      da =-pi-da
                   endif                   
                   !if ( itype .ne. 4)  cycle
                   !
                endif
                !
                da      = da - dirbnd(1)
                iang(1) = min( floor(   da / dang ) + 1, ndirbnd )
                iang(2) = min( iang(1) + 1 , ndirbnd)
                da      = mod( da , dang )
                !
                ! Interpolation in angles
                !
                fac(1)= 1. - da / dang
                fac(2)= da / dang
                !
                ! Interpolation in time (if relevant)
                !
                if (ntbnd > 1) then
                   !
                   dt = data( iTime , iRay, idir,ifdata,ip )
                   it(1) = floor( dtbnd / dt )
                   it(2) = it(1) + 1
                   fac(4) = ( dtbnd /dt ) - real(it(1)-1,kind=rKind)
                   fac(3) = 1 - fac(4)
                   !
                   if ( it(2) > ntbnd ) then
                      !
                      it(1) = ntbnd
                      it(2) = ntbnd
                      fac(3) = 1.
                      fac(4) = 0.
                      !
                   endif
                   !
                   it(1) = ntbnd - it(1)
                   it(2) = ntbnd - it(2)                   
                   !
                   Eintp = &
                        + fac(3) * ( Ebnd( iang(1) , ifout , it(1) ) * fac(1)   &
                        +            Ebnd( iang(2) , ifout , it(1) ) * fac(2) ) &
                        + fac(4) * ( Ebnd( iang(1) , ifout , it(2) ) * fac(1)   &
                        +            Ebnd( iang(2) , ifout , it(2) ) * fac(2) )
                   !
                else
                   ! 
                   Eintp = Ebnd( iang(1) , ifout , 1) * fac(1) + Ebnd( iang(2) , ifout , 1) * fac(2)
                   !
                endif
                !
                if (itype == 4) then
                   !
                   Eintp = 1.
                   !
                endif
                ! 
                ! Calculate Output
                !------------------------------
                !
                ! With Friction
                !
                Eout( iType ,idir ,ifout ,ip) = Eout( iType ,idir ,ifout ,ip) + &
                     Eintp * Jac * dTheta * Fric
                !
                ! Without Friction
                !
                Eout( iType+4,idir ,ifout ,ip) = Eout( iType+4,idir ,ifout ,ip) + &
                Eintp * Jac * dTheta
                !
             enddo
             !
          end do !idir
          !
       end do !ifout
       !
    end do !ip
    !
  end subroutine outputspec
  !

  !
  subroutine outputdepmat( Depmat ,fbnd,dirbnd, dtbnd,nfbnd,ndirbnd, ntbnd, data ,fdata, &
       dirdata, nout , nray,ndir,nf,np,numray,kdleak,cffac) bind(c)
    !
    use modtools, only:disperK
    
    real( kind = rKind )  , intent(out):: Depmat( 8,ndirbnd,ntbnd,ndir ,nfbnd,np )
    integer( kind = iKind), intent(in) :: ntbnd
    real( kind = rKind )  , intent(in) :: dtbnd
    real( kind = rKind )  , intent(in) :: dirbnd( ndirbnd )
    integer(kind=iKind)   , intent(in) :: ndirbnd    
    real( kind = rKind )  , intent(in) :: fbnd( nfbnd )
    integer(kind=iKind)   , intent(in) :: nfbnd
    real( kind = rKind )  , intent(in) :: data( nout , nray , ndir , nf , np )
    real( kind = rKind )  , intent(in) :: fdata( nf )
    real( kind = rKind )  , intent(in) :: dirdata( ndir )
    real( kind = rKind )  , intent(in) :: kdleak(2)
    real( kind = rKind )  , intent(in) :: cffac(2)       
    integer(kind=iKind)   , intent(in) :: nout,nray,ndir,nf,np
    integer(kind = iKind) , intent(in) :: numray( ndir , nf , np )

    integer(kind = iKind)              :: ip , jf, idir, iRay, it(2), jj, kk
    integer(kind = iKind)              :: ifloc(2),iffac(2),df,iType,iang(2),ifdata,ifout,dt
    real( kind = rKind )               :: K , KD,KD2, Dep, Jac, Fric(2), dTheta
    real( kind = rKind )               :: dThetaTotal,ang1,ang2,dang,fac(4),Eintp, Dep2
    real( kind = rKind )               :: cs,sn,da
!    real( kind = rKind )               :: depmat( 8 , ndirbnd , ndirbnd , nfbnd , ntbnd , np ) !The Dependency Matrix...   
    !
    dang = dirbnd(2) - dirbnd(1)
    Depmat = 0.
   ! call log_real( 'kdleak1',kdleak(1),0 )
   ! call log_real( 'kdleak2',kdleak(2),0 )    
    !
    do ip = 1 ,np
       !
       do ifout = 1 , nfbnd          
          !
          ifdata = minloc( abs(fdata - fbnd(ifout)),1 )
          !
          do idir = 1 , ndir
             !
             dThetaTotal = abs( data( 9 , 1 , idir , ifdata , ip ) &
                  - data( 9 , numray(idir ,ifdata ,  ip) , idir , ifdata , ip ) )
             !
             do iRay = 1 , numray( idir , ifdata , ip )
                !
                Dep  = data( 22 , iray , idir, ifdata , ip )
                Dep2 = data( 4 , iray , idir, ifdata , ip )                
                !
                ! Relative depth of the forward traced ray (i.e. did the
                ! ray exit the domain at deep water?)
                !
                K  = disperK( fbnd(ifout), data(22,iray,idir,ifdata,ip) )
                KD = K * data(22,iray , idir , ifdata ,ip )
                !
                ! Relative depth of the backward traced ray (i.e., did
                ! the ray come from deep water?)
                !
                K  = disperK( fbnd(ifout), data(4,iray,idir,ifdata,ip) )
                KD2 = K * data(4,iray , idir , ifdata ,ip )                
                !
                ! Did the ray originate nearshore ( irot = 1 ) or off-
                ! shore ( irot = 0 ), and did it end nearshore (irotf = 1)
                ! or offshore ( irotf = 0 )
                !
                if (       data( irot , iray,idir,ifdata,ip ) == 1 &
                     .and. data( irotf, iray,idir,ifdata,ip ) == 1 ) then
                   !
                   ! Trapped waves
                   !
                   iType = 1
                   !
                elseif (   data( irot , iray,idir,ifdata,ip ) == 1 &
                     .and. data( irotf, iray,idir,ifdata,ip ) == 0 &
                     .and. ( ( kd  > kdleak(1) ) .or. ( dep  > kdleak(2) ) ) ) then
                   !
                   ! "Free" or "Leaky" waves
                   !
                   iType = 2
                   !
                elseif (   data( irot , iray,idir,ifdata,ip ) == 0 &                     
                     .and. ( ( kd2 > kdleak(1) ) .or. ( Dep2 > kdleak(2) ) ) ) then
                   ! .and. data( irotf, iray,idir,ifdata,ip ) == 1 &
                   ! "Incident" waves from distant shores
                   !
                   itype = 4
                   !                   
                else
                   !
                   ! Can't determine waves, two options, I can't determine where they come from,
                   ! or where they go to (in case of trapped)
                   iType = 3
                   !
                   ! If I do not know where they come from, I can't assign a energy density,
                   ! so just skip the entry, else, we know it came from the coast, so lets
                   ! just log it under three -> either trapped or leaky radiated energy
                   !
                   if (.not. (data( irot , iray,idir,ifdata,ip ) == 1) ) then
                      !
                      cycle
                      !
                   endif
                   !
                endif
                !
                ! Jacobian
                !
                Jac =  data( icg1 , iray,idir,ifdata,ip ) * data( ic1 , iray,idir,ifdata,ip ) &
                   / ( data( icg0 , iray,idir,ifdata,ip ) * data( ic0 , iray,idir,ifdata,ip ) )
                !
                ! Friction
                !
                Fric = exp( cffac * data( iFric , iRay , idir , ifdata , ip ) )                
                !
                ! Integration factor
                !
                dTheta = 0.
                !
                if  (Iray > 1) then
                   !
                   dTheta =  abs( data( 9 , iray - 1 , idir , ifdata , ip ) &
                  - data( 9 , iray , idir , ifdata , ip ) ) / 2.
                   !
                endif
                !
                
                !
                if  (Iray < numray( idir ,ifdata ,  ip ) ) then
                   !
                   dTheta = dTheta + abs( data( 9 , iray + 1 , idir , ifdata , ip ) &
                  - data( 9 , iray , idir , ifdata , ip ) ) / 2.
                   !
                endif
                !
                dTheta = dTheta / dThetaTotal
                !
                !  Interpolate in directions...
                !-------------------------------

                !
                ! Calculate the angle with the normal...
                !
                ang1 = data( inormal , iRay, idir ,ifdata , ip )
                ang2 = data( ith1    , iRay, idir ,ifdata , ip )

                cs  = cos( ang1 ) * cos( ang2 ) + sin( ang1 ) * sin( ang2 )
                sn  = sin( ang1 ) * cos( ang2 ) - cos( ang1 ) * sin( ang2 )
                da  = atan2( sn , cs )
                !
                ! find surrounding points
                !
                if ( da < dirbnd(1) .or. da > dirbnd(ndirbnd) ) then
                   !
                   if ( itype .ne. 4)  cycle
                   !
                end if
                !
                if ( da > 0.5 *pi .or. da < -0.5*pi ) then
                   !
                   if (da > 0.5 *pi) then
                      da = .5*pi - da;
                   else
                      da = -.5*pi - da;
                   endif
                   !if ( itype .ne. 4)  cycle
                   !
                endif
                !
                da      = da - dirbnd(1)
                iang(1) = min( floor(   da / dang ) + 1, ndirbnd )
                iang(2) = min( iang(1) + 1 , ndirbnd)
                da      = mod( da , dang )
                !
                ! Interpolation in angles
                !
                fac(1)= 1. - da / dang
                fac(2)= da / dang
                !
                ! Interpolation in time (if relevant)
                !
                if (ntbnd > 1) then
                   !
                   dt = data( iTime , iRay, idir,ifdata,ip )
                   it(1) = floor( dt / dtbnd )
                   it(2) = it(1) + 1
                   fac(4) = ( dt /dtbnd ) - real(it(1),kind=rKind)
                   fac(3) = 1. - fac(4)
                   !
                   if ( it(2) == ntbnd ) then
                      !
                      it(1) = ntbnd-1
                      it(2) = ntbnd-1
                      fac(3) = 1.
                      fac(4) = 0.
                      !
                   endif
                   !
                   it(1) = ntbnd - it(1)
                   it(2) = ntbnd - it(2)                   
                   !
                 !  Eintp = &
                 !       + fac(3) * ( Ebnd( iang(1) , ifout , it(1) ) * fac(1)   &
                 !       +            Ebnd( iang(2) , ifout , it(1) ) * fac(2) ) &
                 !       + fac(4) * ( Ebnd( iang(1) , ifout , it(2) ) * fac(1)   &
                 !       +            Ebnd( iang(2) , ifout , it(2) ) * fac(2) )
                   !
                else
                   !
                   it     = ntbnd
                   fac(3) = 1.
                   fac(4) = 0.
!                   Eintp = Ebnd( iang(1) , ifout , 1) * fac(1) + Ebnd( iang(2) , ifout , 1) * fac(2)
                   !
                endif
                !
                if (itype == 4) then
                   !
                   Eintp = 1.
                   !
                endif
                ! 
                ! Calculate Output
                !------------------------------
                !
                ! With Friction
                !
                do jj = 1 , 2
                   !
                   do kk = 1 , 2
                      !                      
                      Depmat( iType , iang(jj) , it(kk) , idir , ifout , ip ) = &
                         Depmat( iType , iang(jj) , it(kk) , idir , ifout , ip ) + &
                         fac(2 + kk ) * fac( jj ) * Jac * dTheta * Fric(1)

                      Depmat( iType +4 , iang(jj) , it(kk) , idir , ifout , ip ) = &
                         Depmat( iType +4 , iang(jj) , it(kk) , idir , ifout , ip ) + &
                         fac(2 + kk ) * fac( jj ) * Jac * dTheta * Fric(2)
                      !
                   enddo
                   !
                enddo
                !
             enddo
             !
          end do !idir
          !
       end do !ifout
       !
    end do !ip
    !
  end subroutine outputdepmat
  !

  !
  subroutine getspec(E , Depmat, Ebnd, Ebndo ,nfbnd,ndirbnd, ntbnd, ndir, np) bind(c)
    !

    real( kind = rKind )  , intent(out):: E( 8, ndir , nfbnd, np )
    real( kind = rKind )  , intent(in) :: Depmat( 8,ndirbnd,ntbnd,ndir , nfbnd,np )        
    real( kind = rKind )  , intent(in) :: Ebnd( ndirbnd , nfbnd,ntbnd )
    real( kind = rKind )  , intent(in) :: Ebndo( ndirbnd , nfbnd,ntbnd )
    integer( kind = iKind), intent(in) :: ntbnd
    integer(kind=iKind)   , intent(in) :: ndirbnd
    integer(kind=iKind)   , intent(in) :: nfbnd
    integer(kind=iKind)   , intent(in) :: np
    integer(kind=iKind)   , intent(in) :: ndir    

    integer(kind = iKind)              :: ip , jf, idir,idirbnd, itype, ifout, it
    real(kind=rKind)                   :: fac(8)
    !
    fac = [ 1. , 1. , 1.,0.,1.,1.,1.,0. ]
    !
    E = 0.

        !
    !
 !  !$OMP PARALLEL PRIVATE( ip,ifout,idir,it,idirbnd,itype,fac,Ebnd,Ebndo )
 !   !$OMP DO
    !
    do ip = 1 ,np
       !
       do ifout = 1 , nfbnd          
          !
          do idir = 1 , ndir
             !
             do it = 1 , ntbnd
                !
                do idirbnd = 1 , ndirbnd             
                   !
                   do itype = 1,8
                      !

                      E( iType, idir , ifout, ip ) = E( iType, idir , ifout, ip ) + &
                           fac(itype) * Ebnd(idirbnd,ifout,it)*depmat(iType,idirbnd,it,idir,ifout,ip) + &
                           ( 1. - fac(itype) ) * Ebndo(idirbnd,ifout,it)*depmat(iType,idirbnd,it,idir,ifout,ip)
                        
                      !
                   enddo
                   !
                end do
                !
             enddo
             !
          end do !idir
          !
       end do !ifout
       !
    end do !ip
   ! !$OMP END DO
   ! !$OMP END PARALLEL!
    
    !Ebnd( ndir , nfbnd,ntbnd )
   ! call log_real( 'E',E( 1,1, nfbnd,5),0 )
    !call log_real( 'Ebnd',Ebnd( ndir , nfbnd,ntbnd ),0 )
  end subroutine getspec
  
end module raytracing
!
