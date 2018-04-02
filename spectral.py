import code
import numpy as np
#==========================================================================================
# This file contains the definitions of the different spectral classes:
# 
# tSpectrum  :: base class
# spectrum   :: 2d spectrum class
# spectrum1d :: 1d spectrum class
#
# In addition it contains some helper routines
#
#
#
#==========================================================================================
class tSpectrum:
    #
    #
    def __init__( self , var ):
        #
        import scipy
        from scipy import io
        #
        if type(var) == str:
            #
            # Load spectrum from a file
            #
            var = scipy.io.loadmat( var )
        #    
        #
        
        self.E   = np.squeeze(np.array(var['E']))
        if 'f' in var:
            #
            self.freqaxis = 'f'
            self.f   = np.squeeze(np.array(var['f']))
            self.f2w     = 2. * np.pi
            self.f2f     = 1.
            #
        elif 'w' in var:
            #
            self.freqaxis = 'w'                
            self.f   = np.squeeze(np.array(var['w']))
            self.f2w     = 1
            self.f2f     = 1/ (2. * np.pi)
            #
        else:
            #
            self.freqaxis = None
            self.f      = []            

        
        if ('ang' in var) or ('deg' in var):
            #
            self.diraxis = 'ang'     
            self.ang  = np.squeeze(np.array(var['ang']))
            self.dir2rad = np.pi/180.
            self.dir2deg = 1.
            #
        elif 'rad' in var:
            #
            self.diraxis = 'rad'                
            self.ang   = np.squeeze(np.array(var['rad']))
            self.dir2rad = 1.
            self.dir2deg = 180./np.pi
            #
        else:
            #
            self.diraxis = None
            self.ang     = []
        #

        maxdim = 3
        self.type = '2d'
        if self.diraxis == None or self.freqaxis==None:
            #
            self.type = '1d'
            maxdim = 2
            #
        
            
        ndim = len(np.shape(self.E))
        if ndim == maxdim:
            #
            self.nloc= np.shape(self.E)[0]
            #
        else:
            #
            self.nloc = 1
            #         
            
        if 'loc' in var:
            #
            self.loc = np.array(var['loc'])
            #
        else:
            #
            self.loc = list( range(0,self.nloc) )
            #
            
        self.nf  = len(self.f)
        self.ndir= len(self.ang)
        #
    # end init

    def save( self , filename ):
        #
        # Save a set of spectral observations
        #
        import scipy
        from scipy import io

        if self.diraxis==None:
            #
            dic = { 'E':self.E , self.freqaxis:self.f, 'loc':self.loc,'a1':self.a1,'a2':self.a2,'b1':self.b1,'b2':self.b2 }
            #
        elif self.freqaxis==None:
            #
            dic = { 'E':self.E , self.diraxis:self.ang, 'loc':self.loc }
            #
        else:
            #
            dic = { 'E':self.E , self.freqaxis:self.f, self.diraxis:self.ang, 'loc':self.loc }
            #
        #
        scipy.io.savemat( filename , dic )
        
    def todeg(self):
        #
        # Change directional units of the spectrum to degrees
        # - if already in degrees, nothing happens
        self.E   = self.E / self.dir2deg
        self.ang = self.ang * self.dir2deg
        self.diraxis = 'ang'
        self.dir2rad = np.pi/180.
        self.dir2deg = 1.
        #


    def torad(self):
        #
        # Change directional units of the spectrum to radians
        # - if already in degrees, nothing happens
        self.E   = self.E   / self.dir2rad
        self.ang = self.ang * self.dir2rad
        self.diraxis = 'rad'
        self.dir2rad = 1.
        self.dir2deg = 180/np.pi
        #
        
    def tow(self):
        #
        # Change frequency units of the spectrum to angular frequency
        # - if already in angular frequency, nothing happens
        self.E   = self.E / self.f2w
        self.f   = self.f * self.f2w
        self.freqaxis = 'w'
        self.f2w = 1.
        self.f2f = 1/ (2. * np.pi)
        #

    def tof(self):
        #
        # Change frequency units of the spectrum to Hz
        # - if already in Hz, nothing happens
        self.E   = self.E / self.f2f
        self.f   = self.f * self.f2f
        self.freqaxis = 'f'
        self.f2w = 2. * np.pi
        self.f2f = 1.
        #  
        
    def angrad(self):
        #
        return( self.ang * self.dir2rad )

    def angdeg(self):
        #
        return( self.ang * self.dir2deg )

    def fw(self):
        #
        return( self.f * self.f2w )

    def ff(self):
        #
        return( self.f * self.f2f )  
        
    def interpLoc(self, coordinate, method='linear'):
        #
        return(interpolateSpec( self , coordinate, method ))
        #

    def interpFreq(self,f):
        #
        import scipy
        from scipy import interpolate
        import numpy

        if self.type=='2d':
            interpolant = scipy.interpolate.interp1d( self.f , self.E, axis=-1,bounds_error=False,assume_sorted=True )
            self.E  = interpolant( f )
            self.E[numpy.isnan(self.E) ]=0.
            self.f  = f
            self.nf = len(f)
        else:
            interpolant = scipy.interpolate.interp1d( self.f , self.E, axis=-1,bounds_error=False,assume_sorted=True )
            self.E  = interpolant( f )
            self.E[numpy.isnan(self.E) ]=0.

            interpolant = scipy.interpolate.interp1d( self.f , self.a1, axis=-1,bounds_error=False,assume_sorted=True )
            self.a1  = interpolant( f )
            self.a1[numpy.isnan(self.a1) ]=0.
            
            interpolant = scipy.interpolate.interp1d( self.f , self.a2, axis=-1,bounds_error=False,assume_sorted=True )            
            self.a2  = interpolant( f )
            self.a2[numpy.isnan(self.a2) ]=0.
            
            interpolant = scipy.interpolate.interp1d( self.f , self.b1, axis=-1,bounds_error=False,assume_sorted=True )            
            self.b1  = interpolant( f )
            self.b1[numpy.isnan(self.b1) ]=0.
            
            interpolant = scipy.interpolate.interp1d( self.f , self.b2, axis=-1,bounds_error=False,assume_sorted=True )            
            self.b2  = interpolant( f )
            self.b2[numpy.isnan(self.b2) ]=0.            
            
            self.f  = f
            self.nf = len(f)            

    def Tm0n(self,n=1):
        #
        #
        # Calculate Tm0n, Tm0n, etc.
        # 
        import numpy
        #
        if n==1:
            #
            m0 = self.m(0)
            m1 = self.m(1)
            return( m0 / m1 )
        #
        elif n==2:
            #
            m0 = self.m(0)
            m2 = self.m(2)
            return( numpy.sqrt( m0 / m2) )
            #
        #endif
        #
        
    def Hm0(self):
        #
        import numpy
        #
        return( 4.*numpy.sqrt(self.m(0)) )
        #

    def Tm01(self):
        #
        return( self.Tm0n(1) )
        #    

    def Tm02(self):
        #
        return( self.Tm0n(2) )
        #           
        
#END tSpectrum

#==========================================================================================
class spectrum(tSpectrum):
    #
    def __init__( self , var ):
        #
        tSpectrum.__init__( self , var )      
        #
    #end def 

    def correctDim(self):
        #
        # Add leading location dimension for calculations...
        #
        E = self.E
        #
        if self.nloc == 1:
            #
            E = E[None,:,:]
            #
        #
        return(E)
        #
    #end correctDim

    def constraint( self, ang0, length ):
        # set zero part of spectrum
        self.center

        ang = self.ang - ang0
        ang[ ang < -np.pi / self.dir2rad ] =  ang[ ang < -np.pi / self.dir2rad ] + 2 * np.pi /self.dir2rad

        msk = np.logical_and( ang >= 0 , ang < length )
        
        E= self.correctDim()
        E[ : , ~msk , : ] = 0.
        self.E = np.squeeze( E )
        #
        
    def regularize( self, n ):
        #
        # 
        #
        import scipy
        from scipy import interpolate
        
        self.center
        
        ang = np.linspace( -np.pi/self.dir2rad,  np.pi/self.dir2rad,n , endpoint=False )
        E = self.E
        interpolant = scipy.interpolate.interp1d( self.ang , self.E,
                                                      axis=-2,bounds_error=False,assume_sorted=True,fill_value=0. )
        self.E  = interpolant( ang )
        self.E[np.isnan(self.E) ]=0.
        self.ang  = ang
        self.ndir = len(ang)
        #
        
    def center(self):
        #
        # Centre the spectrum to -180 (pi) < angle < 180 (pi)
        #
        self.ang = np.arctan2( np.sin( self.ang * self.dir2rad ) , np.cos( self.ang * self.dir2rad ) ) / self.dir2rad
        indices = np.argsort( self.ang )
        self.ang = self.ang[indices]
        self.E   = np.squeeze(self.correctDim()[:,indices,:])
        #
    #end def

    def bulkSprd(self):
        #
        # Calculate the bulk directional spreading according to Kuik (see also the description in the SWAN 
        # manual, http://swanmodel.sourceforge.net/online_doc/swanuse/node35.html, under DSPR)

        #Convert to radians (if needed)        
        ang = self.ang * self.dir2rad
        E   = self.correctDim() / self.dir2rad

        #calculate cos/sin factors as arrays of size E
        sin1 = np.sin(      ang )
        cos1 = np.cos(      ang )
        sin1 = np.tile( sin1[None,:,None] , [self.nloc,1,self.nf] )
        cos1 = np.tile( cos1[None,:,None] , [self.nloc,1,self.nf] )

        #Directional integration factors
        dang = ang[2] - ang[1]
        df = deltaf(self.f)
        df = np.tile( df , [self.nloc,self.ndir,1] )

        #calculate frequency weighted Lon. Hig. bulk directional moments using midpoint rule for
        #numerical integration
        a1 = dang * np.sum(  np.sum( E * cos1 * df,  axis=-1) ,  axis=-1 )
        b1 = dang * np.sum(  np.sum( E * sin1 * df,  axis=-1) ,  axis=-1 )
        m0 = self.m(0)        
        m1  = np.sqrt( (a1/m0)**2 + (b1/m0)**2 )

        # return the directional with (Kuik, eq. 34), and convert back to directional units of spectrum
        return( np.sqrt( 2. * ( 1. - m1 ) ) / self.dir2rad )
        #
    #end def
        
    def bulkDir(self):
        #
        # Calculate the bulk direction of the spectrum
        #
        import numpy

        #Convert to radians (if needed)
        ang = self.ang * self.dir2rad
        E   = self.correctDim() / self.dir2rad

        #calculate cos/sin factors as arrays of size E
        sin1 = numpy.sin(      ang )
        cos1 = numpy.cos(      ang )
        sin1 = numpy.tile( sin1[None,:,None] , [self.nloc,1,self.nf] )
        cos1 = numpy.tile( cos1[None,:,None] , [self.nloc,1,self.nf] )

        #Directional integration factor
        dang = ang[2] - ang[1]

        #Calculate a1/b1 by integration using midpoint rule for directions,
        # and trapezoidal rule for frequencies
        a1 = dang * numpy.sum(  numpy.trapz( E * cos1, x=self.f, axis=-1) ,  axis=-1 )
        b1 = dang * numpy.sum(  numpy.trapz( E * sin1 ,x=self.f, axis=-1) ,  axis=-1 )  
        
        return( numpy.arctan2( b1 , a1 ) * 1./self.dir2rad )
        #
    #end bulkDir
        
    def meanDir(self):
        #
        # Calculate the mean direction per frequency
        # of the spectrum from the directional moments a1,b1
        # 
        import numpy

        #Convert to radians (if needed)
        ang = self.ang * self.dir2rad
        E   = self.correctDim() / self.dir2rad

        #calculate cos/sin factors
        sin1 = numpy.sin(      ang )
        cos1 = numpy.cos(      ang )
        sin1 = numpy.tile( sin1[None,:,None] , [self.nloc,1,self.nf] )
        cos1 = numpy.tile( cos1[None,:,None] , [self.nloc,1,self.nf] )

        #Directional integration factor
        dang = ang[2] - ang[1]

        #Calculate a1/b1 by integration using midpoint rule
        a1 = dang * numpy.sum(  E * cos1 ,  axis=-2 )
        b1 = dang * numpy.sum(  E * sin1 ,  axis=-2 )
        
        return( numpy.arctan2( b1 , a1 )*1./self.dir2rad )
        #
    #end meanDir

    def Ef( self , kind=None):
        #
        # return a 1d spectral object
        #
        import numpy
        import copy
        
        dang = self.ang[2] - self.ang[1]        
        Ef = dang * numpy.sum( self.E , axis=-2 )        
        Ef = numpy.squeeze(Ef)

        a1,b1,a2,b2 = self.directionalMoments()

        if kind==None:
            #
            # Return 1d spectrum with dimensions in place
            #
            if self.f2f == 1.:
                return( spectrum1d( {'E':Ef,'f':copy.deepcopy(self.f),'loc':copy.deepcopy(self.loc),
                                 'a1':a1,'b1':b1,'a2':a2,'b2':b2} ) )
            else:
                return( spectrum1d( {'E':Ef,'w':copy.deepcopy(self.f),'loc':copy.deepcopy(self.loc),
                                 'a1':a1,'b1':b1,'a2':a2,'b2':b2} ) )
        elif (kind=='f'):
            # return frequency spectrum
            return( spectrum1d( {'E':Ef/self.f2f,'f':copy.deepcopy(self.ff),'loc':copy.deepcopy(self.loc),
                                 'a1':a1,'b1':b1,'a2':a2,'b2':b2} ) )
        elif (kind=='w'):
            # return frequency spectrum
            return( spectrum1d( {'E':Ef/self.f2w,'w':copy.deepcopy(self.fw),'loc':copy.deepcopy(self.loc),
                                 'a1':a1,'b1':b1,'a2':a2,'b2':b2} ) )        
        #

    def Ed( self ):
        #
        # return a 1d spectral object
        #
        import numpy
        import copy

        df = deltaf(self.f)
        df = numpy.tile( df , [self.nloc,self.ndir,1] )
        E = self.correctDim()

        Ed = numpy.sum( self.E * df , axis=-1 )
        Ed = numpy.squeeze(Ed)
        #
        return( spectrum1d( {'E':Ed,self.diraxis:copy.deepcopy(self.ang),'loc':copy.deepcopy(self.loc)} ) )
        #
    #def
    #

    #
    def plot(self  , iloc=0, fig=None,subplot=None):
        #
        # A simple function to visually inspect the spectrum
        #
        import matplotlib.pyplot as plt
        import numpy as np

        # make sure it as a leading iloc dimension
        E  = np.squeeze(self.correctDim()[iloc,:,:])
        Ef = np.squeeze(self.Ef().correctDim()[iloc,:])
        Ed = np.squeeze(self.Ed().correctDim()[iloc,:])


        ii = np.argmax( self.Ef().correctDim()[iloc,:]  )        
        # Plot directional spectrum

        if fig == None:
            #
            fig = plt.figure()            
            #

        if subplot==None:
            #
            subplot = [1,3,1,2,3]
            #
        
        ax1=fig.add_subplot( subplot[0],subplot[1],subplot[2])
        ax1.pcolor( self.f, self.ang, E )

        # Plot frequency spectrum
        ax2=fig.add_subplot(subplot[0],subplot[1],subplot[3])
        ax2.plot( self.f, Ef )

        # Plot direction spectrum
        ax3=fig.add_subplot(subplot[0],subplot[1],subplot[4])
        ax3.plot( self.ang, E[:,ii] )        
        fig.canvas.draw()
        fig.show()
        return(fig)
        
    def directionalMoments( self ):
        #
        # Calculate the directionalMoments
        #
        import numpy
        #
        #
        #convert angles to radians (if appropriate)
        ang = self.ang * self.dir2rad
        E   = self.correctDim() / self.dir2rad

        #get frequency spectrum for normalization...
        #Ef  = self.Ef().E
        dang = ang[2] - ang[1]        
        Ef = dang * numpy.sum( E , axis=-2 )
        #print(np.max(Ef))
        #... but first ensure that no div. by zero occurs ...
        msk = Ef <= 0.
        Ef[msk] = 1.
        Ef  = numpy.tile( Ef[:,None,:] , [1,self.ndir,1] )
        
        # ... and then normalize directional distributions.
        En  = E / Ef

        #calculate the sin/cos factors, and expand ndir to ndir/nfreq for
        #pointwise multiplication
        sin1 = numpy.sin(      ang )
        cos1 = numpy.cos(      ang )
        sin2 = numpy.sin( 2. * ang )
        cos2 = numpy.cos( 2. * ang )        

        sin1 = numpy.tile( sin1[None,:,None] , [self.nloc,1,self.nf] )
        cos1 = numpy.tile( cos1[None,:,None] , [self.nloc,1,self.nf] )
        sin2 = numpy.tile( sin2[None,:,None] , [self.nloc,1,self.nf] )
        cos2 = numpy.tile( cos2[None,:,None] , [self.nloc,1,self.nf] )

        #
        # Calculate directional Fourier moments using midpoint rule
        #
        dang = ang[2] - ang[1]
        a1 = dang * numpy.sum(  En * cos1 ,  axis=-2 )
        b1 = dang * numpy.sum(  En * sin1 ,  axis=-2 )
        a2 = dang * numpy.sum(  En * cos2 ,  axis=-2 )
        b2 = dang * numpy.sum(  En * sin2 ,  axis=-2 )
        #
        # Return coef.
        #
        return( a1 , b1 , a2 , b2 )
        
    def m( self , n = 0):
        #
        # Calculate the spectral moments
        #
        import numpy
        #
        fn = self.ff()**n
        df = deltaf(self.f)        
        fn = numpy.tile( fn , [self.nloc,self.ndir,1] )
        df = numpy.tile( df , [self.nloc,self.ndir,1] )        
        E  = self.correctDim() * fn
        
        dang = self.ang[2] - self.ang[1]
        m = dang * numpy.sum( numpy.sum( E * df, axis=-1) , axis=-1 )
 
        if len(m) == 1:
            m = m[0]
        return( m )
        #
        
    #

#==========================================================================================
class spectrum1d(tSpectrum):
    #
    def __init__( self , var=None ):
        #
        import numpy

        tSpectrum.__init__( self , var )      
        #

        if 'a1' in var:
            self.a1 = numpy.array(var['a1'])
        else:
            self.a1 = np.nan + self.E
            
        if 'a2' in var:
            self.a2 = numpy.array(var['a2'])
        else:
            self.a2 = np.nan + self.E
            
        if 'b1' in var:
            self.b1 = numpy.array(var['b1'])
        else:
            self.b1 = np.nan + self.E
            
        if 'b2' in var:
            self.b2 = numpy.array(var['b2'])
        else:
            self.b2 = np.nan + self.E            
        #
    #end init
    
    def correctDim(self):
        #
        # Add leading location dimension for calculations...
        #
        E = self.E
        if self.nloc == 1:
            E = E[None,:]
        return(E)     

    def bulkDir(self):
        #
        # Calculate the bulk direction of the spectrum from the directional moments according to KUIK
        #
        import numpy
        df = deltaf(self.f)
        df = numpy.tile( df , [self.nloc,1] )

        m  = self.m()
        a1 = numpy.sum( self.E * self.a1 * df , axis=-1 ) /m
        b1 = numpy.sum( self.E * self.b1 * df , axis=-1 ) /m
        
        return(numpy.arctan2( b1,a1) )

    def bulkSprd(self):
        #
        # Calculate the bulk spread of the spectrum from the directional moments according to KUIK
        #
        import numpy
        
        df = deltaf(self.f)
        df = numpy.tile( df , [self.nloc,1] )

        m  = self.m()
        E = self.correctDim()


        a1 = numpy.sum( self.E * self.a1 * df , axis=-1 ) /m
        b1 = numpy.sum( self.E * self.b1 * df , axis=-1 ) /m

        return(   numpy.sqrt(  2*( 1-numpy.sqrt(a1**2+b1**2) )  )   )

    def dirSprd(self):
        #
        # Calculate the bulk spread of the spectrum from the directional moments according to KUIK
        #
        import numpy
        

        m  = self.m()
        a1 = self.a1
        b1 = self.b1

        
        return(   numpy.sqrt(  2*( 1-numpy.sqrt(a1**2+b1**2) )  )   )

    def meanDir(self):
        #
        # Calculate the bulk spread of the spectrum from the directional moments according to KUIK
        #
        import numpy
        

        m  = self.m()
        a1 = self.a1
        b1 = self.b1

        
        return(   numpy.arctan2( b1 , a1 )   )      

    
    def plot(self  , iloc=0, fig=None):
        #
        # A simple function to visually inspect the spectrum
        #
        import matplotlib.pyplot as plt
        import numpy as np

        # make sure it as a leading iloc dimension
        Ef = np.squeeze(self.correctDim()[iloc,:])
        # Plot directional spectrum
        plt.figure    
        plt.plot( self.f, Ef )        
        plt.show(block=False)
        return()
    
    def Tp( self ):
        import numpy

        ii = numpy.argmax( self.E,axis=-1 )


        try:
            res = numpy.zeros((len(ii)))
        except:
            ii = (ii,)
            res = numpy.zeros((len(ii)))
            
        for jj in range( 0,len(ii)):
            #
            res[jj] = 1./self.f[ii[jj]] #self.E[jj,ii[jj]]
            #
        return(res)           

    def peakDir( self ):
        import numpy

        ii = numpy.argmax( self.E,axis=-1 )
        if self.nloc == 1:
            #
            res = numpy.arctan2( self.b1[ii],self.a1[ii] )
            #
        else:
            #
            res = numpy.zeros((len(ii)))

            for jj in range( 0,len(ii)):
                #
                res[jj] = numpy.arctan2( self.b1[jj,ii[jj]],self.a1[jj,ii[jj]] )
                #
            #endfor
            #
        #endif
        #
        return(res)
        #
    #enddef

    def peakSprd( self ):
        import numpy

        ii = numpy.argmax( self.E,axis=-1 )
        if self.nloc == 1:
            #
            res = numpy.sqrt(  2*( 1-numpy.sqrt(self.a1[ii]**2+self.b1[ii]**2) )  )
            #
        else:
            res = numpy.zeros((len(ii)))
            
            for jj in range( 0,len(ii)):
                res[jj] = numpy.sqrt(  2*( 1-numpy.sqrt(self.a1[jj,ii[jj]]**2+self.b1[jj,ii[jj]]**2) )  )
                #
            #endfor
            #
        return(res)
        #
    #enddef
    
    def m( self , n = 0):
        #
        import numpy
        #

        fn = self.ff()**n
        df = deltaf(self.f)

        fn = numpy.tile( fn , [self.nloc,1] )
        df = numpy.tile( df , [self.nloc,1] )
        E  = self.correctDim() * fn
        E[ numpy.isnan(E) ] = 0.
        m = numpy.sum( E * df  , axis=-1)
        #
        if len(m) == 1:
            m = m[0]
        return( m )
        #
    #enddef
    #
    def salt( self , uncertainty , seed=None ):
        #
        # this definition "salts" the spectrum with random noise to mimick observational
        # uncertainty
        #
        import copy
        np.random.seed(seed)
        std = self.E * uncertainty / 2. 

        spec = copy.deepcopy( self )
        #
        for iloc in range( 0 , self.nloc ):
            #
            for ifreq in range( 0, self.nf ):
                #
                std    = self.E[iloc,ifreq] * uncertainty / 2.
                errors = np.random.normal( loc = 0., scale = std  )
                spec.E[iloc,ifreq] = spec.E[iloc,ifreq] + errors

                if spec.E[iloc,ifreq] < 0.:
                    spec.E[iloc,ifreq] = 0.
                    
                r1 = np.sqrt(spec.a1[iloc,ifreq]**2 +  spec.a1[iloc,ifreq]**2)
                r2 = np.sqrt(spec.a2[iloc,ifreq]**2 +  spec.a2[iloc,ifreq]**2)                    

                std    = np.abs( r1 * uncertainty / 2.)
                errors = np.random.normal( loc = 0., scale = std  )
                spec.a1[iloc,ifreq] = spec.a1[iloc,ifreq] + errors

                std    = np.abs(r1 * uncertainty / 2.)                        
                errors = np.random.normal( loc = 0., scale = std  )
                spec.b1[iloc,ifreq] = spec.b1[iloc,ifreq] + errors
                
                std    = np.abs(r2 * uncertainty / 2.)
                errors = np.random.normal( loc = 0., scale = std  )
                spec.a2[iloc,ifreq] = spec.a2[iloc,ifreq] + errors

                std    = np.abs( r2 * uncertainty / 2.)
                errors = np.random.normal( loc = 0., scale = std  )
                spec.b2[iloc,ifreq] = spec.b2[iloc,ifreq] + errors                

                r1 = np.sqrt(spec.a1[iloc,ifreq]**2 +  spec.a1[iloc,ifreq]**2)
                r2 = np.sqrt(spec.a2[iloc,ifreq]**2 +  spec.a2[iloc,ifreq]**2)

                if r1 > 1.:
                    spec.a1[iloc,ifreq] = spec.a1[iloc,ifreq]/r1/1.01
                    spec.b1[iloc,ifreq] = spec.b1[iloc,ifreq]/r1/1.01
                if r2 > 1.:
                    spec.a2[iloc,ifreq] = spec.a2[iloc,ifreq]/r2/1.01
                    spec.b2[iloc,ifreq] = spec.b2[iloc,ifreq]/r2/1.01
            #end for ifreq
            #
        #end for iloc
        #
        return( spec )
    #enddef
#end class spectrum1d
#
#==========================================================================================
def deltaf( f ):
    #
    import numpy
    delta = numpy.diff( f )
    forward  = numpy.append( delta,delta[-1] )
    backward = numpy.insert( delta , 0,  delta[0] )
    return( (forward + backward) * 0.5 )


def interpolateSpec( A , c0, method='linear' ):
    #
    # routine to interpolate spectrum
    #
    import numpy as np
    import copy

    c = A.loc
    it = (np.abs(c-c0)).argmin()
    if c[it] - c0 < 0:
        it2 = np.min([it + 1,len(c)-1])
        it1 = it        
    else:
        it1 = np.max([it - 1,0])
        it2 = it
    
        
    
    dist = c[it2] - c[it1]
    #
    if it1 == it2:
        #
        dist = 1
        fac  = ( 1, 0)
    else:
        #
        fac  = ( c0/dist - c[it1]/dist , c[it2]/dist - c0/dist)

    if type(A)==spectrum:
        #
        # Create new spectral object with same angles/directions as A/B
        C = spectrum( {'E':np.zeros((A.ndir,A.nf)), 'f':copy.deepcopy(A.f),
                'ang':copy.deepcopy(A.ang)} )
        #
        if method.lower() == 'linear':
            #
            C.E = fac[1] * A.E[it1,:,:] + fac[0] * A.E[it2,:,:] 
            #
        #endif
        #
    elif type(A)==spectrum1d:
        #
        #
        # Create new spectral object with same angles/directions as A/B
        C = spectrum1d( {'E':np.zeros((A.nf)), 'f':copy.deepcopy(A.f)})
        #
        if method.lower() == 'linear':
            #
            C.E = fac[1] * A.E[it1,:] + fac[0] * A.E[it2,:]
            C.a1 = fac[1] * A.a1[it1,:] + fac[0] * A.a1[it2,:]
            C.a2 = fac[1] * A.a2[it1,:] + fac[0] * A.a2[it2,:]
            C.b1 = fac[1] * A.b1[it1,:] + fac[0] * A.b1[it2,:]
            C.b2 = fac[1] * A.b2[it1,:] + fac[0] * A.b2[it2,:] 
            #
        #endif
        #
    #endif
    #
    return( C )
    #
#enddef

def load( filename ):
    #
    # Load a spectral object from disk and return a spectrum object
    #
    import scipy
    from scipy import io

    var = scipy.io.loadmat( filename )

    isfreq = ('f' in var or 'w' in var)
    isdir = ('ang' in var or 'rad' in var)
    #
    if isfreq and isdir:
        #
        # 2d Spectrum
        #
        return( spectrum( var ) )
        #
    else:
        #
        # 1d Spectrum
        #
        return( spectrum1d( var ) )
        #
    #endif
    #
#enddef
#

def parSpec( freq  ,Hs=1., Tp=10., ang=None, freqshape='Jonswap' , dirshape='cosn',
                 cosPower=20. ,meanAng =0., gamma=3.3 , sb = 0.09 , sa=0.07 , g = 9.81 , dirWidth = 10 * np.pi / 180. ):

    from scipy.stats import norm    
    #
    # Define a parametric spectrum based on frequency
    #
    if freqshape.lower() == 'jonswap':
        #
        Ef = parJonswap( freq ,1./Tp, Hs,gamma,sb,sa,g )
        #
    if freqshape.lower() == 'white':
        #
        # White spectrum with given
        #
        Ef = np.ones( len(freq) ) * ( Hs/4. ) **2 / np.sum( freq )
        #
    #endif
    #
    if ang is not None:
        #
        if dirshape.lower() == 'cosn':
            #
            D = parcosn( ang , meanAng=meanAng , cosPower=cosPower )
            #
        elif dirshape.lower() == 'guassian':
            #
            D = norm.pdf( angdif( ang , meanAng ) , loc = 0. , scale = dirWidth )
            #
        #endif
        #
        E = D[:,None] @ Ef[None,:]
        E = spectrum( { 'E':E , 'f':freq , 'rad':ang } )
            
        #
    else:
        #
        E = spectrum1d( { 'E':Ef , 'f':freq } )
        #
    #end if
    #
    return( E )
    #
#enddef parSpec

def parJonswap( freq , peakfreq, Hs, gamma=3.3 , sb = 0.09 , sa=0.07 , g = 9.81 ):
    #
    # JONSWAP spectral shape (See, e.g. Holthuijsen waves in coastal and oceanic water)
    #
    #tail
    Et = g**2 / ( ( np.pi * 2 )**4 ) * freq**(-5)

    #cut-off
    Ec = np.exp( - 5./4. * (freq/peakfreq)**(-4))

    #peak enhancement
    sigma = np.zeros( len(freq) )
    sigma[ freq <= peakfreq ] = sa
    sigma[ freq >  peakfreq ] = sb
    Ep = gamma**(    np.exp(   -.5 * (  ( (freq/peakfreq - 1.)/sigma )**2  )   )    )

    #raw spec
    E = Et * Ec * Ep
    E[ np.isnan(E) ] = 0.

    #Scaling
    E = E * ( Hs/4. ) **2 / np.trapz( freq , E ) 
    return( E )
    #
#enddef parJonswanp

def angdif( ang1 , ang2 ):
    #
    # Return the smallest mutual angle between two angles
    # i.e. ang1 - ang2
    ang1 = np.arctan2( np.sin( ang1 ) , np.cos( ang1 ) )
    ang2 = np.arctan2( np.sin( ang2 ) , np.cos( ang2 ) )

    diff = ang1 - ang2
    msk = diff>np.pi
    diff[msk] = -2*np.pi + diff[msk]
    msk = diff<-np.pi
    diff[msk] = 2*np.pi + diff[msk]
    #
    return( diff )
    #
#enddef angdif
        
def parcosn( ang , meanAng=0. , cosPower=2. ):
    #
    # Raised cosine (cos^n distribution) directional distribution (See, e.g. Holthuijsen waves
    # in coastal and oceanic water)
    #
    ang = angdif( ang , meanAng )
    #
    D = np.zeros( len( ang ) )
    msk = ( ang < np.pi/2 ) & ( ang > -np.pi/2 )
    D[msk] = np.cos( ang[msk] )**cosPower

    D = D / np.trapz( ang , D )
    return(D)

#enddef parcosn
    
