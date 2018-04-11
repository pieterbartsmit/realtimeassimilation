import numpy as np
import ctypes
import swan

class raytracer:
    #
    def __init__( self , cgrid , lon , lat , dep, bndlist , dcont, workdir='./',save=True,radangles=None,angfreq=None ):
        #
        import platform
        
        #
        # A rudimentary check to see if I am on my local mac or on the remote server
        # mostly convinient for development; should be removed later
        #
        if platform.system() == 'Darwin':
            #
            libloc = '/users/pieterbartsmit/Google Drive/repos/DAI/fortran/raytracing.lib'
            #
        else:
            #
            libloc = '/mnt/2ndHDD/fortran/raytracing.lib';
            #        
        #endif

        #load the external library
        self.lib = ctypes.cdll.LoadLibrary( libloc )

        #Set the domain        
        self.cgrid   = cgrid
        doublearray2 = ctypes.c_double * 2        
        xlim = doublearray2( self.cgrid.domain[0] , \
                                     self.cgrid.domain[0]+ self.cgrid.domain[2] )
        ylim = doublearray2( self.cgrid.domain[1] , \
                                     self.cgrid.domain[1]+ self.cgrid.domain[3] )
        #
        self.lib.setdom.restype = None
        self.lib.setdom( ctypes.byref( xlim) , ctypes.byref( ylim ) )

        
        #Set/Load Bathymetry in library
        self.lon  = lon
        self.lat  = lat
        self.dep  = dep
        self.nlon = len(lon)
        self.nlat = len(lat)
        self.n    = self.nlon * self.nlat
        self.workdir =  workdir
        self.save    = save

        lonarr = ctypes.c_double * self.nlon
        latarr = ctypes.c_double * self.nlat
        deparr = ctypes.c_double * self.n

        d=np.reshape(self.dep,(self.n,) )

        depi = (ctypes.c_double * self.n   )(*d.tolist())
        loni = (ctypes.c_double * self.nlon)(*lon.tolist())
        lati = (ctypes.c_double * self.nlat)(*lat.tolist())
        nloni = ctypes.c_int(self.nlon)
        nlati = ctypes.c_int(self.nlat)        
        
        self.lib.setbat.restype = None
        self.lib.setbat( ctypes.byref( depi  ) ,
                         ctypes.byref( loni  ) ,
                         ctypes.byref( lati  ) ,
                         ctypes.byref( nloni ) ,
                        ctypes.byref( nlati ) )                       
        #
        self.bndlist = []
        for bnd in bndlist:
            #
            if bnd in ['W','w',1]:
                #
                self.bndlist.append( 1 )
                #
            elif bnd in ['E','e',2]:
                #
                self.bndlist.append( 2 )
                #
            elif bnd in ['S','s',3]:
                #
                self.bndlist.append( 3 )
                #
            elif bnd in ['N','n',4]:
                #
                self.bndlist.append( 4 )                
                #
            #
        #
        bndlisti = (ctypes.c_int * len(self.bndlist) )(*self.bndlist)
        nbndlisti = ctypes.c_int( len(self.bndlist) )
        self.lib.setbndlist.restype = None        
        self.lib.setbndlist( ctypes.byref(bndlisti) , ctypes.byref(nbndlisti) )

        self.maxstep = 10000
        self.frac    = 1
        self.lib.setnum.restype = None
        self.lib.setnum( ctypes.byref( ctypes.c_double( self.frac    ) ),
                         ctypes.byref( ctypes.c_int(    self.maxstep ) ) )

        self.lib.setcont.restype = None
        self.lib.setcont( ctypes.byref(ctypes.c_double( dcont ) ) )
        self.radangles = radangles
        self.angfreq   = angfreq
        self.nsub = 20 #This seems to be unused - no time to check, but be aware!
    #end __init__

    def getmatrix( self , xp , yp , radangles,nsub, angfreq ):
        #
        # Wrapper to calculate the dependency matrix using the raytracer code
        # NOTE: angles and frequencies are LISTS in terms of RADIANS and
        #       angular frequency!
        #
        # Is there already a workfile...?
        #
        import os
        import spectral

        self.nsub = nsub
        
        workfile =  self.workdir + 'matrix.npy'
        fileexists = os.path.isfile(workfile)
        #
        if fileexists and self.save:
            #
            # ...if so load the transfermatrix...
            #
            matrix = np.load(workfile)
            #
        else:
            #
            # ... else we calculate the transfermatrix; for each input
            # transform to the appropriate ctype
            #
            nang   = len(radangles)
            angles = (ctypes.c_double * nang)(*radangles)
            matrix = np.zeros( nang * nang * len(angfreq) * len(xp) )
            matrix = (ctypes.c_double * len( matrix ) )(*matrix)
            xp     = (ctypes.c_double * len( xp ) )(*xp)
            yp     = (ctypes.c_double * len( yp ) )(*yp)
            npo    = ctypes.c_int( len(xp) )
            nfreq = ctypes.c_int( len(angfreq) )
            nangi = ctypes.c_int( nang )
            nsubi = ctypes.c_int( nsub )        
            freq  = (ctypes.c_double * len( angfreq ) )(*angfreq)
            #
            # ... Call the library ...
            #
            self.lib.calc.restypes = None
            self.lib.calc( ctypes.byref( matrix ),
                               ctypes.byref( xp ),
                               ctypes.byref( yp ),
                               ctypes.byref( npo),
                               ctypes.byref( freq ),                       
                               ctypes.byref( nfreq ),
                               ctypes.byref( angles ),
                               ctypes.byref( nangi ),
                               ctypes.byref( nsubi ) )

            #
            # ...Reshape to correct dimensions and return transfer matrix...
            #
            matrix = np.array( matrix ).reshape( ( len(angfreq),len(xp),nang,nang) )
            #
            #
            # ... and save it to a file.
            #
            np.save(workfile,matrix)
            #
        #endif
        #
        return(matrix)        
        #
    #endDef

    def prediction( self , dirspec , xp , yp , nsub  ):
        #
        # Predict spectrum at location based on boundary spectrum
        #        
        import os
        import spectral

        dirconv = False
        freqconv = False

        # ensure that input is in ang. freq and rad
        if dirspec.freqaxis == 'f':
            freqconv = True
            dirspec.tow()

        if dirspec.diraxis in ['ang','deg']:
            dirconv = True
            dirspec.torad()            

        self.radangles = dirspec.angrad()
        self.angfreq   = dirspec.fw()
        matrix = self.getmatrix(  xp , yp , dirspec.angrad().tolist() , nsub, dirspec.fw().tolist() )        
      
        #
        # Calculate the local spectrum
        #
        E0 = np.zeros( (len(xp), dirspec.nf, dirspec.ndir) )
        #
        for jf in range(0, dirspec.nf):
            #
            # For each frequency
            #
            for ip in range( 0 , len(xp) ):
                #
                # and at each point, multply dependency matrix with offshore spectrum
                # to obtain local spectrum
                #
                E0[ip,jf,:]  = np.matmul( matrix[ jf,ip , :, : ] , dirspec.correctDim()[ 0,: , jf ] )
                #
            #end for
            #
        #end for
        #
        # Finally, reorder output...        
        E0 = np.transpose( E0 , [0,2,1] )

        # and return a spectrum object
        spec = spectral.spectrum( { 'E':E0 , 'w':dirspec.f , 'rad':dirspec.ang } )

        # Ensure that a spectrum with the same units is returned
        if freqconv:
            dirspec.tof()
            spec.tof()
            
        if dirconv:
            dirspec.todeg()
            spec.todeg()
            
        return( spec )
        #
    #end prediction
#
