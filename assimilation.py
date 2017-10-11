import numpy as np

class assimilation:
    #
    # The assimilation object that performs the data assimilation...
    #
    def __init__( self , raytracerobject , xp , yp, halfCircleStart=None,
                      momentsPerSpotter = 5,kind='energy', nyquistPenalty=0 ):
        #
        # INPUT:
        #
        # raytracerobject  :: raytracerobject instance
        # angles           :: angles      in radians
        # frequencies      :: frequencies in Hertz
        # numberOfSpotters :: number of spotters in the optimization
        #
        # DEPENDENCIES:
        #
        import spotter
        import raytracer
        import copy

        #
        # SETTINGS:
        #==============================================================================================
        #
        self.method = 'exact' #'lsmr'  #Solving method choose between:  'lsmr','exact' ,'COBYLA','SLSQP',    
        self.momentsPerSpotter = momentsPerSpotter + 2 * nyquistPenalty
        self.penaltyWeights = [1., 1.,1.,1.,1.]
        for ii in range(0,nyquistPenalty*2):
            self.penaltyWeights.append(1. )
        
        self.nyquistPenalty = nyquistPenalty
        self.kind = kind           

        # Some properties
        #==============================================================================================
        #
        # Associate raytracer object with data assimilation
        self.raytracer = raytracerobject
        #
        # deepcopy angles/frequencies to ensure they remain the same
        angles = raytracerobject.radangles
        self.ang = np.arctan2( np.sin( angles ) , np.cos(angles) ) #atan2 to ensure domain is -pi to pi        

        self.f   = copy.deepcopy( raytracerobject.angfreq )
        self.nf = len(self.f)
        self.nang  = len(self.ang)
        self.xp = xp
        self.yp = yp        
        self.numSpot = len(xp)
        self.nsub = raytracerobject.nsub
        self.svd = None
        self.max_iter = 400 #maximum number of iterations of solvers
        
        #
        # Determine the inward halfsector of the 360 degrees directional distribution
        #==============================================================================================
        if halfCircleStart== None:
            #
            halfCircleStart = -np.pi/2
            #

        # Determine which angles point inwards, we first rotate by the start of the inward halfsector,
        # and make sure that the angles are in the expected range again. After this, the inward angles
        # lie between 0 and pi
        temporary = np.arctan2( np.sin( self.ang - halfCircleStart  ) ,  np.cos( self.ang - halfCircleStart ) )
        temporary[ temporary < 0. ] = temporary[ temporary < 0. ] + 2. * np.pi
        self.inwardAngles = np.array(  [ x > 20 * np.pi/180 and x < 180 * np.pi/180 for x in temporary.tolist()]  )

        print( temporary * 180/np.pi)
        print( self.ang* 180/np.pi)
        print( self.inwardAngles)
        #Based on the kind of decomposition of the directional spectrum, the number of variables in the
        #optimization problem is different        
        if self.kind.lower() in ['energy','energysmooth']:
            #
            self.numVar = len(self.ang[self.inwardAngles] )
            #
        #endif        
        #
        # Create the matrix used in data assimilation
        #==============================================================================================
        self.createMatrix()
        
    def createMatrix( self,activeSpotterList=None ):
        #
        # CreatMatrix:
        #
        # Create the optimization matrix which relates an offshore spectrum to
        # nearshore observations. 
        #
        # if no active list is given, all spotters are presumed to be active
        if activeSpotterList == None:
            #
            activeSpotterList = list( range( 0 , self.numSpot ) )
            #
        #endif

        print( 'ASSIMILATION :: Number of spotters in assimilation is, ' + str(len(activeSpotterList)) )
        # Determine number of active spotters
        self.numSpotActive = len(activeSpotterList)

        # Number of equations depends on ACTIVE spotters        
        self.numEq    = self.momentsPerSpotter * self.numSpotActive

        #
        # create the (inverse) moment matrices for each point
        M  = momentMatrix( 5 , self.ang , self.nyquistPenalty)
        iM = imomentMatrix( self.numVar , self.ang , kind=self.kind,validAngles=self.inwardAngles )
        #
        # Create the transfer matrix that relates offshore to nearshore energies
        # (note if a matrix exists on disk that matrix is loaded from disk)
        T = self.raytracer.getmatrix( self.xp , self.yp , self.ang.tolist(), self.nsub, self.f)       

        #
        # Initialize variables
        self.matrix = np.zeros( ( self.nf , self.numEq , self.numVar) )
        self.Qmatrix = np.zeros( ( self.nf , self.numVar , self.numVar) )        
        #
        #do for each frequency...
        for jf in range( 0,len(self.f) ):
            #
            #and for each point...
            for ip in range(0,self.numSpotActive):
                #
                # Built the transfer matrix from optimization variables to rhs
                jp = activeSpotterList[ip]
                self.matrix[ jf , self.momentsPerSpotter*ip:(self.momentsPerSpotter*(ip+1)) , : ] \
                           = M @ T[jf,jp,:,:] @  iM
                #
            #endfor
            #
            
            # Create the "Q" matrix of the quadriatic form x^T Q x, where in this case Q = 1/2 A^T A
            #
            self.Qmatrix[ jf , : , : ] = 0.5 * np.matmul( self.matrix[jf,:,:].T , self.matrix[jf,:,:] )
            #
        #endfor
        #
    #end def
        
    def optimize( self , spec1d , returnAngularDimension='ang',activeSpotterList=None):
        #
        import scipy
        from scipy import optimize
        import spectral
        import copy
        import code
        #import warnings

        #
        #warnings.filterwarnings("error")        
        #create the matrix
        self.createMatrix( activeSpotterList = activeSpotterList )

        #create the rhs
        rhs,pvec = self.rhs( spec1d , activeSpotterList = activeSpotterList)

        #The inverse matrix
        iM = np.array(imomentMatrix( self.numVar , self.ang, kind=self.kind,validAngles=self.inwardAngles )  )

        #Start the optimization
        E = np.zeros( (self.nf , self.nang) )
        for jf in range( 0 , self.nf):
            #
            print(' - asim: frequency ' + str(jf) + ' of ' + str(self.nf) )
            x0 = np.zeros( (self.numVar) )
    
            #Multiply RHS and Vector with the Penalty Matrix
            RHS = pvec[jf,:] * rhs[jf,:]
            MAT = self.matrix[jf,:,:].copy()
            for ii in range( 0 , self.numEq ):
                #
                MAT[ ii , :] = MAT[ ii , :] * pvec[jf,ii]
                #
            #endfor
            #


            keep      = []
            removeCol = []
            #
            for iv in range( 0 , self.numVar ):
                #
                #
                if False: #meas[iv] < 0.00 * np.mean( meas):
                    #
                    removeCol.append(iv)
                    #
                else:
                    #
                    keep.append(iv)
                    #
                #
            #
            inv = scipy.linalg.pinv( MAT )

            dataresmat = np.matmul( MAT , inv )

            fac = np.diag( dataresmat )
            #RHS = RHS   * fac 
            #MAT = MAT.T * fac
            #MAT = MAT.T


            U, s, V = np.linalg.svd(MAT, full_matrices=True)


            
            csum  = np.cumsum( s / np.sum( s ) )
            s[ csum > 0.9 ] = 0.
            S = np.zeros(MAT.shape)
            S[:MAT.shape[1], :MAT.shape[1]] = np.diag(s)

            MAT = np.dot(U, np.dot(S, V))

            if jf==25:
                print( MAT.shape )
                print( MAT[1,1] )
                self.tinv = inv
                self.tdataresmat = dataresmat
                self.MAT = MAT
            # Remove entries that are insignificant
            #MAT = np.delete( MAT , removeCol , 1 )
            
            #Setup the bounds
            bounds = ( np.zeros( MAT.shape[1] ), np.ones( MAT.shape[1] ) * np.inf)

            #Solve the constrained least squares problem
            res = scipy.optimize.lsq_linear( MAT, RHS,bounds=bounds,lsq_solver=self.method,max_iter=400)
                
                
            #
            #Return the solution
            x       = np.zeros( self.numVar )
            x[keep] = res['x'] 
            E[jf , : ] =  np.matmul( iM , x )
            #
        #endfor
        #
        spec = spectral.spectrum( {spec1d.freqaxis:spec1d.f,'rad':self.ang,'E':E.T} )
        return(spec)

        
    def rhs( self , spec1d , activeSpotterList=None ):
        #
        # Function to create the right hand side
        import spectral        
        
        #the active spotter list denotes the number of each of the active
        #spotters; if no list is given, assume all spotters are active
        if activeSpotterList==None:
            #
            activeSpotterList = list( range( 0 , self.numSpot ) )
            #
        #endif
        #
        #number of equations in the problem
        numEq = len(activeSpotterList) * self.momentsPerSpotter        
        #
        # Create the RHS for the optimization problem
        rhs    = np.zeros( (self.nf , numEq) )
        pvec   = np.zeros( (self.nf , numEq) )        
        for jf in range(0,self.nf):
            #
            for jp in range(0, len( activeSpotterList ) ):
                #
                ip = activeSpotterList[jp]

                if spec1d.E[ ip,jf] > 0.:
                    #
                    pvec[ jf , jp * self.momentsPerSpotter + 0 ] = self.penaltyWeights[0]
                    pvec[ jf , jp * self.momentsPerSpotter + 1 ] = self.penaltyWeights[1]
                    pvec[ jf , jp * self.momentsPerSpotter + 2 ] = self.penaltyWeights[2]
                    pvec[ jf , jp * self.momentsPerSpotter + 3 ] = self.penaltyWeights[3]
                    pvec[ jf , jp * self.momentsPerSpotter + 4 ] = self.penaltyWeights[4]

                    rhs[ jf , jp * self.momentsPerSpotter + 0 ] = spec1d.E[  ip , jf ]
                    rhs[ jf , jp * self.momentsPerSpotter + 1 ] = spec1d.a1[ ip , jf ] *  spec1d.E[ ip , jf  ]
                    rhs[ jf , jp * self.momentsPerSpotter + 2 ] = spec1d.b1[ ip , jf ] *  spec1d.E[ ip , jf  ]
                    rhs[ jf , jp * self.momentsPerSpotter + 3 ] = spec1d.a2[ ip , jf ] *  spec1d.E[ ip , jf  ]
                    rhs[ jf , jp * self.momentsPerSpotter + 4 ] = spec1d.b2[ ip , jf ] *  spec1d.E[ ip , jf  ]
                    
                    for iN in range(0,self.nyquistPenalty):
                        #
                        pvec[ jf , jp * self.momentsPerSpotter + 5 + 2*iN  ] = self.penaltyWeights[5+2*iN]
                        pvec[ jf , jp * self.momentsPerSpotter + 6 + 2*iN  ] = self.penaltyWeights[6+2*iN]
                        rhs[ jf , jp * self.momentsPerSpotter + 5 + 2*iN ]  = 0.
                        rhs[ jf , jp * self.momentsPerSpotter + 6 + 2*iN  ] = 0.
                        #
                    #end for
                    #
                #end if >0
                #
            #endfor over active spotters
            #
        #endfor over frequencies
        #
        return(rhs , pvec)
    #edef
    #
#eclass
    
def momentMatrix( numvar , angles, nyquistPenalty=0 ):
    #
    # Create a block moment matrix M - used to calculate the directional moments
    # from the spectrum
    #
    # input: 
    # n      :: order of the moment matrix
    # angles :: angles for which we calculate the matrix
    #
    # output
    #
    # matrix M = 5 by n matrix
    #
    numdir = len( angles )
    dang = angles[2]-angles[1]    
    M = np.ones( (1,numdir) ) * dang

    #
    n = (numvar - 1) // 2
    for i in range(0,n):
        #
        s = np.sin( (i+1)*angles) * dang   
        c = np.cos( (i+1)*angles) * dang   
        M = np.concatenate( ( M , c[None,:] ) , axis=0 )
        M = np.concatenate( ( M , s[None,:] ) , axis=0 )
        #

    if nyquistPenalty>0:
        #
        for i in range(0,nyquistPenalty):
            iNyq = len(angles) // 2 - i
            s = np.sin( iNyq*angles) * dang
            c = np.cos( iNyq*angles) * dang
            M = np.concatenate( ( M , c[None,:] ) , axis=0 )
            M = np.concatenate( ( M , s[None,:] ) , axis=0 )
        #
    return(  M  )
    #
#enddef

def imomentMatrix( numvar , angles , kind='energy',validAngles=None):
    #
    # The inverse moment matrix - used to calculate the directional
    # spectrum from the moments
    #
    # input: 
    # n      :: order of the moment matrix
    # angles :: angles for which we calculate the matrix
    #
    # output
    #
    # matrix M^-1 = n by 5 matrix
    #
    if kind.lower() == 'energy':
        #
        #
        ifound = 0
        M = np.zeros( (len(angles) , numvar ) )
        for iang in range( 0 , len(angles) ):
           #
           if validAngles[iang]:
               #
               M[iang,ifound] = 1.
               ifound = ifound+1               
           #
        #
    elif kind.lower() == 'energysmooth':
        #
        ifound = 0
        M = np.zeros( (len(angles) , numvar ) )
        for iang in range( 0 , len(angles) ):
           #
           if validAngles[iang]:
               #
               if ifound==0:
                   #
                   M[iang,ifound] = 0.6
                   M[iang,ifound+1] = 0.2
                   #
                   if numvar==len(angles):
                       #
                       M[iang,numvar] = 0.2
                       #
                   #endif                   
                   #
               elif ifound == numvar - 1:
                   #
                   M[iang,ifound  ] =1.
                   M[iang,ifound-1] = 0.
                   #
                   if numvar==len(angles):
                       #
                       M[iang,0] = 0.
                       #
                   #endif
                   #
               else:                
                   #
                   #
                   M[iang,ifound+1] = 0.2
                   M[iang,ifound  ] = 0.6
                   M[iang,ifound-1] = 0.2
                   #
               #endif
               #
               ifound = ifound+1
               #
            #endif
            #
        #endfor
        #
    #endif
    return( M )
    #

def costFunction( x  , Q , A , rhs ):
    #
    vec =  np.matmul( A , x ) - rhs
    
    return( .5 * np.inner( vec , vec ) )
    #

def Jacobian( x , Q , A , rhs ):
    #
    # Jacobian vector of a standard quadratic form 1/2 x^T Q x - x^T A^T rhs + 1/2 rhs^T rhs that follows
    # from a least squares formulation
    #
    return( np.matmul( x , Q.T ) - np.matmul( A.T , rhs ) )

def Hessian( x , Q , A, rhs):
    #
    return(Q)
    #
#

def penaltyMatrix( numeq , weights ):
    #
    # Construct a penalty matrix; basically a higher penalty is associated with energies
    # as with directional moments

    numweights = len(weights)
    if numweights==numeq:
        #
        M = np.diag( weights ) 
        #
    else:
        #
        M = np.zeros( (numeq , numeq) )
        j = 0
        for i in range( 0 , numeq ):
            #
            M[i,i] = weights[j]
            j += 1
            if j == numweights:
                #
                j = 0
                #
            #endif
            #
        #endfor
        #
    #endif
    #
    return( M )
    #
#end def penaltyMatrix
