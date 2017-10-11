import numpy as np

class assimilation:
    #
    # The assimilation object that performs the data assimilation...
    #
    def __init__( self , raytracerobject , xp , yp,
                      momentsPerSpotter = 5,kind='energy' ):
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
        self.momentsPerSpotter = momentsPerSpotter
        self.penaltyWeights = [1., 1.,1.,1.,1.]        
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
        #
        # Create the transfer matrix that relates offshore to nearshore energies
        # (note if a matrix exists on disk that matrix is loaded from disk)
        #==============================================================================================
        self.transferMatrix = self.raytracer.getmatrix( self.xp , self.yp ,
                                                            self.ang.tolist(), self.nsub, self.f)  
        #
    #enddef
    #
        
    def createMatrix( self, jf, activeSpotterList=None ):
        import scipy
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
        self.numVar = len( self.ang )
        
        # Determine number of active spotters
        self.numSpotActive = len(activeSpotterList)

        # Number of equations depends on ACTIVE spotters        
        self.numEq    = self.momentsPerSpotter * self.numSpotActive

        # We start under the assumption that all angles are inward angles
        self.inwardAngles = np.array(  [ True for x in range(0,self.nang) ]  )
        
        #
        # create the (inverse) moment matrices for each point
        #
        M  = momentMatrix( 5 , self.ang )
        iM = imomentMatrix( self.numVar , self.ang , kind=self.kind,validAngles=self.inwardAngles )
     
        #
        # Initialize variables
        matrix  = np.zeros( (self.numEq , self.numVar) )
        scaling = np.zeros( (self.numEq ) )        
        #
        for ip in range(0,self.numSpotActive):
            #
            # Built the transfer matrix from optimization variables to rhs
            jp = activeSpotterList[ip]

            mat = M @ self.transferMatrix[jf,jp,:,:] @  iM
            #
            # Use the maximum matrix entry as a scaling factor
            #
            factor = np.max( np.abs( mat[ 0 , : ] ) )
            scaling[ self.momentsPerSpotter*ip:(self.momentsPerSpotter*(ip+1)) ] = 1. / factor
            
            matrix[ self.momentsPerSpotter*ip:(self.momentsPerSpotter*(ip+1)) , : ] \
                       = (1 / factor ) * mat
            #
        #endfor
        #



        inv = scipy.linalg.pinv( matrix )
        self.inwardAngles[  np.diag(np.matmul(  inv , matrix )) < 0.1 ] = False

        maxconfidence = np.sum( self.inwardAngles ) / self.nang
        
        matrix = matrix[ : , self.inwardAngles ]
        self.numVar = matrix.shape[1]

        uniqueness = matrix[  0:-1:5,: ] @ matrix[0:-1:5,:].T
        scale      =  np.sqrt(np.diag( uniqueness) )
        
        for irow in range( 0 , uniqueness.shape[0] ):
            uniqueness[irow ,:]   = uniqueness[irow ,:]/ scale / scale[irow]
            uniqueness[irow,irow] = np.NaN

        #uniqueness = np.sum( uniqueness , -1 )/self.numSpotActive


        norm =  np.max(np.sum( matrix[0:-1:self.momentsPerSpotter,:] , axis=0))
        confidence = np.sum( matrix[0:-1:5,:] , axis=0) / norm  * maxconfidence
        
        return( matrix , scaling ,confidence,uniqueness )
        #
    #end def

    def optimizeFreq( self , MAT , RHS ):
        #
        import scipy
        from scipy import optimize
        #
        #The inverse matrix
        iM = np.array(imomentMatrix( self.numVar , self.ang, kind=self.kind,validAngles=self.inwardAngles )  )        
        #

        #x0 = np.zeros( (self.numVar) )

        #Setup the bounds
        bounds = ( np.zeros( MAT.shape[1] ), np.ones( MAT.shape[1] ) * np.inf)

        #Solve the constrained least squares problem
        res = scipy.optimize.lsq_linear( MAT, RHS,bounds=bounds,lsq_solver=self.method,max_iter=400)
        x = res['x']
        return( np.matmul( iM , res['x'] ) )
        #
    #enddef
    #

    #
    def optimize( self , spec1d , returnAngularDimension='ang',activeSpotterList=None,
                      removeSingularValues=0.9,method='norm',ModelPriorSpec=None,regval=(1.,1.,1.,1.),jfreq=-1 ):
        #
        import scipy
        from scipy import optimize
        import spectral
        import copy
        import code
        #import warnings

        #
        #Start the optimization

        if activeSpotterList == None:
            #
            activeSpotterList = list( range( 0 , self.numSpot ) )
            #
        #endif
        
        self.numSpotActive = len(activeSpotterList)
        E = np.zeros( (self.nf , self.nang) )
        self.confidenceCurve = np.zeros( (self.nf , self.nang) )
        self.uniqueness = np.zeros( (self.nf , self.numSpotActive , self.numSpotActive) )        

        if jfreq==-1:
            #
            jrange = range(0,self.nf)
            #
        else:
            #
            jrange = range(jfreq,jfreq+1)
            #
        
        #
        for jf in jrange:
            #
            # Create the matrix/rhs
            #
            matrix,scaling,confidence,uniqueness = self.createMatrix( jf , activeSpotterList = activeSpotterList )
            rhs    = self.rhs( spec1d , jf , activeSpotterList = activeSpotterList )
            rhs    = rhs * scaling


            self.confidenceCurve[jf ,  self.inwardAngles] =confidence
            self.uniqueness[jf ,  : , :] = uniqueness

            #Remove Small singular values
            matrix = removeSmallSingularValues( matrix , criterium = removeSingularValues )

            if ModelPriorSpec is not None:
                #
                ModelPrior = ModelPriorSpec.E[ self.inwardAngles , jf ]
                #
            else:
                #
                ModelPrior = (0,)
                #
            #
            
            #Apply regularization
            matrix,rhs = regularization( matrix , rhs , ModelPrior=ModelPrior , methods=method ,
                                             lambda1 = regval[0], lambda2 = regval[1], lambda3 = regval[2],
                                             lambda4 = regval[3],confidence=confidence  )
            
            #
            # Do the inversion at a certain frequency            
            #
            E[jf , : ] =  self.optimizeFreq(  matrix , rhs )
            #
        #endfor
        #
        spec = spectral.spectrum( {spec1d.freqaxis:spec1d.f,'rad':self.ang,'E':E.T} )
        return(spec)
        #
    #end def
    #

    #
    def rhs( self , spec1d ,jf, activeSpotterList=None ):
        #
        #Function to create the right hand side
        
        #the active spotter list denotes the number of each of the active
        #spotters; if no list is given, assume all spotters are active
        #
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
        rhs    = np.zeros( ( numEq) )
        #
        for jp in range(0, len( activeSpotterList ) ):
            #
            ip = activeSpotterList[jp]
            rhs[ jp * self.momentsPerSpotter + 0 ] = spec1d.E[  ip , jf ]
            rhs[ jp * self.momentsPerSpotter + 1 ] = spec1d.a1[ ip , jf ] *  spec1d.E[ ip , jf  ]
            rhs[ jp * self.momentsPerSpotter + 2 ] = spec1d.b1[ ip , jf ] *  spec1d.E[ ip , jf  ]
            rhs[ jp * self.momentsPerSpotter + 3 ] = spec1d.a2[ ip , jf ] *  spec1d.E[ ip , jf  ]
            rhs[ jp * self.momentsPerSpotter + 4 ] = spec1d.b2[ ip , jf ] *  spec1d.E[ ip , jf  ]
            #
        #endfor over active spotters
        #
        return(rhs)
        #
    #end def

    #
    def optimizeRegularizationCoeficients( self , meanAngles , salt = 0.1,seed=0,removeSingularValues=0.9,
                                               method='norm', boundaryArguments=None,searchBounds=None,perFreq=False ):
        #
        #
        #
        if searchBounds is None:
            #
            searchBounds = [ 0 , 10 ]
            #
        #endif
        #

        #
        # Arguments to determine the shape of the boundary spectrum in frequency and directional space
        #
        if boundaryArguments is None:
            #
            boundaryArguments = {'freqshape':'white', 'dirshape':'guassian', 'dirWidth':10 * np.pi / 180.}
            #
        #endif

        if perFreq:
            #
            # Do we want to optimize the reg. coef. for each frequency? If so...
            #
            nf = len(self.f)
            #
        else:
            #
            nf=1
            #
        #endif
        
        par = np.zeros(nf)
        RMS = np.zeros(nf)
        numInt = 5
        for jf in range(0,nf):
            #
            rms     = np.zeros( 11 )            
            #
            for kk in range( 0 , 4):
                #                             
                regvals = np.linspace(searchBounds[0], searchBounds[1],numInt) # [ searchBounds[0], medval ,searchBounds[1] ]
                
                jfreq=-1
                if perFreq:
                    jfreq=jf
                
                for ik,regval in enumerate(regvals):
                    #
                    rms[ik] = self.rmseCostFunction( regval , meanAngles , salt ,seed,
                                                         removeSingularValues,  method , boundaryArguments,jfreq=jfreq )
                    #

                iim   = max( np.argmin( rms ) - 1 , 0 )
                iimin = np.argmin( rms ) 
                iip   = min( np.argmin( rms ) + 1 , numInt-1) 
             
                #
                searchBounds = [  (regval[iim] + regval[iimin]) / 2, (regval[iip] + regval[iimin]) / 2 ]
                par[jf] = regvals[iimin]
                RMS[jf] = rms[iimin]                     
                #
                print( kk,searchBounds,par[jf],RMS[jf])
                #
            #endfor
            #
        #endfor
        #
        return(par,RMS)

    def rmseCostFunction( self , regval ,  meanAngles ,salt ,seed, removeSingularValues,  method , boundaryArguments, jfreq=-1):
        #
        import spectral
        
        rms = np.zeros( len( self.f ) )

        optimizeArguments = {
            'removeSingularValues':removeSingularValues,
            'method':method,
            'jfreq':jfreq
            }        
        if method == 'norm':
            #
            optimizeArguments['regval'] = ( regval , 0,0,0 )
            #
        elif method == 'exposure':
            #
            optimizeArguments['regval'] = ( 0 , 0,0,regval )
            #

            
        
        #
        for meanAngle in meanAngles:
            #
            boundarySpectrum = spectral.parSpec( self.f / np.pi / 2., ang=self.ang,
                                         meanAng=meanAngle,**boundaryArguments)

            error , ___ = self.rmse( boundarySpectrum , salt,seed , **optimizeArguments)
            rms = rms + error
            #
        #end for
        #

        #
        if jfreq > -1:
            #
            return( rms[ jfreq ]  )
            #
        else:
            #
            return( np.sum(rms) )
            #
        #endif
        #
    #end def
    
    def rmse( self , boundarySpectrum , salt=0.1,seed=0 , **optimizeArguments ):
        #
        # this routine calculates the rmse in trying to reconstructing a boundary spectrum
        # if the "observations" are obtained from a forward prediction of a given boundary
        # spectrum which is "salted" with random errors to mimick observational uncertainty

        #Create the observations:
        obs= self.raytracer.prediction(  boundarySpectrum , self.xp , self.yp , 20  )
        obs= obs.Ef()

        #salt the observations with random errors
        if salt>0.:
            #
            obs=obs.salt(salt,seed=seed)
            #

        #do a data-assimilation step
        asim = self.optimize( obs, **optimizeArguments )

        #calculate rms error compared with boundary
        rms = np.sqrt( np.sum( (boundarySpectrum.E - asim.E)**2,axis=0 )/len(asim.ang) )

        # return error and asimiated spectrum
        return( rms , asim )
        #
    #enddef
    #
#eclass
    
def momentMatrix( numvar , angles ):
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

def removeSmallSingularValues( matrix , criterium = 0.9 ):
    #
    # This criterium removes small singular value contributions from the matrix, 
    # and returns the result + rank
    #
    U, s, V = np.linalg.svd( matrix , full_matrices=True)
    #
    # scale the ordered singular values
    #
    relcumsum  = np.cumsum( s )  / np.sum( s )
    #
    #print( np.max(s)/min(s) )

   # F =10
 #   eps = ( np.max(s) - F * np.min(s)  ) / (F - 1)
   # s = s + eps
    
    s[ relcumsum > criterium ] = 0.
    S = np.zeros( matrix.shape)
    S[:len(s), :len(s)] = np.diag(s)
    #
    #return the resulting matrix
    return ( np.dot(U, np.dot(S, V)) )
#end def

def regularization( Mat , Rhs , ModelPrior = (0.,) , methods = 'norm', lambda1 = 1. , lambda2 = 1., lambda3 = 1.,lambda4=1. ,weights=( -1,-1),confidence=0. ):
    #
    numeq = len(Rhs) 
    
    if not type(methods) == list:
        #
        methods = [methods]
        #
    #

    methods = [ x.lower() for x in methods ]

    #
    for method in methods:
        #
        if not ( method in ['norm','slope','curvature','smooth','confidence','none','realconfidence','exposure'] ):
            #
            print( 'WARNING: unknown regularization method: ' + method, ' method is ignored' )
            #
        #
    #endfor
    #
    
    #
    if 'norm' in methods:
        #
        # Norm regularization (energy)
        #
        ncols = Mat.shape[1]
        
        if len(ModelPrior) == 1:
            ModelPrior = np.zeros( ncols )
            
        for i in range( 0 , ncols ):
            #
            vec    = np.zeros( ncols )
            vec[i] = lambda1
            Rhs    = np.append( Rhs , lambda1 * ModelPrior[i] )
            Mat    = np.concatenate( ( Mat , vec[None,:] ) , axis=0 )                
            #
        #endfor
        #
    #endif
    #
        
    #
    if 'slope' in methods:
        #
        # SLope Variance regularization
        #
        ncols = Mat.shape[1]
        
        for i in range( 1 , ncols ):
            #
            vec      = np.zeros( ncols )
            vec[i  ] =  lambda2
            vec[i-1] = -lambda2
            Rhs      = np.append( Rhs , 0. )
            Mat      = np.concatenate( ( Mat , vec[None,:] ) , axis=0 )
            #
        #endfor
        #
    #endif
    #
        
    #
    if 'curvature' in methods:
        #
        # Curvature regularization
        #
        ncols = Mat.shape[1]
        
        for i in range( 1 , ncols-1 ):
            #
            vec      = np.zeros( ncols )
            vec[i  ] =  lambda3
            vec[i-1] = -2*lambda3
            vec[i+1] =  lambda3          
            Rhs      = np.append( Rhs , 0. )
            Mat      = np.concatenate( ( Mat , vec[None,:] ) , axis=0 )                
            #
        #endfor
        #
    #endif                
    #

    #
    if 'smooth' in methods:
        #
        # Curvature regularization
        #
        ncols = Mat.shape[1]

        if np.sum(weights) < 0:
            #
            weights = np.array( [  0.1 , 0.8, 0.1  ] )
            #
        numweights = len(weights)
        
        for i in range( 0 , ncols ):
            #
            vec      = np.zeros( ncols )

            nstart = i - (numweights-1)//2
            nend   = i + (numweights-1)//2 + 1
            jstart = 0
            jend   = numweights
            jhalf  =(numweights-1)//2
            
            if nstart<0:
                nstart = 0
                jstart = jhalf - i

            if nend>ncols:
                nend = ncols
                jend = jhalf + ncols - i

            vec[ nstart : nend ] = weights[ jstart : jend ] / np.sum( weights[ jstart : jend ] ) * lambda4
            #vec[ i ]             = vec[ i ] + lambda4
                
            Rhs      = np.append( Rhs , 0. )
            Mat      = np.concatenate( ( Mat , vec[None,:] ) , axis=0 )                
            #
        #endfor
        #
    #endif                
    #

    if 'exposure' in methods:
        #
        #
        # Get the relative contribution of each observation to the array
        relativeEnergy = Rhs[0:numeq:5]/np.sum(Rhs[0:numeq:5])

        #
        # Weighted sum over all coeficients of the matrix that relate offshore
        # energy at different angles to nearshore observations to create an
        # exposure vector for each direction
        #
        exposureVector = np.zeros( Mat.shape[1] )
        #
        for irow in range( 0 , len(relativeEnergy) ):
            #
            exposureVector = exposureVector + relativeEnergy[irow] * Mat[ irow*5 , : ]
            #
        #
        # normalize to sum==1
        exposureVector = exposureVector / np.sum(exposureVector)
        #
        #
        # Add regularization factor
        for i in range( 0 , Mat.shape[1] ):
            #        
           vec      = np.zeros( Mat.shape[1] )
           vec[i]   = ( 1. - exposureVector[i] ) * lambda4
           
           Rhs      = np.append( Rhs , 0. )
           Mat      = np.concatenate( ( Mat , vec[None,:] ) , axis=0 )    
        #
    return( Mat , Rhs )
    #
#end def
    
#def optimizeRegularizationCoeficients( asim ,  ):
