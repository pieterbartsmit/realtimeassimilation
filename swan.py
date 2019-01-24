#
# Python package to do swan stuff
#
import code
import numpy as np

class tswanboun:
    #
    #
    def __init__( self , kind='file' , par='W' ):
        #
        import spectral
        #
        if kind.lower() == 'file':
            #
            self.side     = par
            self.spec     = []
            self.kind     = kind.upper()
            self.filename = par+'.spc'
            #
        elif kind.lower() == 'nest':
            #
            self.kind = kind.upper()
            self.filename = par
            #
        #
    
class compgrid:
    #
    def __init__( self ):
        #
        # Computational area [ Lon (deg), Lat (deg), Width (deg) , Height (deg) ]
        #   #rect = [ -121.019 , 34.5,     0.5 , 1 ]
        #domain = [ -122.81 , 37.4,     0.5 , 0.5 ]  
        self.domain = [ -122.81 , 37.4,     0.5 , 0.5 ]    
        self.mxc    = 100  #Number of meshes (Longitude)
        self.myc    = 200  #Number of meshes (Latitude)
        self.alpc   = 0.   #Number of meshes (Latitude)
        self.mdc    = 36   #Directional resolution
        self.flow   = 0.03 #Lowest  resolved frequency (Hz)
        self.fhigh  = 0.5  #Highest resolved frequency (Hz)
        self.name   = ''
    #
#END CLASS COMPGRID

class frame:
    #
    def __init__( self , grid ):
        #
        self.name  = ''
        self.xpfr  = grid.domain[0]
        self.ypfr  = grid.domain[1]
        self.xlenfr  = grid.domain[2]
        self.ylenfr  = grid.domain[3]        
        self.mxfr  = grid.mxc
        self.myfr  = grid.myc
        self.alpfr = 0.        
        #
    #
#END CLASS FRAME

class points:
    #
    def __init__( self , xp , yp, name ):
        #
        self.name  =name
        self.xp  = xp
        self.yp  = yp
        #
    #
#END CLASS FRAME

class tPhysics:
    #
    def __init__( self , default=True):
        #
        self.gen3 = default
        self.wcap = default
        self.quad = default
        self.wind = default
        self.breaking = default
        self.fric = default 

class output:
    #
    def __init__( self , data, name,variables ):
        #
        import spotter
        #        
        if type(data) == compgrid:
            #
            self.outputloc = frame( data )
            self.name = name
            self.variables = variables
            self.filename = name + '.mat'
            #
        elif type(data) == list:
            #
            self.outputloc = []
            #
            for elem in data:
                #
                if type( elem ) == spotter.spotres:
                    #
                    self.outputloc.append( points(elem.lon , elem.lat, elem.id) )
                    #
                #
            #
            self.name = name
            self.variables = variables
            self.filename = name + '.tbl'
            self.specfilename = name + '.spc'
        #
    #
#END CLASS OUTPUT


def writeCurrents( filename , Ux , Uy ):
    #
    # Write a current in a Swanreadable fashion; by default we specify the current
    # using a rectangular array where the first index is located on the left/lower (South-East) corner
    # of the grid, and values in Ux[ ix , iy ] increase along the E-W axis for increasing iy
    #
    vel = np.concatenate( [Ux.T , Uy.T] , axis=0 )
    np.savetxt( filename , vel, fmt='%.4e')
    #
#
def writeInput( filename , cgrid, bot, bounlist, wind, tide, outputlist, nestedgridlist, physics=None, inputFields=None ):

    #
    # PBS - Jul, 2017
    #
    # This function creates a command input text file to be used by SWAN; for details on the commands 
    # see the swan manual at http://swanmodel.sourceforge.net/

    #
    if physics == None:
        #
        # Set physics to default physics if not provided.
        #
        physics = tPhysics
        #

    #
    #================================================================================
    # 1) STARTUP ( http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html)
    #================================================================================
    #
    fid = open( filename , 'w' )
    fid.write( '$SWAN INPUT FILE\n')
    fid.write( "$------------------------------------------------------ \n" )
    fid.write( '\n')        
    fid.write( '$This file was automatically generated\n')    

    #
    # Select mode of the simulation as stationary
    fid.write( "$ INIT \n" )
    fid.write( "$------------------------------------------------------ \n" )        
    fid.write( 'MODE STAT TWOD \n')
    fid.write( '\n')

    #
    # Set the tidal elevation
    fid.write( 'SET level={:5.2f} \n'.format(tide['z']) )
    fid.write( '\n')
    #
    # Select global coordinate system as spherical
    fid.write( "$ COMPUTATIONAL GRID \n" )
    fid.write( "$------------------------------------------------------ \n" )        
    fid.write( 'COORDINATES SPHE CCM\n')

    #
    #================================================================================
    # 2) MODEL DESCRIPTION
    #    ( http://swanmodel.sourceforge.net/online_doc/swanuse/node23.html)
    #================================================================================
    #
    #
    # (2a) COMPUTATIONAL GRID
    # (http://swanmodel.sourceforge.net/online_doc/swanuse/node25.html)
    #--------------------------------------------------------------------------------
    fid.write( 'CGRID REG ')
    #Geographical grid
    fid.write( ' xpc={:8.3f} ypc={:8.3f} alpc={:8.3f} xlenc={:8.3f} ylenc={:8.3f}'.format(
        cgrid.domain[0] , cgrid.domain[1], 0.0, cgrid.domain[2], cgrid.domain[3] ) )
    #Geographical Resolution
    fid.write( ' mxc={:4.0f} myc={:4.0f}'.format( cgrid.mxc , cgrid.myc ) )
    #Directional Resolution
    fid.write( ' & \n CIRcle mdc={:4.0f}'.format( cgrid.mdc ) )
    #Frequency grid
    fid.write( ' flow={} fhigh={}'.format( cgrid.flow, cgrid.fhigh ) )
    fid.write( '\n' )
    fid.write( '\n' )

    #
    # (2b) INPUT GRIDS AND DATA
    # (http://swanmodel.sourceforge.net/online_doc/swanuse/node26.html)    
    #--------------------------------------------------------------------------------
    #BATHYMETRY:
    fid.write( "$ INPUT GRIDS \n" )
    fid.write( "$------------------------------------------------------ \n" )    
    fid.write( 'INP BOT')
    #Geographical grid
    fid.write( ' xpinp={:8.3f} ypinp={:8.3f} alpinp={} mxinp={} myinp={} dxinp={} dyinp={}' \
      .format( bot['xpinp'],bot['ypinp'],0.0,bot['mxinp'],bot['myinp'],bot['dxinp'],bot['dyinp'] ) )
    fid.write( '\n' )    
    #
    fid.write( "READINP BOT 1. '{}' idla=3 FREE".format(bot['filename']) )    
    fid.write( '\n' )

    #
    # INPUTFIELDS
    if inputFields is not None:
        #
        for ip in inputFields:
            #
            # First write the input grid
            #
            fid.write( 'INP {}'.format(ip['kind']))
            fid.write( ' xpinp={:8.3f} ypinp={:8.3f} alpinp={} mxinp={} myinp={} dxinp={} dyinp={}' \
            .format( ip['xpinp'],ip['ypinp'],0.0,ip['mxinp'],ip['myinp'],ip['dxinp'],ip['dyinp'] ) )
            fid.write( '\n' )

            #
            # And then the command to read the data
            #
            fid.write( "READINP {} 1. '{}' idla=3 FREE".format(ip['kind'],ip['filename']) )
            #



    #
    # (2c) BOUNDARY CONDITIONS
    # (http://swanmodel.sourceforge.net/online_doc/swanuse/node27.html)    
    #--------------------------------------------------------------------------------
    #
    fid.write( "\n" )
    fid.write( "$  BOUNDARY \n" )
    fid.write( "$------------------------------------------------------ \n" )    
    for boun in bounlist:
        #
        if boun.kind.upper() == 'FILE':
            #            
            fid.write( "BOUN SIDE {} CON FILE '{}'\n".format(boun.side,boun.filename) )
            #
        elif boun.kind.upper() == 'NEST':
            #
            fid.write( "BOUN NEST '{}.nest'\n".format(boun.filename) )
            #
        #
    #
            
    #
    # (2d) PHYSICS
    # (http://swanmodel.sourceforge.net/online_doc/swanuse/node28.html)    
    #--------------------------------------------------------------------------------
    #
    if physics.wind:
        #
        fid.write( "WIND vel={} dir={}".format(wind['U'],wind['Udir']) )
        #
    #
    fid.write( "\n" )
    fid.write( "$  PHYSICS \n" )
    fid.write( "$------------------------------------------------------ \n" )

    if physics.gen3:
        #
        fid.write( "GEN3 KOMEN \n" )

    if physics.wcap:
        #        
        fid.write( "WCAP KOMEN \n" )
    else:
        fid.write( "OFF WCAP \n" )

    if physics.quad:
        #        
        fid.write( "QUAD \n" )
    else:
        fid.write( "OFF QUAD \n" )

    if physics.fric:
        #        
        fid.write( "FRIC JON \n" )

    if physics.breaking:
        #        
        fid.write( "BREAK CON  \n" )
        
    #fid.write( "NUM STOPC STAT mxitst=2  \n" )    
    #
    # OUTPUT
    #
    fid.write( "\n")
    fid.write( "$  OUTPUT \n" )
    fid.write( "$------------------------------------------------------ \n" )
    #
    for out in outputlist:
        #
        if type(out.outputloc) is frame:
            #
            string = "FRAME '{}' xpfr={:8.3f} ypfr={:8.3f} alpfr={:8.3f} xlenfr={:8.3f} ylenfr={:8.3f} mxfr={:4.0f} myfr={:4.0f}"
            string = string.format( out.name, out.outputloc.xpfr , out.outputloc.ypfr, out.outputloc.alpfr,\
                                    out.outputloc.xlenfr, out.outputloc.ylenfr, out.outputloc.mxfr , out.outputloc.myfr)
            fid.write( string + '\n' )
            fid.write( "BLOCK '{}' NOHEADER '{}'".format(out.name , out.filename) )
            #
            for variable in out.variables:
                #
                fid.write( ' '+variable )
                #
            #end for
            #
            fid.write( '\n' )
            #
        elif type(out.outputloc) is list:
            #
            fid.write("POINTS '{}'".format( out.name ))
            #
            for point in out.outputloc:
                #code.interact( local=locals() )
                fid.write(" & \n ")
                fid.write(" {:8.3f} {:8.3f}".format( point.xp , point.yp ) )
                #
            #
            fid.write( '\n')            
            fid.write( '\n')
            fid.write("TAB '{}' HEAD '{}' ".format( out.name , out.filename ) )
            #
            for variable in out.variables:
                #
                fid.write( ' '+variable )
                #
            #end for
            #
            fid.write( '\n' )
            fid.write("SPEC '{}' SPEC2D ABS '{}' \n".format( out.name , out.specfilename ) )            
        #end if
        #
    #end for
    #
    # WRITE NESTED GRIDS
    #
    if len(nestedgridlist) > 0:
        #
        fid.write( "\n")
        fid.write( "$  OUTPUT \n" )
        fid.write( "$------------------------------------------------------ \n" )    
        for nest in nestedgridlist:
            #
            fid.write( "NGRID '{}' ".format( nest.name ) )
            fid.write( ' xpn={:8.3f} ypn={:8.3f} alpn={:8.3f} xlenn={:8.3f} ylenn={:8.3f}'. \
                           format(nest.domain[0] , nest.domain[1], nest.alpc, nest.domain[2], nest.domain[3]) )
            fid.write( ' mxn={:4.0f} myn={:4.0f} \n'.format( nest.mxc , nest.myc ) )
            fid.write( "NEST '{}' '{}.nest' \n".format( nest.name , nest.name ) )
            #
        #endfor
        #
    #endif
    #
    # Start computation
    #
    fid.write( "\n")
    fid.write( "$  COMPUTE \n" )
    fid.write( "$------------------------------------------------------ \n" )
    fid.write( "COMPUTE \n" )
    fid.write( "STOP \n" )        
    fid.close()
    
def writeDirectionalSpectrum( spec  , workingdirectory='./', filename='spec.swn', x=None, y=None , lonlat=False ):
    #---------------------------------------------------------------------------    
    #
    # This routine writes a directional spectrum and writes it to a textfile
    # according to the SWAN format guidlines
    #
    # The swan format is described at (as of July, 2017):
    #
    # http://swanmodel.sourceforge.net/online_doc/swanuse/node50.html 
    # 
    # INPUT (required):
    #
    # E   :: numpy array [ndir , nfreq]  Variance density spectrum (units: m2/Hz/Deg)
    # f   :: numpy array [nfreq]         frequency                 (units: Hz)
    # ang :: numpy array [ndir]          directions                (units: degrees)
    # 
    # INPUT (Optional):
    #    
    # filename :: output filename
    # loc      :: numpyarray [2] Location [default=(0. ,0.)] of spectrum (units see latlon)
    # lonlat   :: [=True / False(DEFAULT)]  coordinates of location given as [x,y] or [long,lat] (in that order)  
    #
    # PBS - Juli, 2017
    #---------------------------------------------------------------------------
    #
    # Code description:
    # 1) sanity checking of input
    # 2) grab data from server
    # 3) reformat data to usuable format
    # 4) return frequency/direction spectrum
    #

    if spec is None:
        return

    if x == None:
        #
        x = [0.]
        y = [0.]
        #
    
    E = spec.correctDim()
    f = spec.f
    ang = spec.ang
    
    nfreq = len(f)
    ndir  = len(ang)
    
    fid = open( workingdirectory + '/' + filename , 'w' )
    fid.write( 'SWAN 1\n')
    #
    # LOCATION
    #
    nloc = len(x)
    if lonlat:
        #
        fid.write( 'LONLAT\n')
        #
    else:
        #                
        fid.write( 'LOCATIONS\n')
        #

    fid.write(' {} \n'.format( len(x) ) )
    for i in range( 0, len(x) ):
        #
        fid.write(' {} {} \n'.format(x[i], y[i]) )       
    #
    # WRITE FREQUENCIES
    #
    fid.write('RFREQ\n' )
    fid.write(' {} \n'.format( nfreq ) )
    for i in range( 0, nfreq ):
        #
        fid.write(' {} \n'.format(f[i]) )

    #
    # WRITE DIRECTIONS
    #
    fid.write('CDIR\n' )
    fid.write(' {} \n'.format( ndir ) )
    for i in range( 0, ndir ):
        #
        fid.write(' {} \n'.format(ang[i]) )

    fid.write('QUANT\n' )
    fid.write(' 1\n' )
    fid.write('VaDens\n' )
    fid.write('m2/Hz/degr\n' )
    fid.write(' {}\n'.format(-99.0) )

    for iloc in range( 0 , nloc ):
        #
        fid.write('FACTOR\n' )
        fid.write(' 1.0\n' )
        #
        for i in range( 0, nfreq ):
            #    
            for j in range( 0, ndir ):
                #
                fid.write(' {:.3e}'.format(E[iloc,j,i] )  )
                #
            #
            #end for directions
            #
            fid.write('\n' )
            #        
        #end for frequencies
        #
    #end for locations
    fid.close()
    #

def swanReadSpec(filename):
    import numpy
    import spectral
    #
    # Read a swan spectral file as documented in
    # "Spectrum files, input and output" of the SWAN user manual
    #
    with open( filename , 'r' ) as fid:
        #
        lines  = fid.readlines()
        lines  =[line.replace('\n','') for line in lines]
        lines  =[line.strip() for line in lines]
        #
    iline = -1
    try:
        #
        # read locations
        #
        [line,iline] = swanReadSpec_newLine( iline , lines,  ['LOCATIONS','LONLAT'] )
        [nloc,iline] = swanReadSpec_val( iline , lines  )
        nloc         = int(nloc)
        [loc,iline]  = swanReadSpec_val( iline , lines , 2, nloc )
        #
        # read frequency
        #
        [line,iline]  = swanReadSpec_newLine( iline , lines,  ['RFREQ','AFREQ'] )
        [nfreq,iline] = swanReadSpec_val( iline , lines  )
        nfreq         = int(nfreq)
        [freq,iline]  = swanReadSpec_val( iline , lines , 1, nfreq )
        #
        # read direction
        #        
        [line,iline]  = swanReadSpec_newLine( iline , lines,  ['NDIR','CDIR'] )
        [ndir,iline]  = swanReadSpec_val( iline , lines  )
        ndir          = int(ndir)
        [ang,iline]   = swanReadSpec_val( iline , lines , 1, ndir )
        #
        # read quant
        #      
        [line,iline]  = swanReadSpec_newLine( iline , lines,  ['QUANT'] )
        [nquant,iline]= swanReadSpec_val( iline , lines  )
        nquant        = int(nquant)
        #
        # read unit
        #
        [line,iline]  = swanReadSpec_newLine( iline , lines,  ['VADENS'] )
        [line,iline]  = swanReadSpec_newLine( iline , lines,  ['m2/Hz/degr'] )
        [excepval,iline] = swanReadSpec_val( iline , lines  )

        E = np.zeros([ nloc , ndir , nfreq ])
        #
        for j in range(0,nloc):
            #
            try:
                [line,iline]    = swanReadSpec_newLine( iline , lines,  ['FACTOR'] )
                [factor,iline]  = swanReadSpec_val( iline , lines , n=1 , m=1 )
            except:
                iline = iline+1
                continue
            [temp,iline]    = swanReadSpec_val( iline , lines , n=ndir , m=nfreq )
            temp            = temp.T
            E[ j , : , : ]  = temp * factor
            #
        #
        return( spectral.spectrum( {'E':E,'f':freq,'ang':ang } ) )
        #
    #
    except:
        #
        print( 'error reading ' + filename )
        raise
        #

        
def swanReadSpec_val( iline , lines , n=1 , m=1 ):
    #
    import numpy

    array = np.zeros( [m , n] )
    for j in range( 0 , m ):
        #
        [line,iline] = swanReadSpec_newLine( iline , lines, None )
        vals = [ float(line[index]) for index in range(0,n) ]          
        array[ j , : ] = np.array( vals )
        #
    #
    #
    if n > 1 or m > 1:
        #
        return( array , iline)
        #
    else:
        #
        return( array[0,0] , iline)
        #            
    #
        
def swanReadSpec_newLine( iline , lines, expectedKeywords=None ):
    #
    eof = False
    while (not eof):
        #
        iline = iline + 1
        if iline == len(lines):
            #
            eof = True
            line = None
            raise Exception('unexpected end of file')
            break
            #
        #
        line = lines[iline]
        #
        # Commands we ignore (for now)
        #
        if line[0] == '$':
            #
            continue
            #
        elif line[0:4].upper() == 'SWAN':
            #
            continue
            #
        line = line.split()
        break
    #
    if (not expectedKeywords==None):
        #
        expectedKeywords = [word.upper() for word in expectedKeywords]
        if not (line[0].upper() in expectedKeywords):
            #
            str = 'unexpected keyword:' + line[0] + ' expected: '
            for word in expectedKeywords:
                str = str + word + ' '            
            raise Exception( str )

                
    
    return(  line, iline )
    
def swanRunSuccess(workingdirectory):
    import os
    #
    # Check for errfile
    #
    if os.path.isfile(workingdirectory + '/' +  'Errfile'):
        #
        return(False)
        #
    #
    # check for print + norm_end files
    #
    with open( workingdirectory + '/' + 'PRINT' ) as fid:
        #
        lines  = fid.readlines()[-1]
        #
    #

    #
    if lines.strip() == 'STOP':
        #
        if os.path.isfile(workingdirectory + '/' +  'norm_end'):
            #
            return(True)
            #
        #
    else:
        #
        return(False)
        #
    #
#
    
def readSwanBlock( filename ):
    #
    # Open a swan block data file written as a binary mat file
    #
    import scipy.io as sio
    #
    data = sio.loadmat(filename)
    
    return(data)
#

#
def cleanSwanFiles( directory ):
    #
    import os.path
    import os
    #
    swanfiles = ['INPUT', 'Errfile','norm_end','PRINT']
    #
    for name in swanfiles:
        #
        filename = directory + '/' + name
        if os.path.isfile(filename):
            #
            os.remove(filename)
            #
        #
    #

def readSwanTable( filename ):
    #
    with open( filename , 'r' ) as fid:
        #
        lines  = fid.readlines()
        lines  =[line.replace('\n','') for line in lines]
        lines  =[line.replace('%','') for line in lines]
        header = lines[4].split()
        lines  = lines[7:]
        #
        data = dict()
        for name in header:
            #
            data[name] = []
            #
        #
        
        for line in lines:
            #
            variables = line.split()
            #
            for [index,name] in enumerate(header):
                #
                data[name].append( variables[index] )
                #
            #
        #        
    #
    return(data)
#
