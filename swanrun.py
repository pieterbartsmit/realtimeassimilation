def genDirs(directory,dataAssimilation=False):
    #
    import os
    import platform

    #
    # A rudimentary check to see if I am on my local mac or on the remote server
    # mostly convinient for development; should be removed later
    #
    if platform.system() == 'Darwin':
        #
        swanloc = '/Users/pieterbartsmit/Source/swan/41.10.A/OMP/swan.exe'
        #
    else:
        #
        swanloc = '/mnt/2ndHDD/swan/41.10.A/OMP/swan.exe'
        #
        
    if dataAssimilation:
        #
        prefix = 'asim' + os.path.sep
        #
    else:
        #
        prefix = 'pred' + os.path.sep
        #
    #
    dirs = { 'work':'work',   \
             'outp':prefix+'output', \
             'err':prefix+'err', \
             'inp':prefix+'inp', \
             'outputblock':prefix+'output/block/', \
             'outputpointsbulk':prefix+'output/points/bulk/', \
             'outputpointsspec':prefix+'output/points/spec/', \
             'log':prefix, \
             'bnd2d':prefix+'output/boundary/2d/', \
             'bnd1d':prefix+'output/boundary/1d/', \
             'ray':'raytracer/'
           }
    #
    # Check if the root exists - if so, do nothing
    #
    root = os.path.abspath(directory)
    #
    for key, value in dirs.items():
        #
        dirs[key] = root + os.path.sep + value
        #
    #End for            
    return( dirs , root, swanloc )
    #

def getOutputFilenames( kind , directory,gridname=None,dataAssimilation=False ):
    import os

    if gridname==None:
        gridname = directory
    
    [dirs,root,swanloc] = genDirs(directory,dataAssimilation)
    out = None
    if kind.lower() == 'bnd2d':
        #
        files = os.listdir(dirs['bnd2d'])
        out = {'N':[],'S':[],'E':[],'W':[]}
        #
        for filename in files:
            #
            if filename[-3:].lower() == 'spc':
                #
                side = filename[-5]
                out[side].append( dirs['bnd2d'] + filename )
                #
    elif kind.lower() in ['hspswa','hspray','hspobs']:
        #
        print( kind.lower() )
        files = os.listdir(dirs['outputpointsbulk'])
        #        
        if kind.lower() in ['hspswa']:
            #
            pre = 'hs'
            post = 'sim'
            #
        elif kind.lower() in ['hspray']:
            #
            pre = 'hs'
            post = 'ray'
            #
        elif kind.lower() in ['hspobs']:
            #
            pre = 'hs'
            post = 'obs'
            #  
            
        for filename in files:
            #
            if filename[-3:].lower() == post and filename[0:2].lower() == pre: 
                #
                out = dirs['outputpointsbulk'] + filename
                break
                #
    return(out)
    #
#end def getOutputFilenames


    

class swanrun:
    #
    # Public Methods:
    #
    # __init__     :: constructor
    # run()        :: run a swan simulation
    # setup()      :: setup a swan simulation (must be called before run)
    # addNest(...) :: add a nested simulation
    # update(...)  :: check if spotters in the grid have new data
    #
    #===============================================================================================     
    def __init__( self ,name , directory, cgrid, source, outputlist, spotterlist=None,bounlist=None,
                      dataAssimilation = False , swanopt = None, rayTracer = None, token = None ):
    #===============================================================================================         
        #
        import spotter
        import swan
        import calendar
        import time
        import os
        import raytracer
        import utils
        
        print('- Initialize model: ' + name)

        self.dataAssimilation = dataAssimilation
        self.regval = 1.
        self.asimLimitingAngles = None
        
        if bounlist == None:
            bounlist = []

        if swanopt == None:
            #
            self.physics = swan.tPhysics(default=False)
            self.breaking = True
            #

        self.token=token #The token is used to communicate with the Spoondrift backend
        #
        # Generate Set of working directories
        #        
        [self.dirs,self.root,self.swanloc] = genDirs(directory,dataAssimilation)
        #
        for key, value in self.dirs.items():
            #
            if not os.path.isdir(self.dirs[key]):
                #
                os.makedirs(self.dirs[key])
                #
            #End if
            #
        #End for
        #        

        #
        self.cgrid            = cgrid
        self.source           = source
        self.outputlist       = outputlist


        if spotterlist == None:
            #
            self.spotters = None
            self.spotterlist = []
            #
        else:
            #
            self.spotterlist      = spotterlist.getIds()
            self.spotters         = spotterlist
            #
        #

        self.bounlist = []
        for boun in bounlist:
            #
            if boun.upper() in ['N','S','E','W']:
                #
                self.bounlist.append( swan.tswanboun( par=boun.upper() ) ) 
                #
            
        self.simulationtime   = -1
        self.nestedrunlist = []
        self.nestedgridlist = []        
        self.nestdepth = 0
        self.nestnum   = 0        
        self.parentrun = None
        self.name   = name
        self.commandfilename   = name + '.swn'
        self.data      = []
        self.blockdata = []

        #
        # RETRIEVE Bathymetry
        #-----------------------------------------------------------------------
        workingdirectory = self.dirs['work'] + '/'           
        self.bat = self.retrieveBathymetry( self.cgrid.domain, \
                        workingdirectory=workingdirectory, reloaddata=False)
        
        #
        # Setup raytracer object
        #-----------------------------------------------------------------------
        #
        self.dcont = 5.
        self.rt    = False
        #
        if len( bounlist ) > 0:
            #
            if rayTracer == None:
                #
                self.raytracer = raytracer.raytracer( self.cgrid ,
                                self.bat['lon'][0,:] , self.bat['lat'][:,0] ,
                                self.bat['dep']      , bounlist        ,
                                self.dcont           , workdir=self.dirs['ray'])
                
                #
            else:
                #
                self.raytracer=rayTracer
                #
            #end
            #
            self.rt = True
            # 
            #
        #
        # Generate contour lines
        #
        self.zlevels = [ 0 , 10,20,30,40,50,100,200,300,400,500,600,700,800,900,1000]
        self.contourlines = utils.get_contours( self.bat['lon'] , self.bat['lat'] , self.bat['dep'] ,self.zlevels )
        #
    #End CONSTRUCTOR
    #
    #===============================================================================================   
    def update( self , interval=0):
    #===============================================================================================
        #
        # FUNCTION DESCRIPTION
        #---------------------------------------------------------------------------    
        # This functions checks if updates are available from the spotters,
        # retrieves them, and if updates were available returns true to indicate
        # that we can do a new simulation with new data if desired. The interval
        # parameter (default=0) allows us to wait with a new simulation if the new spotter
        # data doesnt change the mean time all that much (i.e. there is not significantly
        # new data). Note, if the last simulation (infered from self.simulationtime) was
        # at t=0, this means that we havent done any simulations yet, it that case, also
        # return True.
        #
        #---------------------------------------------------------------------------
        upd = False
        #
        if ( ( len( self.spotters.update() )  > 0 ) and 
           ( -self.simulationtime + self.spotters.meanTime > interval ) ) or \
           (self.simulationtime < 0):
            #
            # There is an updated spotter - or we havent done any simulations yet
            #
            self.simulationtime = self.spotters.meanTime
            self.spotters.appendcsv( self.dirs['outputpointsbulk'] + '/' )
            
            upd = True
            #
        #end if
        #        
        return( upd )
        #
    #end def    

    #===============================================================================================    
    def setup( self  ):
    #===============================================================================================
        #
        # FUNCTION DESCRIPTION
        #---------------------------------------------------------------------------    
        # This functions creates the necessary swan input files to perform a new run.
        # 
        # 1) Retrieve bathymetry
        # 2) Retrieve boundary spectra/wind data/tide data
        # 3) Write boundary conditions
        # 4) Generate swan input file
        # 5) call the setup() function for any nested runs
        #
        #---------------------------------------------------------------------------
        import swan
        import numpy
        import calendar
        import time
        import copy
        import spectral
        #
        workingdirectory = self.dirs['work'] + '/'        
        print( '* Setup input files for model: ' + self.name )
    
        # 2) RETRIEVE DATA
        #---------------------------------------------------------------------------

        #
        dat = self.retrieveDataFromServers( self.simulationtime , \
                                                workingdirectory=workingdirectory )
        #self.currents = dat['currents']
        #
        # 3) BOUNDARY CONDITIONS
        #---------------------------------------------------------------------------
        #
        for index,boun in enumerate(self.bounlist):
            #
            if (not (boun.kind.upper() == 'NEST')):
                #
                # assign spectrum
                #    
                self.bounlist[index].spec =dat['dirspec']
                #
                # write spectrum
                #
                swan.writeDirectionalSpectrum(self.bounlist[index].spec, workingdirectory,\
                                                  self.bounlist[index].filename )
                #
            #
        #endif
        #
        # 4) Generate Swan Input
        #---------------------------------------------------------------------------
        #
        outputlist = copy.deepcopy(self.outputlist)
        name       = 'points'
        outputlist.append( swan.output( \
            self.spotters.spotters , name , ['XP','YP','HS','TM01','DIR','DSPR','DEP','PDIR','RTP'] ) )

        if dat['currents'] is not None:
            inputFields = [dat['currents']]
        else:
            inputFields = None

        name = workingdirectory + '/' + self.commandfilename
        swan.writeInput( name       , self.cgrid                           , \
                         self.bat   , self.bounlist                        , \
                         dat['wind'], dat['tide']                          , \
                         outputlist , self.nestedgridlist, physics=self.physics,
                             inputFields=inputFields )                                     
        #
        #
        # 5) Call any nested runs (if applicable)
        #---------------------------------------------------------------------------
        #
        if  (len(self.nestedrunlist) > 0 ):
            #
            for nest in self.nestedrunlist:
                #
                nest.simulationtime = self.simulationtime
                nest.setup()                
                #
            #end for
            #
        #end if nest
        #        
    #end def setup
    #
    
    #
    #===============================================================================================     
    def retrieveBathymetry( self , rect ,workingdirectory='./',reloaddata=False):
    #===============================================================================================
        #
        # FUNCTION DESCRIPTION
        #---------------------------------------------------------------------------
        #
        # This function retrieves the bathymetry from the NDGC servers and saves the
        # results locally as a grid in text files. Note that if these grids already
        # exist, the bathymetry is simply retrieved from local files
        #
        # INPUT:
        #--------
        # 
        # rect = [ xpinp , ypinp,  width , heigth ] :: list describing the bounding box
        #                                              of the region
        # workingdirectory                          :: where to save bat files
        # reloaddata                                :: force reloading from online
        #                                              databased if True
        #
        # OUTPUT
        #--------
        #
        # bat :: dictionary containg the bathymetry and description thereof (see 
        #        below)
        #
        # DEPENDENCIES
        #---------------------------------------------------------------------------
        import ngdc
        import numpy
        import code
        import os.path
        #
        fileexists = os.path.isfile(workingdirectory+'/bat.dep')
        #
        # Get bathymetry from ngdc server
        #---------------------------------------------------------------------------
        if reloaddata or not fileexists:
            #
            print('   - Downloading Bathymetry from ngdc')
            bat=ngdc.getCoastalReliefData( rect , filename=workingdirectory+'bat.tmp')    
            numpy.savetxt( workingdirectory+'bat.dep' , bat['dep'], fmt='%.4e')
            numpy.savetxt( workingdirectory+'lon.coo' , bat['lon'], fmt='%.7e')
            numpy.savetxt( workingdirectory+'lat.coo' , bat['lat'], fmt='%.7e')
            numpy.save( workingdirectory+'bat.npy', bat['dep'])
            numpy.save( workingdirectory+'lon.npy', bat['lon'])
            numpy.save( workingdirectory+'lat.npy', bat['lat'])
            #
        else:
            #
            print('   - Retrieving Bathymetry from disk')
            #
            bat = dict()
            #bat['dep'] = numpy.loadtxt( workingdirectory+'bat.dep' )
            #bat['lon'] = numpy.loadtxt( workingdirectory+'lon.coo' )
            #bat['lat'] = numpy.loadtxt( workingdirectory+'lat.coo' )
            bat['dep'] = numpy.load( workingdirectory+'bat.npy' )
            bat['lon'] = numpy.load( workingdirectory+'lon.npy' )
            bat['lat'] = numpy.load( workingdirectory+'lat.npy' )
            #
        #end if
        #
        # Infer the bathymetry characteristics from the files
        # (see swanmanual, input grid description for clarification of names)
        #
        bat['xpinp'] = bat['lon'][1,1]
        bat['ypinp'] = bat['lat'][1,1]
        bat['nlat']  = numpy.shape( bat['dep'] )[0] 
        bat['nlon']  = numpy.shape( bat['dep'] )[1]
        bat['mxinp'] = numpy.shape( bat['dep'] )[1] - 1
        bat['myinp'] = numpy.shape( bat['dep'] )[0] - 1
        bat['dxinp'] = (bat['lon'][1,2]-bat['lon'][1,1])
        bat['dyinp'] = (bat['lat'][2,1]-bat['lat'][1,1])    
        bat['filename'] = 'bat.dep'
        #
        return( bat )
        #
        #
        #
    #===============================================================================================     
    def retrieveDataFromServers( self ,epochtime , workingdirectory='./' ):
    #===============================================================================================         
        #
        # (2) WW3 Spectra
        # Get directional spectra from buoy
        #

        import wavewatch3
        import coops
        import swan
        import numpy
        import ndbc
        import dataAssimilation

        #
        # SOURCE OF BOUNDARY SPECTRA
        #

        if self.dataAssimilation:
            #
            # Create assimilation object
            #
            x,y = self.spotters.getLoc()
            asim = dataAssimilation.assimilation( self.raytracer , x, y,
                            limitingAngles=self.asimLimitingAngles)
            #
            # Do the data assimilation - and make sure the spectrum is in /deg/Hz
            #
            spec = self.spotters.getSpec()

            if self.source['dirspec'][0].lower() == 'ww3':
                #
                print('   - Downloading directional spectrum from NOAA ww3'
                          +' as model Prior')
                prior = wavewatch3.getSpectrum(
                    epochtime,self.source['dirspec'][1] )
                prior.torad()
                prior.interpFreq( spec.f )
                #dirspec.constraint( -90 , 180)          
                #
            else:
                #
                print('   - No Source Spectrum available as a model Prior')
                prior=None
                #
                # end if

            if len( self.spotters.fullSpecIndices ) < 1:
                #
                # If no recent observed full spectra data is available, use ww3
                # as a fallback
                #
                dirspec = prior
            elif len( self.spotters.fullSpecIndices ) < 3:
                #
                # Assimilation with model prior
                #
                dirspec=asim.optimize(
                    spec , activeSpotterList = self.spotters.fullSpecIndices,
                    removeSingularValues=.9, ModelPriorSpec=prior,
                    method=['norm'],regval=(1.,1.,1.,1.)
                    )
                #                
            else:
                #
                # Assimilation
                #
                dirspec=asim.optimize(
                    spec , activeSpotterList = self.spotters.fullSpecIndices,
                    removeSingularValues=.9,
                    method=['norm'],regval=(self.regval,1.,1.,1.)
                    )
                #
            #endif
            #
            if prior is not None:
                #
                if ( (dirspec.Hm0() > 2. * prior.Hm0() ) and
                         ( dirspec.Hm0() > 1.) ):
                    #Sanity check, if there is a large difference in
                    #predictions fall back to ww3 prior prediction
                    print("Warning, assimilation waveheights conflict")
                    print("  asim waveheight at boundary:"
                              + str(dirspec.Hm0() ) )
                    print("  ww3 waveheight at boundary:"
                              + str(prior.Hm0() ) )
                    print("Using ww3 boundary as failsafe")
                    dirspec = prior
                    #
                #
            #
            dirspec.todeg()
            dirspec.tof()                            
            #
        else:
            #
            try:
                #
                if self.source['dirspec'][0].lower() == 'ww3':
                    #
                    print('   - Downloading directional spectrum from NOAA ww3')
                    dirspec = wavewatch3.getSpectrum(epochtime,self.source['dirspec'][1],workingdirectory)
                    #dirspec.constraint( -90 , 180)                    
                    #
                else:
                    #
                    print('   - No Source Spectrum')
                    #
                # end if
                #
            except:
                #
                self.log('Error downloading boundary spectrum - halting simulation')
                raise Exception('Cannot retrieve ww3 boundary spectrum')
                #
            #end try
            #
        #end else assimilation
        #
        # TIDES
        #
        try:
            #
            if self.source['tide'][0].lower() == 'coops':
                #
                # get tide data
                #
                print('   - downloading tide data from coops')
                tide = coops.getTideData( epochtime, station=self.source['tide'][1],workingdirectory=workingdirectory)
                #
            else:
                #
                tide = dict()
                tide['z'] = 0.
                print('   - No Tidal data')
                #
        except:
            tide = dict()
            tide['z'] = 0.
            self.log('Error downloading tide data - msl set to 0')            
            
        #
        #get wind data
        #
        if self.source['wind'][0].lower() == 'coops':
            #
            print('   - downloading wind data from coops')
            #
            # NOT IMPLEMENTED
            #
            #wind = coops.getData( '' , kind=1 )
            #
        elif self.source['wind'][0].lower() == 'ndbc':
            #
            try:
                #
                print('   - downloading wind data from ndbc')
                wind = ndbc.getMeanWindData( epochtime , length=1800 , \
                        stationID=self.source['wind'][1],workingdirectory=workingdirectory  )
            except:
                #
                wind = {'U':0.,'Udir':0.}
                self.log('Error downloading wind data - wind set to 0')
            #
        else:
            #
            wind = {'U':0.,'Udir':0.}
            #
        #Endif

        if 'currents' in self.source:
            #
            currents = self.source['currents']( epochtime )
            outfile = workingdirectory + currents['filename']
            swan.writeCurrents( outfile , currents['U'] , currents['V'] )
            #
        else:
            #
            currents = None
            #
        #
        return( {'wind':wind , 'tide':tide , 'dirspec':dirspec, 'currents':currents } )
        #


    #End retrieveDataFromServers
    #

    #===============================================================================================     
    def run(self):
    #===============================================================================================         
        #
        import os
        import shutil
        import utils
        import os.path
        import time
        import calendar
        import swan
        import spotter
        import code
        import copy
        import scipy
        import numpy as np
        #
        self.runstart = calendar.timegm( time.gmtime() )
        #
        cwd = os.getcwd()
        #
        # Change working directory
        #
        os.chdir(  self.dirs['work'] )
        #
        # make sure no previous run files exist
        #
        swan.cleanSwanFiles('./')
        
        #
        # run swan
        #        
        shutil.copyfile( self.commandfilename , 'INPUT')
        os.system( self.swanloc )

        
        #
        # - add time stamp to identify simulation
        t    = time.strftime("%Y%m%d%H%M", time.gmtime( self.simulationtime ) )        
        #
        # Was the run succesfull?
        #
        if (not swan.swanRunSuccess('./') ):
            #
            # If not succesfull...
            #
            #copy input file
            srce = 'INPUT'
            dest = self.dirs['err'] +'/'+ t + '.prt'
            shutil.move( srce , dest)
            #copy errorfile
            srce = 'Errfile'
            dest = self.dirs['err'] +'/'+ t + '.err'
            shutil.move( srce , dest)
            #change directory back
            os.chdir( cwd )
            #clean remaining files
            swan.cleanSwanFiles('./')
            #raise exception
            self.log('Error during swan simulation')
            raise Exception('Error trying to run the swan model')
            #
        #endif
        

        outputlist = copy.deepcopy(self.outputlist)
        name       = 'points'
        outputlist.append( swan.output( \
            self.spotters.spotters , name , ['XP','YP','HS','TM01','DIR'] ) )        
        #
        # Copy result files to appropriate directories...
        #
        for outputfile in outputlist:
            #
            # Source file
            #
            srce = outputfile.filename
            #
            # Point or block data?
            #
            if type(outputfile.outputloc) is list:
                #
                # TABLES
                #--------
                # If point data, load the data from the swan table files
                self.data = swan.readSwanTable(srce)
                #
                # Append the data to output files (if exist)
                ids = self.spotters.getIds()
                #
                for key in self.data:
                    #
                    filename = self.dirs['outputpointsbulk'] \
                      + '/' + key + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , self.data[key] , ids , '{}' )
                    #                    
                #
                # ... transmit point data...
                #
                
                #
                # delete file
                #
                os.remove(srce)
                #
                # SPECTRAL
                #----------
                srce = outputfile.specfilename
                #load as swan spectrum
                spec =  swan.swanReadSpec(srce)

                E = spec.Ef().E
                [a1,b1,a2,b2] = spec.directionalMoments()
                header = [ str(f) for f in spec.f]
                #
                # Write output spectral quantities
                #
                #
                if self.spotters.fake:
                    #
                    # Update fake spotter energies
                    #
                    spec.interpFreq(self.spotters.f)
                    self.spotters.fakeData = spec.Ef()
                    self.spotters.update()
                    #
                #endif
                #
                spc1d = spec.Ef()
                self.data['peakPeriod'] = spc1d.Tp()
                self.data['peakDir'] = spc1d.peakDir()*180/np.pi
                self.data['peakSprd'] = spc1d.peakSprd()*180/np.pi
                    
                for index,pointid in enumerate(ids):
                    #
                    #                   
                    #write spectra
                    filename = self.dirs['outputpointsspec'] \
                      + '/spec_' + pointid + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , E[index,:] , header , '{}' )

                    #write directional moments (a1)
                    filename = self.dirs['outputpointsspec'] \
                      + '/a1_' + pointid + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , a1[index,:] , header , '{}' )
                    
                    #write directional moments (b1)                    
                    filename = self.dirs['outputpointsspec'] \
                      + '/b1_' + pointid + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , b1[index,:] , header , '{}' )
                    
                    #write directional moments (a2)
                    filename = self.dirs['outputpointsspec'] \
                      + '/a2_' + pointid + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , a2[index,:] , header , '{}' )

                    #write directional moments (b2)                    
                    filename = self.dirs['outputpointsspec'] \
                      + '/b2_' + pointid + '_' + self.name + '.sim'
                    utils.writePoints(  self.simulationtime , filename , b2[index,:] , header , '{}' )
                    #
                #
            else:
                #
                # Block data we load...
                self.blockdata = swan.readSwanBlock( srce )

                # If currents are present push those as well
                #if self.currents is not None:
                    #
                #    self.blockdata['U'] = self.currents['U']
                #cat     self.blockdata['V'] = self.currents['V']
                    #
                # ... transmit...
                
                # ... archive
                dest = self.dirs['outputblock'] + '/' + t + outputfile.filename                
                shutil.move( srce , dest)
                #
            #endif
            #
        #endfor
        #
        # Input files...
        #
        srce = 'INPUT'
        dest = self.dirs['inp'] +'/'+ t + self.commandfilename
        shutil.move( srce , dest)        

        for boun in self.bounlist:
            #
            # Copy the boundary files
            #
            if boun.kind.lower()=='file':
                #
                srce = boun.filename

                #load as swan spectrum
                spec =  swan.swanReadSpec(srce)
                spec = spec.Ef()
                header = [ str(f) for f in spec.f]
                E = spec.E.tolist()
                dest = self.dirs['bnd1d'] +'/spec' + boun.filename

                #write 1d Spectra
                utils.writePoints(  self.simulationtime , dest , E , header , '{}' )

                #write Tp, Hm01
                dest = self.dirs['bnd1d'] +'/Tm01' + boun.filename
                utils.writePoints(  self.simulationtime , dest , [spec.Tm0n()] , ['Tm01'] , '{}' )
                dest = self.dirs['bnd1d'] +'/Hm0' + boun.filename
                utils.writePoints(  self.simulationtime , dest , [spec.Hm0()] , ['Hm0'] , '{}' )

                #Move 2d spectrum
                dest = self.dirs['bnd2d'] +'/'+ t + boun.filename
                shutil.move( srce , dest)
                #
            #
        
        #
        # cleanup and change working directory back to old
        #
        swan.cleanSwanFiles('./')
        os.chdir( cwd )
        #
        # Call any nested runs (if applicable)
        #
        if  (len(self.nestedrunlist) > 0 ):
            #
            for nest in self.nestedrunlist:
                #
                nest.run()
                #
            #end for
            #
        #end if nest
        #
        # write run timestamp to file to show succes, also include duration
        #
        self.runend = calendar.timegm( time.gmtime() )        
        if self.parentrun == None:
            #
            dest = self.dirs['log'] + '/' + 'runlist.txt'
            header = ['Start','End','Duration']
            timestamp    = [ time.strftime("%Y%m%d%H%M", time.gmtime( self.runstart ) ), \
                             time.strftime("%Y%m%d%H%M", time.gmtime( self.runend   ) ), \
                             self.runend - self.runstart ]
            utils.writePoints(  self.simulationtime , dest , timestamp , header , '{}' )
            #
    #end def run
    #

    def cleanup(self):
        #
        print('nothing to see here, move along')
        #
    #end def
                
    #
    def addNest( self , nestedrun ):
        #
        # Add a nested simulation to the current run
        #
        import swan
        # 
        # Advance the nest depth
        nestedrun.nestdepth  = self.nestdepth + 1
        
        # Advence the nestnumber (i.e. what is the running number of this grid)
        nestedrun.nestnum    = len( self.nestedrunlist ) + 1
        
        # Set the parentgridname of the nested grid
        nestedrun.parentrun  = self.name

        nestedrun.bounlist.append( swan.tswanboun( kind='NEST' , par=nestedrun.name ) )
        
        # Set the parentgridname of the nested grid
        nestedrun.cgrid.name = nestedrun.name
                
        # Add the nested cgrid to the nested cgrids list
        self.nestedgridlist.append( nestedrun.cgrid )
        
        # Add the nested run to the nested run list
        self.nestedrunlist.append( nestedrun )      
        #
    #end def addNest
    #
    #

    def runRayTracer(self):
        #
        # run the forward raytracer
        #
        import utils
        import copy
        if (self.rt):
            #
            # RUN RAYTRACER
            #
            print('Forward prediction using raytracer')
            [xp, yp] = self.spotters.getLoc()
            spec = copy.deepcopy( self.bounlist[0].spec )
            spec.tof()
            spec.interpFreq( self.spotters.f )
            self.pred= self.raytracer.prediction( spec , xp , yp , 100  )
            
            #
            # PROCESS OUTPUT
            #
            #
            # Append the data to output files (if exist)
            ids = self.spotters.getIds()
            #            
            filename = self.dirs['outputpointsbulk'] \
                + '/' + 'Hsig' + '.ray'
            utils.writePoints(  self.simulationtime , filename , self.pred.Hm0().tolist() ,ids , '{}' )
            
            filename = self.dirs['outputpointsbulk'] \
                + '/' + 'Tm01' + '.ray'
            utils.writePoints(  self.simulationtime , filename , self.pred.Tm0n().tolist() ,ids , '{}' )

            filename = self.dirs['outputpointsbulk'] \
                + '/' + 'Dir' + '.ray'
            utils.writePoints(  self.simulationtime , filename , self.pred.bulkDir().tolist() ,ids , '{}' )


            E = self.pred.Ef().E
            [a1,b1,a2,b2] = self.pred.directionalMoments()
            header = [ str(f) for f in self.pred.f]
            #
            # Write output spectral quantities
            #
            #                    
            for index,pointid in enumerate(ids):
                #
                #                   
                #write spectra
                filename = self.dirs['outputpointsspec'] \
                  + '/spec_' + pointid + '_' + self.name + '.ray'
                utils.writePoints(  self.simulationtime , filename , E[index,:] , header , '{}' )

                #write directional moments (a1)
                filename = self.dirs['outputpointsspec'] \
                  + '/a1_' + pointid + '_' + self.name + '.ray'
                utils.writePoints(  self.simulationtime , filename , a1[index,:] , header , '{}' )
                    
                #write directional moments (b1)                    
                filename = self.dirs['outputpointsspec'] \
                  + '/b1_' + pointid + '_' + self.name + '.ray'
                utils.writePoints(  self.simulationtime , filename , b1[index,:] , header , '{}' )
                    
                #write directional moments (a2)
                filename = self.dirs['outputpointsspec'] \
                  + '/a2_' + pointid + '_' + self.name + '.ray'
                utils.writePoints(  self.simulationtime , filename , a2[index,:] , header , '{}' )

                #write directional moments (b2)                    
                filename = self.dirs['outputpointsspec'] \
                  + '/b2_' + pointid + '_' + self.name + '.ray'
                utils.writePoints(  self.simulationtime , filename , b2[index,:] , header , '{}' )
                #            
            #
        #endif
        #
    #enddef
        
    def log(self,message):
        #
        import time
        #
        timestamp    = time.strftime("%Y%m%d%H%M", time.gmtime( self.simulationtime ) )
        with open( self.dirs['log'] + '/' + 'log.txt','a' ) as fid:
            #
            fid.write(timestamp + ' ' + message + '\n')
            #
        #
    #
    
    def pushToServer(self,asim=None):
        #
        # Push data to backend
        #
        import backend
        #
        if asim is None:
            #
            asim = self.dataAssimilation
            #
            
        if asim:
            #
            modelSource = 'assimilated'
            #
        else:
            #
            modelSource = 'firstGuess'
            #
        #endif
        #
        
        results = backend.pushToServer( self.simulationtime , self.blockdata,
                                        self.spotterlist    , self.data,
                                        self.name           , modelSource,
                                        self.contourlines   , token=self.token)

        if not results:
            #
            message = 'Pushing data to server failed'
            self.log(message)
            print(message)
            #
        #
        #endif
        #
    #enddef
    #
#CLASS
#
