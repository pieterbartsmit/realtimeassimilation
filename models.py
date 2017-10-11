def hmb(dorun=False):
    #
    # This function sets up the hmb model
    #
    import swanrun as sr
    import swan
    import time
    import calendar
    import spotter
    #

    #
    # Coarse Model grid
    #
    cgrid        = swan.compgrid()
    cgrid.domain = [ -122.81 , 37.2 , 0.4 , 0.8 ]
    cgrid.mxc    = 80
    cgrid.myc    = 160
    cgrid.mdc    = 36
    cgrid.flow   = 0.03
    cgrid.fhigh  = 1.

    #
    # Fine Model grid
    #
    cgridf        = swan.compgrid()
    cgridf.domain = [ -122.6 , 37.3 , 0.1 , 0.4 ]
    cgridf.mxc    = 80
    cgridf.myc    = 320
    cgridf.mdc    = 36
    cgridf.flow   = 0.03
    cgridf.fhigh  = cgrid.fhigh
    
    #
    # Model sources
    #
    source             = dict()
    source['tide']     = ['N/A'          , 'N/A'  ]
    source['wind']     = ['ndbc'         , '46012']
    source['dirspec']  = ['ww3'          , '46012']
    #
    # Spotters
    #
    spotterIDs =  ['SPOT-0017','SPOT-0018']
    spotters = spotter.tSpotterlist(spotterIDs)
    #
    # at which boundaries to enforce the spectrum?
    #
    bounlist = ['W','S']
    #
    # output
    #
    output = swan.output(cgrid,  'coarse' , ['HS','TM01','DIR','DSPR','XP','YP','DEP']  )
    outputlist = [output]
    outputf = swan.output(cgridf,  'fine' , ['HS','TM01','DIR','DSPR','XP','YP','DEP']  )
    outputlistf = [outputf]    
    #
    # Create run objects
    #
    run   = sr.swanrun('halfmoonbay', 'hmb', cgrid  , source, outputlist ,spotters,bounlist )
    nest1 = sr.swanrun('fine'  , 'hmb', cgridf , source, outputlistf,spotters )    

    #
    # Get source data and create swan input files
    #
    run.addNest( nest1 )

    if dorun:
        #
        run.setup()    
        run.run()
        #
    #endif
    #
    return( run )
    #

def innershelf(assimilation, raytracer=None, spotters=None):
    #
    # This function sets up the hmb 
    # 
    import swanrun as sr
    import swan
    import numpy as np
    import spotter
    #
    # Coarse Model grid
    #
    cgrid = swan.compgrid()
    cgrid.domain =  [ -121.019 , 34.5,     0.5 , 1 ]
    cgrid.mxc    = 200
    cgrid.myc    = 400
    cgrid.mdc    = 36
    cgrid.flow   = 0.03
    cgrid.fhigh  = 0.5

    ocgrid = swan.compgrid()
    ocgrid.domain =  [-120.9690 , 34.5,     0.4 , 0.8 ]
    ocgrid.mxc    = 100
    ocgrid.myc    = 200
    ocgrid.mdc    = 36
    ocgrid.flow   = 0.03
    ocgrid.fhigh  = 0.5
    
    #
    # Fine grid model
    #    
    source             = dict()
    #source( Field) = [ Online source , stationid ]
    source['tide']     = ['coops'        , '9412110']
    source['wind']     = ['ndbc'         , '46011']
    source['dirspec']  = ['ww3'          , '46011']

    #
    # Spotters
    #

    
    if spotters == None:
        spotterIDs =  []
        #for num in range(0,18):
            #
            #spotterIDs.append( 'Fake'+str(num) )
            #
        #endif

        spotterIDs = [
            'SPOT-0009',
            'SPOT-0010',
            'SPOT-0011',
            'SPOT-0012',
            'SPOT-0014',
            'SPOT-0015',
            'SPOT-0019',
            'SPOT-0020',
            'SPOT-0021',
            'SPOT-0022',
            'SPOT-0023',
            'SPOT-0024',
            'SPOT-0025',
            'SPOT-0026',            
            'SPOT-0028',
            'SPOT-0029',
            'SPOT-0030',
            'SPOT-0031' ]
            
        lat = [35.13438, 35.19234, 35.17114, 35.06744, 34.81299, 35.10238, 34.95907,
                   34.99926, 34.76484, 34.81584, 35.06465, 35.12991, 34.86176, 34.9084, 35.027, 34.85364, 34.91913, 34.95188]
        lon = [-120.80777, -120.86745, -120.90178, -120.82498, -120.79141, -120.7461, -120.79741, -120.66176,
                   -120.67688, -120.6991, -120.65216, -120.67119, -120.71636, -120.72218, -120.72375, -120.63095, -120.68178, -120.72712]


        spotters = spotter.tSpotterlist(spotterIDs,lat=lat,lon=lon )
        #
    #endif
    #
    # at which boundaries to enforce the spectrum?
    #
    bounlist = ['W','S']
    #
    # output
    #
    output = swan.output(ocgrid,  'coarse' , ['HS','TM01','DIR','DSPR','XP','YP','DEP']  )
    outputlist = [output]
    
    #
    # Create run object
    #('coarse', 'hmb', cgrid  , source, outputlist ,spotters,bounlist )
    #
    run = sr.swanrun( 'vandenberg','vandenberg', cgrid , source, outputlist,spotters,bounlist,dataAssimilation = assimilation,
                          rayTracer = raytracer )
    #
    return(run , spotters)
