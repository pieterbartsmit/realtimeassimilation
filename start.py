import time
import models
import spotter
import calendar
import sys
import matplotlib
matplotlib.use('Agg') 
#
interval           = 60

#
# Setup a swanrun object for the region of interest
#
#model              = models.hmb() #hmb model
[model,spotters]   = models.innershelf(False)
[asim,spotters]    = models.innershelf(True , model.raytracer,spotters=spotters )

#lastsimulationtime = 0.
#
lastrun = 0
while True:
    #
    # Check if there are updated spotter observations...
    #
    try:
        #
        updateAvailable = model.update(1800)
        #
    except Exception as e:
        #
        print( 'Could not update spotters' )
        print( str(e) )
        model.log( str(e) )        
        #    

    if updateAvailable:
        #
        # ...if so, retrieve wind/tide/boundary condition data ...
        #
        #
        try:
            #
            # ... and run the model
            #
            print( 'Running the first-Guess Model set!' )
            model.setup(  )
            model.runRayTracer()
            model.run()

            # synchronize the models
            asim.simulationtime = model.simulationtime
            asim.spotters       = model.spotters

            #model.pushToServer()
            #model.pushToServer(True)
            print( 'Running the Assimilated model!' )
            asim.setup(  )
            asim.runRayTracer()
            asim.run()

            lastrun =  calendar.timegm ( time.gmtime() ) 
            #
            #
            #
        except Exception as e:
            #
            # Many reasons the model could fail - in particular online data
            # sources are sometimes not available, this keeps the model going
            # until they do again - clutch solution, needs to be investigated
            #
            print( 'model failure' )
            print( str(e) )
            model.log( str(e) )
            #
        #
    else:
        #
        # ... if no data available, print interval since last simulation ...
        #
        string = str( calendar.timegm ( time.gmtime() ) - lastrun ) + ' seconds since last model run'
        print( string )
        #
    #end if
    #
    # ... and sleep for the desired interval
    #
    time.sleep( interval )
    #
#end while
#
