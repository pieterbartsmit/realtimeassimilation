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
model              = models.innershelf() #hmb model

#lastsimulationtime = 0.
#
while True:
    #
    # Check if there are updated spotter observations...
    #
    try:
        #
        updateAvailable = model.update()
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
            print( 'Simulation Start!' )
            model.setup(  )
            model.runRayTracer()
            model.run()
            #model.pushToServer()
            #model.pushToServer(True)
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
        string = str( calendar.timegm ( time.gmtime() ) - model.simulationtime ) + ' seconds since last update'
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
