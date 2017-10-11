def points( kind , directory , asim=True):
    #
    import swanrun
    import utils
    import matplotlib.pyplot as plt
    
    plt.close('all')
    if kind.lower() == 'hs':
        #
        files = swanrun.getOutputFilenames( 'hspswa' , directory,gridname=None,dataAssimilation=asim )
        model = utils.loadPoints( files )
        files = swanrun.getOutputFilenames( 'hspray' , directory,gridname=None,dataAssimilation=asim )        
        ray   = utils.loadPoints( files )
        files = swanrun.getOutputFilenames( 'hspobs' , directory,gridname=None,dataAssimilation=False )        
        obs   = utils.loadPoints( files )
    #

    plt.figure(1 )
    plt.plot( model )
    plt.plot( ray )
    #plt.plot( obs)    
    plt.show()



