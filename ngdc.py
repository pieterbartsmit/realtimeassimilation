def getCoastalReliefData( rect , resolution = [3,3] , filename = 'bat.dep' ):
    #
    import code
    import numpy as np
    import urllib.request        
    import pandas as pd

    if type(rect) in [list,tuple]:
        #
        rect = np.array( rect )
        #       
    #
    # Convert resolution from arcseconds to degrees
    # - check if numpy array, else convert
    #
    if type(resolution) in [list,tuple]:
        #
        resolution = np.array( resolution )
        #
        
    resolution = resolution / 3600.

    #The API is basically reverse engineered form going through the interface at:
    #    https://maps.ngdc.noaa.gov/viewers/wcs-client/
    #and copying the link to download the data. 
    #
    url = 'https://maps.ngdc.noaa.gov/mapviewer-support/wcs-proxy/wcs.groovy?' \
        + 'filename=crm.xyz&request=getcoverage&version=1.0.0&service=wcs&coverage=crm&CRS=EPSG:4326&format=xyz'

    #Base URL
    #Longt. res. (dec. degrees)
    #Lat.   res. (dec. degrees)
    #Set up the geographic bounding box
    #West  boundary latitude  (sign. dec. degrees)
    #South boundary longitude (sign. dec. degrees)
    #East  boundary latitude  (sign. dec. degrees)
    #North boundary longitude (sign. dec. degrees)        
    url = url \
        + '&resx=' + '{:10.8f}'.format( resolution[0] )           \
        + '&resy=' + '{:10.8f}'.format( resolution[1] )           \
        + '&bbox='                                                \
        + '{:20.14f}'.format( rect[0]         ).strip() + ','     \
        + '{:20.14f}'.format( rect[1]         ).strip() + ','     \
        + '{:20.14f}'.format( rect[0]+rect[2] ).strip() + ','     \
        + '{:20.14f}'.format( rect[1]+rect[3] ).strip()          

    #
    # Get the data and store on local disk at filename
    #
   
    [filename , header ] = urllib.request.urlretrieve(url, filename)

    #
    # Load data and do clever stuff
    #
    data = pd.read_csv( filename,engine='python',sep=' ',header=None)
    
    lon = data[0].values
    lat = data[1].values
    dep = -1.0 * data[2].values

    nlon = int( round( rect[2] / resolution[0] , 0) + 1 )
    nlat = int( round( rect[3] / resolution[1] , 0) + 1 )

    n  = len(lon)

    if n == nlon*nlat:
        #        
        dep = ( np.reshape( dep , ( nlat , nlon ) ) )
        lon = ( np.reshape( lon , ( nlat , nlon ) ) )
        lat = ( np.reshape( lat , ( nlat , nlon ) ) )

        if ( lon[ 1,2 ] - lon[ 1,1 ] < 0. ):
            #
            # Make sure that arrays are in ascending lon. order
            #
            lat = np.flip(lat,1)
            lon = np.flip(lon,1)
            dep = np.flip(dep,1)            
        
        if ( lat[ 2,1 ] - lat[ 1,1 ] < 0. ):
            #
            # Make sure that arrays are in ascending lat. order
            #
            lat = np.flip(lat,0)
            lon = np.flip(lon,0)
            dep = np.flip(dep,0)
        
        return( {'lat':lat,'lon':lon,'dep':dep } )            
        #
    else:
        #
        #code.interact( local=locals() )
        raise Exception('Dimensions of the grid are not correct')
        #
    
