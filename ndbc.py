baseurl = 'http://www.ndbc.noaa.gov/data/realtime2/'

def getMeanWindData( epochtime , length=1800 , stationID='46011',workingdirectory='./'  ):
    #
    #
    #
    import numpy as np
    import code
    import pandas as pd

    fname = workingdirectory+'ndbcwind.csv'
    data = getWindData( stationID , filename=fname)
    
    for jj in range(0,4):
        #
        # Find the closest hour with data
        #
        epoch = epochtime - 3600 * jj
        msk      = ( data['date'] >= epoch - length//2 ) & \
                   ( data['date'] <= epoch + length//2 )
        if np.sum(msk) > 0:
            break
        #
        #endif
        #
    #endfor
    ux       = (np.cos( ( 270. - data[msk]['WDIR'] ) * np.pi/180. ) * data[msk]['WSPD']).mean()
    uy       = (np.sin( ( 270. - data[msk]['WDIR'] ) * np.pi/180. ) * data[msk]['WSPD']).mean()

    Udir     = np.arctan2( uy , ux ) * 180. / np.pi 
    U        =( ux**2 + uy**2 ) ** 0.5
    return( {'U':U,'Udir':Udir} )
    
def getWindData( stationID='46011' , filename='ndbcwind.csv' ):
    #
    import code
    import numpy as np
    import urllib.request        
    import pandas as pd
    import calendar
    import time
    from datetime import datetime
    from pytz import timezone    
    #
    global baseurl
    url = baseurl + stationID + '.cwind'

    [filename , header ] = urllib.request.urlretrieve(url, filename)
    data = pd.read_csv( filename,engine='python', delim_whitespace=True,skiprows=[1],header=0 )
    dates = list()
    #
    for ind,row in data.iterrows():
        #
        date = datetime( int(row['#YY']) , int(row['MM']) , int(row['DD']), \
                         int(row['hh']) , int(row['mm']) ,0,0, tzinfo=timezone('UTC'))
        date = calendar.timegm( date.utctimetuple() )
        dates.append( date )
        #
    #
    data['date'] = dates
    return(data)
#    time.time
    
def getSpectralData( url , workingdirectory = None,kind=0 ):
    import code
    import numpy as np
    import urllib.request        
    import pandas as pd
    import calendar
    import time
    import os
    from datetime import datetime
    from pytz import timezone
    import spectral
    
    if workingdirectory == None:
        #
        workingdirectory = '.' + os.path.sep
        #
    #
    filename = workingdirectory + 'temp.csv'
    [filename , header ] = urllib.request.urlretrieve(url, filename)
    data = pd.read_csv( filename,engine='python', delim_whitespace=True,skiprows=[0],header=None )
    dates = list()
    
    E    = data.loc[ : , 6-kind : data.shape[1] : 2 ]
    freq = data.loc[ 0 , 7-kind : data.shape[1] : 2 ]
    E=E.values
    freq =freq.apply( lambda x: x.replace('(',''))
    freq = freq.apply( lambda x: x.replace(')',''))
    freq = freq.apply( lambda x: float(x) )
    freq = freq.values

    for ind,row in data.iterrows():
        #
        date = datetime( int(row[0]) , int(row[1]) , int(row[2]), \
                         int(row[3]) , int(row[4]) ,0,0, tzinfo=timezone('UTC'))
        date = calendar.timegm( date.utctimetuple() )
        dates.append( date )
        #

    #
    return( E , freq , dates )


def getSpec( stationID='46011' , workingdirectory = None, epochtime=None ):
    import code
    import numpy as np
    import urllib.request        
    import pandas as pd
    import calendar
    import time
    import os
    from datetime import datetime
    from pytz import timezone
    import spectral

    global baseurl


    files = ['.data_spec' , '.swdir','.swdir2','.swr1','.swr2']
    kind  = [0,1,1,1,1]
    data = []
    for index,file in enumerate(files):
        #
        url = baseurl + stationID + file        
        [tmp,freq,dates] = getSpectralData( url , workingdirectory = workingdirectory,kind=kind[index] )
        data.append(tmp)
        #

    #
    # Convert to lhs moments - see http://www.ndbc.noaa.gov/measdes.shtml
    #
    angle1 = ( 270 -   data[1] ) * np.pi / 180.
    angle2 = ( 540 - 2*data[2] ) * np.pi / 180.

    msk = data[3] == 999.0
    msk = data[4] == 999.0

    data[3][msk] = np.nan
    data[4][msk] = np.nan    

    #D(f,A) = (1/PI)*(0.5+R1*COS(A-ALPHA1)+R2*COS(2*(A-ALPHA2))).
    
    R1     =  data[3]
    R2     =  data[4]

    E = data[0]
    a1 = R1 * np.cos(angle1)
    b1 = R1 * np.sin(angle1)
    a2 = R2 * np.cos(angle2)
    b2 = R2 * np.sin(angle2)
    
    spec = spectral.spectrum1d( {'E':E , 'f':freq, 'loc':dates,'a1':a1,'b1':b1,'a2':a2,'b2':b2 } )
    #data['date'] = dates

    if epochtime is not None:
        #
        if epochtime < 0:
            #
            epochtime = dates[0]
            #
        #
        spec = spec.interpLoc( epochtime )
    else:
        #
        epochtime = dates
        #
    #
    
    lat,lon = getStationInfo( stationID,workingDirectory=workingdirectory )
    return(spec,lat,lon,epochtime)

def getLatestSpec(stationID='46012',workingDirectory=None):
    #
    spec, lat,lon, epochtime = getSpec( stationID=stationID , workingdirectory = workingDirectory, epochtime=-1 )
    return( spec,lat,lon,epochtime)
    #

def getStationInfo( stationID,workingDirectory=None ):
    #
    url = 'http://www.ndbc.noaa.gov/data/stations/station_table.txt'
    import urllib.request        
    import pandas as pd
    import os

    if workingDirectory == None:
        #
        workingDirectory = '.' + os.path.sep
        #    

    filename = workingDirectory + 'ndbc_temp_station_loc.csv'
    [filename , header ] = urllib.request.urlretrieve(url, filename)
    data = pd.read_csv( filename,engine='python', sep='|',skiprows=[0,1],header=None )

    #Set the first column as index (contains buoy ID's)
    data=data.set_index(0)

    #Get the correct buoy, and the 6th index contains lat/lon string
    data = data.loc[stationID][6].split()
    
    lat = float( data[0] ) if data[1]=='N' else -float( data[0] )
    lon = float( data[2] ) if data[1]=='E' else -float( data[2] )
    return( lat, lon)
