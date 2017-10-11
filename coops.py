#
# Python package to retrieve data from coops
# API documented at https://tidesandcurrents.noaa.gov/api/
#
import code
import numpy as np

def getTideData( epochtime, station='',workingdirectory='./'):
    #
    import time

    begin_date = time.strftime("%Y%m%d" , time.gmtime(epochtime) )
    end_date    = begin_date
    
    data = getData( station, kind=0,begin_date='',end_date='',workingdirectory=workingdirectory)

    it = (np.abs(data['t']-epochtime)).argmin()
    
    return( {'z':data['z'][it] })
    #

def getData( station, kind=0,begin_date='',end_date='',datum='MSL',time_zone='GMT',units='metric',workingdirectory='./'):
    #
    # Grab wind (kind=1) or tide (kind=0) data from the coops servers for a given stationid. 
    #
    import urllib.request
    import time
    import calendar
    import re
    import pandas as pd    
        
    if kind == 0:
        #
        product = 'water_level'
        #
    elif kind == 1:
        #
        product = 'wind'            
        #
    #endif
    #
    if station=='':
        #
        station='9412110'
        #
    #fi
            
    if begin_date=='':
        #
        begin_date = time.strftime("%Y%m%d" , time.gmtime() )
        #
    #fi

    if end_date=='':
        #
        end_date = time.strftime("%Y%m%d" , time.gmtime() )
        #
    #fi  
        
    #
    # Generate the url
    #
    url = 'https://tidesandcurrents.noaa.gov/api/datagetter?'    
    tideInfo = 'product='     + product    + '&application=NOS.COOPS.TAC.WL' \
             + '&begin_date=' + begin_date + '&end_date='+end_date  \
             + '&datum='      + datum      + '&station=' +station  \
             + '&time_zone='  + time_zone  + '&units='   +units    \
             + '&format=csv'
    url = url + tideInfo
    #
    # Get Data
    #
    #
    # Note on OSX python 3.6 this can throw an error due to missing certificates:
    #      run Applications/Python 3.6/install Cert...   to fix this
    #
    #
    # Dump to local text file before processing
    #
    [filename , header ] = urllib.request.urlretrieve(url, workingdirectory+'/temp.csv')
    #
    # Process with pandas
    #
    data = pd.read_csv( filename,sep='\s*,\s*',engine='python' )
    #
    t    = data['Date Time']
    t = [ calendar.timegm( time.strptime( x  , "%Y-%m-%d %H:%M"))  for x in t ]
    t = np.array(t)
    #    
    if kind == 0:
        # We only need time and water level...        
        z    = data['Water Level'].values
        return( {'t':t,'z':z} )
        #
    elif kind == 1:
        #
        #code.interact( local=locals() )
        U    = data['Speed'].values
        Udir = data['Direction'].values        
        return( {'t':t,'U':U,'Udir':Udir} )
        #code.interact( local=locals() )           
        #

