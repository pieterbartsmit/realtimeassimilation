#
# Python package to retrieve wavewatch files
#
import code
import numpy as np
import spectral

def getSpectrum( epochtime , stationID , workdir='./' ):
    #
    import spectral
    res = getLatestDirectionalSpectrumFromServer(date=epochtime,buoynum=stationID, workdir=workdir )

    spec = res['spec'].interpLoc( epochtime )
    
    return( spec )
    #
    #
    
def getLatestDirectionalSpectrumFromServer( date='',buoynum=0, workdir='./'  ):
    #
    # Dependencies
    #
    import urllib.request
    import time
    import re
    import calendar
    #---------------------------------------------------------------------------    
    #
    # This routine grabs spectral buoy data from the NOAA nomads server. Input is a 
    # date and a buoynum. Note that only the last few days are stored on the NOMAD
    # server.
    #
    # In no date or buoy are given, it defaults to "today (in UTC!)" and buoy 46011 - which just
    # happens to be the buoy for the ONR Innershelf experiment.
    #
    # PBS - July, 2017
    #
    # extended to grab from data.noda.noaa.gov/pub/data for historical dates
    #
    # PBS - Apr, 2018
    #---------------------------------------------------------------------------
    #
    # Code description:
    # 1) sanity checking of input
    # 2) grab data from server
    # 3) reformat data to usuable format
    # 4) return frequency/direction spectrum
    #

    def getUrl( epoch ):
        #
        import calendar
        import time
        import os

        timeStruc = time.gmtime( epoch )
        yy  = timeStruc.tm_year
        mm  = timeStruc.tm_mon
        dd  = timeStruc.tm_mday
        hh  = timeStruc.tm_hour
        min = timeStruc.tm_min
        sec = timeStruc.tm_sec

        time.strftime("%Y%m%d")
        isFile = False
        epochnow = calendar.timegm( time.gmtime() )
        #https://data.nodc.noaa.gov/ncep/nww3/2017/09/points/multi_1_base.buoys_spec.201709/multi_1.46011.SPEC.201709
        #https://data.nodc.noaa.gov/ncep/nww3/2017/09/points/multi_1_base.buoys_spec.201709/multi_1.46011.spec201709
        #ftp://polar.ncep.noaa.gov/pub/history/waves/multi_1/201710/points/multi_1_base.buoys_part.201710.tar.gz        
        if ( epoch - epochnow > -3600*24*3 ):
            day = time.strftime("%Y%m%d" , time.gmtime( epoch ) )
            url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/wave/prod/multi_1.' \
              + day + '/bulls.t00z/multi_1.' + buoynum + '.spec'
            realtime = True
        else:
            #
            # If requested epoch is the first spectrum of the month, we need to retrieve
            # the previous months data, as it is stored there by wavewatch for some ^%&$%^$ reason
            #

            if ( dd==1 ) and (hh==0 or hh==1):
                #
                print( 'doin it')
                mm = mm - 1
                if mm < 1:
                    #
                    mm = 12
                #

            yy = str( yy )
            if mm < 10:
                mm = '0'+str(mm)
            else:
                mm = str(mm)

            fileexists = os.path.isfile(workdir+'/'+buoynum+'.SPEC.'+yy + mm)

            if  fileexists:
                isFile = True
                url = workdir+'/'+buoynum+'.SPEC.'+yy + mm
            else:                
                url = 'https://data.nodc.noaa.gov/ncep/nww3/' + yy +'/' + mm + \
                  '/points/multi_1_base.buoys_spec.' + yy + mm + '/' +         \
                  'multi_1.' + buoynum + '.SPEC.' + yy + mm

            #url = 'ftp://polar.ncep.noaa.gov/pub/history/waves/multi_1/'+yy+mm+ \
            #  '/points/multi_1_base.buoys_spec.' + yy + mm 201710.tar.gz
            realtime = False
        #
        return( url, realtime, isFile )
    #ENDFUNCTIONDEF
    
    #---------------------------------------------------------------------------
    # 1) sanity check
    #
    if date=='':
        #
        day  = 3600 * 24
        date = [ calendar.timegm( time.gmtime() ),
                 calendar.timegm( time.gmtime() ) -   day ,
                 calendar.timegm( time.gmtime() ) - 2*day ]
        #
    else:
        #
        if type(date) != list:
            #
            # 
            #
            date = [ date ]
            #
        #

        
    #fi
        
    if buoynum==0:
        #
        buoynum = '46011'
        #
    #fi


    
    #---------------------------------------------------------------------------
    # 2) grab data from server
    #
    #
    # URL - hardcoded at this moment
    #


    # Pick the latest existing date in the dates list:
    found = False
    #https://data.nodc.noaa.gov/ncep/nww3/2017/09/points/multi_1_base.buoys_spec.201709/
    for day in date:
        #url = 'http://nomads.ncep.noaa.gov/pub/data/nccf/com/wave/prod/multi_1.' \
        #    + day + '/bulls.t00z/multi_1.' + buoynum + '.spec'
        url, realtime, isFile = getUrl( day )

        if isFile:
            response = open(url,'rt')
            found = True
            print('we are in business')
            break
        
        try:
            #
            # Check if the file exists on server
            #
            response = urllib.request.urlopen( url )
            found = True
            break
            #    
            #
        except:
            pass
        #endtry
        #
    #endfor

    if not found:
        #
        return ( {'exist':False } )
        #
    #fi    

    #        
    #
    # Read all data
    #
    lines = ( response.read() )
    #
    # Decode to ascii
    #
    if not isFile:
        lines = lines.decode('utf-8')
    #
    # Remove newline characters
    #
    lines = lines.replace('\n',' ')
    #
    # Use regular expression to split the string
    #        
    lines = [p for p in re.split("( |\\\".*?\\\"|'.*?')", lines) if p.strip()]
    #
    # Now we have a flat list - 
    #
    numfreq = int( lines[1] )
    numdir  = int( lines[2] )
    numpoint= int( lines[3] )

    freq = np.array( list( map( float, lines[ 5:5+numfreq]) ) )
    ang  = np.array( list( map( float, lines[ 5+numfreq : 5 + numfreq + numdir ] ) ) )
    
    ang = 0.5 * np.pi - ang # + np.pi
    #ang =  ang  + np.pi
    for ii in range(0,len(ang) ):
        #
        if ang[ii] > np.pi:
            ang[ii] = ang[ii] - 2.0*np.pi
        if ang[ii] < -np.pi:
            ang[ii] = ang[ii] + 2.0*np.pi
    
    lines = lines[ 5 + numfreq + numdir :   ]
    dataFrameSize  = numfreq * numdir + 9
    numberOfFrames = len( lines ) // dataFrameSize

    dates = lines[ 0 : -1 : dataFrameSize]
    times = lines[ 1 : -1 : dataFrameSize]
    Hs    = np.array(lines[ 5 : -1 : dataFrameSize])
    Dm    = np.array(lines[ 6 : -1 : dataFrameSize])

    datestimes = [ m + n for m,n in zip(dates,times) ]

    times = [ calendar.timegm( time.strptime( x , '%Y%m%d%H%M%S' ) ) for x in datestimes ]
    times = np.array(times)
    #
    indices = np.argsort( ang )
    ang = ang[indices]
    E = np.zeros( (numberOfFrames , numdir , numfreq  ) )
    #        
    for i in range( 1 , numberOfFrames ):
        #
        istart = (i-1)*dataFrameSize +9
        iend   = i*dataFrameSize
        dat = np.array( lines[istart:iend] )
        dat = np.reshape( dat , ( numdir , numfreq ) )
        dat = dat[indices , :]
        E[ i , : , :  ] = dat
        #
    #


    ang = ang * 180./np.pi
    E   = E * np.pi/180.



    if not realtime:
        #
        # First spectrum returned is full of zeros for some reason (stored in previous months data)
        # This potentially screws up interpolation. Here we set the first entry the same as
        # the second to ensure nearest neighbour interpolation
        E[0,:,:] = E[1,:,:]

    spec = spectral.spectrum( {'E':E , 'f':freq , 'ang':ang, 'loc':times } )
    spec.regularize(36)
        
        
    #plt.pcolor(dir,freq,dat)
    #plt.pcolor(dat)
    #plt.show
    #code.interact( local=locals() )
    return ( {'spec':spec , 'exist':True, 'Hs':Hs, 'Dm':Dm } )
