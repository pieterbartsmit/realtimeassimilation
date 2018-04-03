class spotres:
    #
    def __init__( self , data):    
        #
        import time
        import calendar
        import spectral
        import numpy as np

        if not 'spec' in data:
            #
            
            self.Hsig       = data['waves'][0]['significantWaveHeight']
            self.Tp         = data['bulk']['peakPeriod']
            self.Tm01       = data['bulk']['meanPeriod']
            self.peakDir    = data['bulk']['peakDirection']
            self.peakSprd    = data['bulk']['peakDirectionalSpread']        
            self.Dir        = data['bulk']['direction']
            self.dirSprd    = data['bulk']['directionalSpread']
            self.lat        = data['lat']
            self.lon        = data['lon']
            self.id         = data['spotterId']
            self.name       = data['name']
            self.updatetime = [data['timestampStart'],data['timestampEnd']]
            self.epochtimes  = [ calendar.timegm( time.strptime( x , '%Y-%m-%dT%H:%M:%S.%fZ' ) )
                                 for x in self.updatetime]        
            self.epochtime   = int(( self.epochtimes[0] + self.epochtimes[1] ) // 2)        
            self.fullMessage = False
            self.spec = None

        
            if 'frequencyData' in data:
                #
                data=data['frequencyData']
                #
                if 'f' in data:
                    #
                    self.fullMessage = True
                    # adhoc conversion to densities based on full messsage format
                    # hardcoded now - should disappear once full message returns
                    # densities
                    f = data['f']
                    df = np.ones(39)
                    df[0:33]  = 0.009765625
                    df[33:38] = 3 * 0.009765625
                    df[38]    = 0.751953125-0.517578125                
                    E=data['Szz']/df/1000000

                    var = {'f':data['f'],'E':E,'a1':data['a1'],
                            'a2':data['a2'],'b1':data['b1'],'b2':data['b2']}
                    self.spec = spectral.spectrum1d( var )
                    #
                #end if
                #
            #end if
            #
        else:
            #
            
            self.Hsig       = data['spec'].Hm0()
            self.Tp         = data['spec'].Tp()
            self.Tm01       = data['spec'].Tm0n()
            self.peakDir    = 270.-data['spec'].peakDir()*180./np.pi
            if self.peakDir < 0.:
                self.peakDir = self.peakDir + 360.

            if self.peakDir > 360.:
                self.peakDir = self.peakDir - 360.
            
            self.peakSprd    = data['spec'].peakSprd()*180./np.pi
            self.Dir        = 270.-data['spec'].bulkDir()*180./np.pi
            #
            if self.Dir < 0.:
                #
                self.Dir = self.Dir + 360.
                #
            #

            if self.Dir > 360.:
                #
                self.Dir = self.Dir - 360.
                #
            #
            
            self.dirSprd    = data['spec'].bulkSprd()*180./np.pi
            self.lat        = data['lat']
            self.lon        = data['lon']
            self.id         = data['spotterId']
            self.name       = data['name']
            
            self.updatetime = ['0','0']
            self.epochtimes  = [ data['epochtime'],data['epochtime'] ]        
            self.epochtime   = data['epochtime']

            if not np.any( np.isnan( data['spec'].E ) ):
                self.fullMessage = True
            else:
                self.fullMessage = False
                
            self.spec = data['spec']
            #
        #endif
        #
    #ebd def init
                
    def makeNaN(self):
        #
        # Set all data to nan/none
        #
        import numpy as np
        #
        self.Hsig       = np.NaN
        self.Tp         = np.NaN
        self.Tm01       = np.NaN
        self.peakDir    = np.NaN
        self.peakSprd   = np.NaN
        self.Dir        = np.NaN
        self.dirSprd    = np.NaN

        self.fullMessage = False
        self.spec = None        
    #
#end class

class tSpotterlist:
    #
    #
    def __init__ ( self , spotterIDs, fake=None,fakeinterval=60,endpoint=None,
                   fakeData=None,lat=None,lon=None,token=None,
                   isaliveTime=7200,matlabSpotter=False,targetTime=0,
                   matlabSpotterDir=None):
        #
        # Create a spotter list from a set of spotterIDs. Note that
        # when an array of xy coordinates is provided to fake, the
        # spotters are "virtual" spotters.
        #
        import numpy as np

        #PROPERTIES
        self.fakeinterval = fakeinterval
        self.numSpot = 0
        self.spotters = []
        self.numFullSpec = 0
        self.fullSpecIndices = []
        self.available = []
        self.f = np.zeros(39)
        self.isaliveTime = isaliveTime
        import sys
        import platform


        #
        #Are these Spotters 'matlab spotters'; i.e. spectra created from sdcard
        #(Innershelf assimilation)
        self.matlabSpotter=matlabSpotter #Switch to toggle to matlabspotter
        self.targetTime=targetTime       #the targetTime to get data
        if matlabSpotterDir is None:
            #
            if platform.system() == 'Darwin':
                #
                matlabSpotterDir='/Users/pieterbartsmit/Google Drive/Workingmanuscripts/innershelf/dataProcessing/filtered/2share'
                #
            else:
                #
                appdir = '/mnt/2ndHDD/models/innershelf/obs/'
                #            
            
        self.matlabSpotterDir = matlabSpotterDir
        
        if token is not None:
            #
            self.token = token
            #
        else:
            #
            # NOTE: THIS IS THE TOKEN FOR THE ONR INNERSHELF MODEL - only here for backwards
            # compatibility with that model which had token hardcoded. Statement MUST be removed after innershelf
            # experiment ends (Dec 2017) e (PBS Sep 2017)
            #
            self.token = 'a8db94573ddd8ab23a470493091866'
            #
        #

        self.lat=lat
        self.lon=lon

        #Start Initialize
        df         = np.ones(39)
        df[0:33]   = 0.009765625
        df[33:38]  = 3 * 0.009765625
        df[38]     = 0.751953125-0.517578125
        self.f[0] = 0.009765625 * 3
        for i in range(1,39):
            #
            self.f[i] = self.f[i-1] + df[i]/2 + df[i-1]/2
            #
        #
        
        if endpoint==None:
            #
            self.endpoint = 'https://wavefleet.herokuapp.com/api/latest-data'
            #
        else:
            #
            self.endpoint = 'https://wavefleet.herokuapp.com/api/latest-data'
            #
        
        if (not (type(spotterIDs) is list ) ):
            #
            spotterIDs = [spotterIDs]
            #
        #endif
        #
        if not (type(fake)==list or type(fake)==np.ndarray):
            #
            self.fake = False
            self._xy = None
            self.fakeData = None
            #
        else:
            #
            self.fake = True
            self._xy = fake
            self.fakeData = fakeData
            #

        print( 'Fetching Spotter information' )            
        for index, spotterID in enumerate(spotterIDs):
            #
            self.addSpotter(spotterID,index)
            #
        #end for
        #
        #
        # Update the mean and latest update time of
        # all spotters in the list
        #
        print( 'total spotters: ' + str(self.numSpot) + ' , available spotters: ' + str(len(self.available)) )
        print('===============================================================')
        self.updateTime()
        #
    #end __init__
    #

    def getSpec(self):
        #
        # get all the spectra from the spotters as one
        # convinient spectral object. Note that for 
        # spotters without full message nans are returned
        #
        import numpy as np
        import spectral
        
        ind  = self.fullSpecIndices[0]
        nf   = self.spotters[ind].spec.nf
        f    = self.spotters[ind].spec.f
        nloc = len( self.spotters )
        #
        E  = np.zeros( ( self.numSpot , nf ) ) + np.nan
        a1 = np.zeros( ( self.numSpot , nf ) ) + np.nan
        b1 = np.zeros( ( self.numSpot , nf ) ) + np.nan
        a2 = np.zeros( ( self.numSpot , nf ) ) + np.nan
        b2 = np.zeros( ( self.numSpot , nf ) ) + np.nan
        j = 0
        #
        for index,spotter in enumerate( self.spotters ):
            #
            if index in self.fullSpecIndices:
                #
                E[ index  , :]  = spotter.spec.E
                a1[ index , :]  = spotter.spec.a1
                b1[ index , :]  = spotter.spec.b1
                a2[ index , :]  = spotter.spec.a2
                b2[ index , :]  = spotter.spec.b2
                #
            #endif
            #     
        #endfor
        #
        spec = spectral.spectrum1d( {'E':E,'f':f,'a1':a1,'b1':b1,'a2':a2,'b2':b2} )
        return(spec)
            
    def getLoc(self):
        xp = []
        yp = []
        for spotter in self.spotters:
            #
            xp.append( spotter.lon )
            yp.append( spotter.lat )
            #
        #
        return( xp , yp )
    #
    
    #
    def updateTime(self):
        #
        import time
        self.meanTime   = 0
        self.lastUpdate = 0
        n               = len(self.available)
        #
        for index in self.available:
            spotter = self.spotters[ index ]
            #
            #print(spotter.id)
            #print(time.gmtime(spotter.epochtime))
            self.meanTime = self.meanTime + spotter.epochtime/n
            #
            if spotter.epochtime > self.lastUpdate:
                #
                self.lastUpdate = spotter.epochtime
                #
            #endif
            #
        #end for
        #
    #end def
    #
    def addSpotter( self , spotterID,num ):
        #
        #spotter = self.getData( spotterID,endpoint=self.endpoint,num=num )
        if True: #try:
            spotter = self.getData( spotterID,endpoint=self.endpoint,num=num )
        
            self.spotters.append( spotter )
            self.available.append( self.numSpot )
            #
            if spotter.fullMessage:
                #
                self.numFullSpec += 1
                self.fullSpecIndices.append( self.numSpot )
                self.f = spotter.spec.f
                #
            #        
            self.numSpot = self.numSpot + 1
        else: #except Exception as e:
            #
            print( '- could not add spotter ' + spotterID )
            print( str(e) )
            if not self.lat is None:
                #getFakeData( self , spotterID, num ):
                # including a fake spotter instead
                #
                spotter = self.getFakeData( spotterID, num )
                spotter.lat = self.lat[num]
                spotter.lon = self.lon[num]
                spotter.makeNaN()

                self.spotters.append( spotter )
                self.numSpot = self.numSpot + 1
            #
        #
    #enddef

    def getFakeData( self , spotterID, num ):
        #
        import time
        import spectral
        
        data  = dict()
        datab = dict()
        datab['significantWaveHeight'] = 0.
        datab['peakPeriod'] = 0.
        datab['meanPeriod'] = 0.
        datab['direction']  = 0.
        datab['peakDirection']  = 0.
        datab['peakDirectionalSpread']  = 0.        
        datab['directionalSpread'] = 0.
        data['bulk'] = datab
        #
        if not (self._xy) == None:
            #
            data['lat'] = self._xy[num,1]
            data['lon'] = self._xy[num,0]
            #
        else:
            #
            data['lat'] = 0.
            data['lon'] = 0.
            #
        #endif
        #
        data['spotterId'] = spotterID
        data['name'] = 'fake'

        timestamp    = time.strftime("%Y-%m-%dT%H:%M:00.000Z", time.gmtime() )
        data['timestampStart'] = timestamp
        data['timestampEnd'] = timestamp

        spotter =  spotres( data )
        spotter.fullMessage = False
        #
        if not (self.fakeData == None):
            #
            spotter.fullMessage = True
            spotter.spec = spectral.spectrum1d( { 'f':self.fakeData.ff(),
                                               'E':self.fakeData.E[num,:],
                                               'a1':self.fakeData.a1[num,:],
                                               'b1':self.fakeData.b1[num,:],
                                               'a2':self.fakeData.a2[num,:],
                                               'b2':self.fakeData.b2[num,:] } )
            #
        #endif
        #
        return(spotter)
        #
    #enddef
    #

    def updateFake( self , spotter):
        #
        # Update the time interval on a fake spotter
        #
        import time
        import calendar

        curtime = calendar.timegm(time.gmtime())

        updated = False
        if curtime - spotter.epochtime > self.fakeinterval:
            #
            spotter.epochtime   = curtime
            spotter.epochtimes  = [curtime,curtime]
            spotter.updatetime  = [ time.strftime("%Y-%m-%dT:%H:%M00.000Z", time.gmtime(x) ) for x in spotter.epochtimes ]
            updated = True            
            #
        #endif
        #
        return( updated )
        #
    #end updateFake
        
    #
    def update( self  ):
        #
        # This function tries to update all the spotters in the spotterlist with 
        # the latest data from the backend. If no data is available for a given
        # spotter, or the update fails, the spotter is removed from the available
        # list and we check if the spotter has spectral data or not, and keep
        # a seperate list of spotters with spectral data.
        #
        # Order of calculations:
        #
        # 1) FOR EACH SPOTTER DO
        #    1a) updata data in spotters from backend
        #           if succes, continue
        #           if fails, remove spotter from availability lists, go to next spotter
        #    1b) update the spectral availability list
        #    1c) if the data is more recent, add to the updated list
        #    1d) set spotter data to retrieve data
        # 2) update the mean time of the spotterlist
        # 3) update availability list; check if the latest data was recent or not
        #
        import time
        import calendar
        
        updated = []
        # 1) For each spotter do...
        for index, spotter in enumerate(self.spotters):
            #
            #
            #1a) update data
            try:
                #
                # Try and get data from the backend for the given spotter
                # 
                latest = self.getData( spotter.id,endpoint=self.endpoint,num=index )
                #
                if not (index in self.available):
                    #
                    # Back alive! Add this spotter to the available spotters
                    #
                    self.available.append( index )
                    self.available.sort()
                    #
                #
                #
            except:
                #
                # We failed to retrieve data..
                #
                print('Warning: cannot update spotter data')
                spotter.makeNaN()
                latest = self.spotters[ index ]
                #
                if index in self.available:
                    #
                    # MIA, presumed dead :(, Remove this spotter from the available spotters
                    #
                    self.available.remove(index)
                    #

                if index in self.fullSpecIndices:
                    #
                    # MIA, presumed dead :(, Remove this spotter from the available spotters with spectral data
                    #
                    self.fullSpecIndices.remove(index)
                    #                    
                continue
                #
            #
            # 
            # 1b) update spectral list
            if index in self.fullSpecIndices:
                #
                if not latest.fullMessage:                    
                    #
                    self.fullSpecIndices.remove( index )
                    self.numFullSpec += -1                        
                    #
                #endif
                #
            else:
                #
                if latest.fullMessage:
                    #
                    self.fullSpecIndices.append( index )
                    self.fullSpecIndices.sort()
                    self.numFullSpec += 1
                    #
                #endif
                #
            #endif
            #
            #
            # 1c) add to updated list if spotter data is new
            if latest.epochtime > spotter.epochtime:
                #
                updated.append( index )
                #
            #end if
            #
            # 1d) set spotter data to latest retrieved data
            self.spotters[ index ] = latest            
            #
        #end for
        #
            
        #
        # Update the mean and latest update time of
        # all spotters in the list
        #
        self.updateTime()
        #
        #loop over all spotter - and remove those for which
        #the last update was longer than 1.5 hours ago
        #from the last update of the spotter list
        #
        for index, spotter in enumerate(self.spotters):
            #
            if self.lastUpdate - spotter.epochtime > self.isaliveTime or \
                spotter.Tp > 34:
                #
                self.spotters[index].makeNaN()
                #
                if (index in self.available):
                    #
                    self.available.remove(index)
                    #
                #endif
                #
                if (index in self.fullSpecIndices):
                    #
                    self.fullSpecIndices.remove(index)
                    #
                #endif
                #
            #endif
            #
        #endfor
        #
        self.updateTime()        
        #
        # Return a list of all updated spotters            
        return( updated )
        #
    #end update
    #
    
    def getData(self,spotterID ='SPOT-ST0009',
                    endpoint='https://wavefleet.herokuapp.com/api/latestData',num=0):
        #
        import urllib.request
        import json
        import ndbc
        import numpy as np
        import h5py
        import spectral

        #If this is a fake, or virtual spotter, do update from fake
        if self.fake:
            spot = self.getFakeData( spotterID , num )
            return(spot)

        if spotterID[0:4].lower() == 'ndbc':
            #
            #
            ndbcID = spotterID[4:]
            spec,lat,lon,epochtime = ndbc.getLatestSpec(stationID=ndbcID)
            data = {'spec':spec,'lat':lat,'lon':lon,'epochtime':epochtime,'name':spotterID,'spotterId':spotterID}
            spot = spotres( data )
            #
        elif self.matlabSpotter:
            #
            # This is a matlab spotter, read from dir
            #
            fileName = self.matlabSpotterDir + '/' + spotterID + '.mat'

            dataFile = h5py.File( fileName , 'r' )

            time = np.array(dataFile.get('time'))[0,:]

            #Convert from matlab to Unix epoch, 719529 is the number of days
            #at 1/1/1970 according to matlab datenum( 1970 , 1 , 1)
            time = (time - 719529) * 3600. * 24

            #Find minimum index in arrays
            minimumIndex = np.argmin( np.abs( self.targetTime - time ) )
            epochtime = time[minimumIndex]
            #
            E = np.array(dataFile.get('E'))[:,minimumIndex]
            f = np.array(dataFile.get('f'))
            a1 = np.array(dataFile.get('a1'))[:,minimumIndex]
            b1 = np.array(dataFile.get('b1'))[:,minimumIndex]
            a2 = np.array(dataFile.get('a2'))[:,minimumIndex]
            b2 = np.array(dataFile.get('b2'))[:,minimumIndex]
            lat = dataFile.get('latitude')[0][0]
            lon = dataFile.get('longitude')[0][0]
            spec = spectral.spectrum1d({'f':f,'E':E,'a1':a1,'b1':b1,'a2':a2,'b2':b2} )
            data = {'spec':spec,'lat':lat,'lon':lon,'epochtime':epochtime,'name':spotterID,'spotterId':spotterID}
            spot = spotres( data )            
            #
        else:
            #
            req = urllib.request.Request(endpoint+ '?spotterId=' + spotterID)
            if self.token is not None:
                #
                req.add_header('token', self.token)
                #
            #endif
            #
            #...decode json...
            response = urllib.request.urlopen( req )            
            data = json.loads( response.read().decode('utf-8') )
            #
            #... and return a spot datatype
            print(data)
            spot = spotres( data['data'] )
            #
        #endif
        return( spot )
        #
    #

    def getIds( self ):
        #
        ids = []
        #
        for spotter in self.spotters:
            #
            ids.append( getattr( spotter , 'id' ) )
        #
        return(ids)
    
    #
    def appendcsv( self,directory ):
        #
        import utils
        import numpy as np
        #
        varnames = ['Hsig','Tm01','Dir','lat','lon','dirSprd','peakSprd','peakDir','Tp']
        #
        for varname in varnames:
            #
            filename = directory + varname + '.obs'
            #
            #
            ids = []
            data = []
            #
            for spotter in self.spotters:
                #
                ids.append( getattr( spotter , 'id' ) )

                tmp = getattr( spotter , varname )
                data.append( tmp )
                #
            #endfor
            #
            utils.writePoints( self.meanTime , filename , data , ids , '{}' )           
            #
            
            #
        #endfor
        #
        # Log the availability of spotters in the analysis
        #

            
        available = np.zeros( self.numSpot )
        fullSpecAvailable = np.zeros( self.numSpot )
        available[ self.available ] = 1.
        fullSpecAvailable[ self.fullSpecIndices ] = 1.         

        filename  = directory + 'available' + '.obs'
        utils.writePoints( self.meanTime , filename , available , ids , '{}' )
        filename  = directory + 'fullSpecAvailable' + '.obs'
        utils.writePoints( self.meanTime , filename , fullSpecAvailable , ids , '{}' )
        #
    #enddef
    #
#end class
