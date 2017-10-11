pushUrls = ['https://wavefleet.herokuapp.com/api/prediction-data'] #,'http://b59826ed.ngrok.io/api/predictionData']
def pushToServer( epochtime , block, ids, pointdata,modelName,modelSource, contours=None, nretry=3,interval=10,token=None ):
    #
    import time
    #
    iretry = 0
    print( 'Pushing data to server:' )
    while( iretry < nretry ):
        #
        data   = createJSONblock( epochtime , block, ids, pointdata,modelName,modelSource , contours )
        status = post( data,modelSource, token )
        #
        if status.ok:
            #
            print( ' - success; attempt {} of {}'.format(iretry+1,nretry) )
            break
            #
        else:
            #
            print( ' - failed; attempt {} of {}'.format(iretry+1,nretry) )
            time.sleep( interval )
            iretry += 1
            #
        #
    return( status.ok )
    #
#enddef

def createJSONblock( epochtime , block, ids, pointdata,modelName,modelSource , contours=None ):
    #
    # create json block to push
    #
    import json
    import time
    import numpy as np

    timestamp    = time.strftime("%Y-%m-%dT%H:%M:00.000Z", time.gmtime( epochtime ) )

    pointSamples = dict()

    #block = removeNaN(blockdata)
    
    #
    for [index,spotID] in enumerate(ids):
        #
        # Convert directions
        #
        direction = 270 - float(pointdata['Dir'][index])
        #
        if direction < 0.:
            #
            direction += 360
            #
        #
        peakdirection = 270 - float(pointdata['peakDir'][index])
        #
        if peakdirection < 0.:
            #
            peakdirection += 360
            #
        #            
        pointSample = {
            'latitude' :float(pointdata['Yp'][index]),
            'longitude':float(pointdata['Xp'][index]),
            'meanPeriod':float(pointdata['Tm01'][index]),
            'peakPeriod':float(pointdata['peakPeriod'][index]),
            'height':float(pointdata['Hsig'][index]),
            'meanDirection':direction,
            'meanDirectionalSpread':float(pointdata['Dspr'][index]),
            'peakDirection':peakdirection,
            'peakDirectionalSpread':float(pointdata['peakSprd'][index])
            }
            
        pointSamples[spotID] = pointSample
        #
    #end loop creating pointsamples
    #

    size  = block['Yp'].size
    shape = block['Yp'].shape

    meshSamples = dict()
    
    directions = block['Dir'].reshape( ( size ) )


    notnan    = ~ np.isnan( directions)
    notnandir = directions[ notnan ]
    
    notnandir = 270 - notnandir    
    notnandir[ notnandir < 0 ]   = notnandir[ notnandir < 0 ] + 360
    notnandir[ notnandir > 360 ] = notnandir[ notnandir > 360 ] - 360
    directions[ notnan ] = notnandir

    hs = block['Hsig'].reshape( (size) )
    notnan    = ~ np.isnan( hs)
    notnanhs  = hs[ notnan ]
     
    notnanhs[ notnanhs < 0.0 ]   = 0.
    hs[ notnan ] = notnanhs
    
    meshSamples['latitude']   = block['Yp'].reshape( ( size ) ).tolist()
    meshSamples['longitude']  = block['Xp'].reshape( ( size ) ).tolist()
    meshSamples['height']     = hs.tolist() #block['Hsig'].reshape( ( size ) ).tolist()
    meshSamples['period']     = block['Tm01'].reshape( ( size ) ).tolist()
    meshSamples['direction']  = directions.tolist()

    #
    if ( 'Vel_x' in block ):
        #
        meshSamples['U']   = block['Vel_x'].reshape( ( size ) ).tolist()
        #
    #

    #
    if ( 'Vel_y' in block ):
        #
        meshSamples['V']   = block['Vel_y'].reshape( ( size ) ).tolist()
        #
    #
        
    for key in meshSamples:
        #
        for index,data in enumerate(meshSamples[key]):
            #
            if np.isnan(data):
                #
                meshSamples[key][index] = "NaN"
                #
            #
        #
    #

    meshSamples['meshStride'] = shape[0]    
    model = {'modelName':modelName, 'source':modelSource}

    if not contours==None:
        #
        data = {"data":
                    {"spotterIds":ids,          
                         "timestamp":timestamp,     
                         "pointSamples":pointSamples,
                         "meshSamples":meshSamples,
                         "modelParams":model,
                         "contourlines":[]
                    }
                }        
        #
    else:
        #
        data = {"data":
                    {"spotterIds":ids,          
                         "timestamp":timestamp,     
                         "pointSamples":pointSamples,
                         "meshSamples":meshSamples,
                         "modelParams":model,
                         "contourlines":[]
                    }
                }
        #
    #
    return( json.dumps(data) )
    #

def post( post_data,outname,token=None ):
    #
    import gzip
    import requests
    import zlib
    headers = dict()
    headers['Content-Encoding'] = 'gzip'
    headers['Content-Type'] = 'application/gzip'

    #
    if token is not None:
        #
        headers['token'] = token
        #
    else:
        #
        headers['token'] = 'a8db94573ddd8ab23a470493091866'
        #
    #
    
    data = gzip.compress( post_data.encode('utf8') )

    #write data to file
    fout = gzip.open( outname+'.gz','wb')
    fout.write( post_data.encode('utf8') )
    fout.close()
    success = False
    for pushUrl in pushUrls:
        #
        res = requests.post(pushUrl, data=data, headers=headers)
        #
        if res.ok:
            #
            success=True
            #
        #
    #

    return(res)

def removeNaN(data):
    #
    import numpy as np
    #
    for key in data:
        #
        if type( data[key] ) == np.ndarray:
            #
            msk = data[key] == -999.
            data[key][msk] = np.NaN
            #
        #
    #
    return(data)
#
