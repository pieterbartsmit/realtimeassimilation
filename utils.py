def loadPoints( filename , epochtime=None ):
    #
    import os.path
    import calendar
    import time
    import pandas
    import numpy

    if ( not os.path.isfile(filename) ):
        #
        return( [] )

    data = pandas.read_csv( filename,engine='python', header=0 )
    data.columns = data.columns.str.strip()

    if epochtime==None:
        #
        data = numpy.array( data.iloc[-1].tolist() )[2:]
        #
    return( data)
    
def writePoints(  epochtime , filename , data , dataNames , formatting ):
    #
    # Apped data to a text file
    #
    import os.path
    import calendar
    import time
    #

    #
    # First, check if the file already exists...
    #
    append = False
    if ( not os.path.isfile(filename) ):
        #
        # ...if not, we create it ...
        #
        fid = open( filename , 'w' )
        #
        # ...and write the header ...
        #
        fid.write( '  date/time , epochtime ')
        for name in dataNames:
            #
            fid.write( ', ' + str(name) )
            #
        #endif
        #
        fid.write('\n')
        #
        # and close
        #
        fid.close()
        append = True
        #
    #endif
    #

    #
    # Read the existing contents
    #
    with open( filename , 'r' ) as fid:
        # read contents
        lines = fid.readlines()
        lines  =[line.replace('\n','') for line in lines]

    #
    # Check the last line of the existing contents
    #
    lastline = lines[-1]
    datetime     = lastline.split(',')[0].strip()
    timestr      = time.strftime( "%Y%m%d%H%M" , time.gmtime( int(epochtime) ) )

    #
    # Did we already write data for this time slot? (can happen if we e.g. recover from crash)
    #
    if datetime == timestr:
        #
        # if so we are going to overwrite
        #
        append = False
        #
    else:
        #
        # else we add
        #
        append = True

    #
    # Write the datastring
    # 
    datastr = timestr +', ' + str(int(epochtime))
    for dat in data:
        #
        datastr = datastr + ', ' + formatting.format( dat )
        #
    #endfor

    #
    # add or append to the line list
    #
    if append:
        lines.append(datastr)
    else:
        lines[-1]=datastr    

    #
    # and write to disk again
    #
    with open( filename , 'w' ) as fid:
        #
        for line in lines:
            fid.write( line + '\n' )
        #
    #end with
    #
#end def

def get_contours( x,y,z,levellist):
    #
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax1 = fig.add_subplot(1,1,1)

    cn  = ax1.contour(x,y,z,levellist)
    matplotlib.pyplot.close(fig)
    contours = []
    #
    # for each contour line
    for [index,cc] in enumerate(cn.collections):
        #
        #
        for [index2,pp] in enumerate(cc.get_paths()):
            #
            xy =pp.vertices
            if len(xy) > 10:
                #
                path = dict()
                path['lon'] = xy[:,0].tolist()
                path['lat'] = xy[:,1].tolist()
                path['depth'] = levellist[index]
                contours.append( path )
                #
            #end if
            #
        #end for
        #
        #
    #end
    #
    return contours
    #
#end def
#

