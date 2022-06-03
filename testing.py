#%%
import datetime

# from matplotlib import pyplot as plt 
# import numpy as np
from music import musicArray,musicDataObj
# from ptrFitacf import PtrFit
from radDataRead import radDataOpen


bmnum = None
# cp = None
fileType = 'fitacf'
# filtered = False
# src = None
rad='gbr'
syear = 2011
smonth = 12
sday = 9
shour = 8
sminute = 0

eyear = 2011
emonth = 12
eday = 9
ehour = 19
eminute = 41

# local_dirfmt = f'/media/fran/Expansion/PydarnW3usr/{syear}/{fileType}/data_vt_{rad}'
local_dirfmt = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/{syear}/{fileType}/data_vt_{rad}/'
# fileName = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/{syear}/{fileType}/data_vt_{rad}/20101119.1401.00.bks.fitacf.bz2'
sTime = datetime.datetime(syear,smonth,sday,shour)
eTime = datetime.datetime(eyear,emonth,eday,ehour)

myPtr = radDataOpen(sTime=sTime,radcode=rad, eTime=eTime, channel=None, bmnum=bmnum, cp=None,
                fileType=fileType, filtered=False, src=None, fileName=None,
                noCache=False, local_dirfmt=local_dirfmt, local_fnamefmt=None,
                local_dict=None, remote_dirfmt=None, remote_fnamefmt=None,
                remote_dict=None, remote_site=None, username=None,
                password=None, port=None, tmpdir=None, remove=False,
                try_file_types=True)

dataObj     = musicArray(myPtr,fovModel='GS')
# dataObj.get_data_sets()
# dataObj.DS000_originalFit.printMetadata()
# dataObj.DS000_originalFit.printHistory()
# dataObj.active.printHistory()
# from musicPlot import musicRTI
# fig = musicRTI(dataObj)
from musicPlot import musicFan
plotTime = datetime.datetime(2011,12,9,14)
fig = musicFan(dataObj,time=plotTime)
# from music import defineLimits, filterTimes
# defineLimits(dataObj,gateLimits=[30,45])
# fig = musicRTI(dataObj)

# We also to restrict the amount of time we process.  Before actually running the music algorithm, we
# will be filtering the data using a FIR filter that will eat up data at the beginning and end of the filter.
# We can calculate exactly how much time that will be if we know some of the filter characteristics and
# the time resolution of the data.  This can then be used to give us new start and end times.

# For now, let's say we are going to use a filter with 101 taps and a dataset with 120 s resolution.
# numtaps = 101
# timeres = 120

#Let's also say that we are interested in the MSTID feature between 1400 and 1600 UT.
# sTime_of_interest = datetime.datetime(2011,11,21,15)
# eTime_of_interest = datetime.datetime(2011,11,21,23)

#Now calculate the new start and end times...
# new_times = filterTimes(sTime_of_interest, eTime_of_interest, timeres, numtaps)
# defineLimits(dataObj,timeLimits=new_times)

# fig = musicRTI(dataObj)


# Now we apply the limits and replot once more.
# Note that many of the processing routines will automatically call applyLimits() before they
# run the processing algorithm.

# dataObj.active.applyLimits()
# fig = musicRTI(dataObj)
# %%
import datetime
from music import musicArray,musicDataObj
from radDataRead import radDataOpen


local_dirfmt = '/home/fran/code/SuperdarnW3usr/ForGitRepo/data_vt_bks/sd-data/2011/fitacf/bks/'
sTime = datetime.datetime(2011,11,21,15)
eTime = datetime.datetime(2011,11,21,23)

myPtr = radDataOpen(sTime=sTime,radcode='bks', eTime=eTime, channel=None, bmnum=7, cp=None,
                fileType='fitacf', filtered=False, src=None, fileName=None,
                noCache=False, local_dirfmt=local_dirfmt, local_fnamefmt=None,
                local_dict=None, remote_dirfmt=None, remote_fnamefmt=None,
                remote_dict=None, remote_site=None, username=None,
                password=None, port=None, tmpdir=None, remove=False,
                try_file_types=True)
dataObj     = musicArray(myPtr,fovModel='GS')
from musicPlot import musicFan
plotTime = datetime.datetime(2011,11,21,15)
fig = musicFan(dataObj,time=plotTime)
# %%
import datetime

from matplotlib import pyplot as plt 
import numpy as np
from music import musicArray,musicDataObj
from musicPlot import musicRTI
from radDataRead import radDataOpen
from music import defineLimits, filterTimes,beamInterpolation,timeInterpolation,determineRelativePosition
from musicPlot import musicFan,timeSeriesMultiPlot,plotRelativeRanges

#Choose the radar and time of the data that we want to look at.
#Then, connect to the data source using the standard DaViTPy
#pydarn.sdio.radDataOpen() routine.
# channel = 'a'
bmnum = 7
# cp = None
fileType = 'fitacf'
# filtered = False
# src = None
rad='bks'
syear = 2010
smonth = 11
sday = 19
shour = 14
sminute = 0

eyear = 2010
emonth = 11
eday = 19
ehour = 14
eminute = 1

# local_dirfmt = '/home/fran/code/SuperdarnW3usr/ForGitRepo/data_vt_bks/sd-data/2010/fitacf/bks/'
local_dirfmt = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/{syear}/{fileType}/data_vt_{rad}/'
fileName = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/{syear}/{fileType}/data_vt_{rad}/20101119.1401.00.bks.fitacf.bz2'
sTime = datetime.datetime(syear,smonth,sday,shour,sminute)
eTime = datetime.datetime(eyear,emonth,eday,ehour,eminute)

myPtr = radDataOpen(sTime=sTime,radcode=rad, eTime=eTime, channel=None, bmnum=bmnum, cp=None,
                fileType=fileType, filtered=False, src=None, fileName=fileName,
                noCache=False, local_dirfmt=None, local_fnamefmt=None,
                local_dict=None, remote_dirfmt=None, remote_fnamefmt=None,
                remote_dict=None, remote_site=None, username=None,
                password=None, port=None, tmpdir=None, remove=False,
                try_file_types=True)
            