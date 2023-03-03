#Import the modules we need.
import datetime

from matplotlib import pyplot as plt 
import numpy as np
from pyDARNmusic import musicRTP,musicFan,timeSeriesMultiPlot,plotRelativeRanges,spectrumMultiPlot,plotFullSpectrum,plotDlm, plotKarr,plotKarrDetected
from pyDARNmusic import load_fitacf
from pyDARNmusic import music
from pyDARNmusic import (getDataSet 
                              ,stringify_signal,stringify_signal_list     
                               ,beamInterpolation         
                               ,defineLimits              
                               ,checkDataQuality          
                               ,applyLimits               
                               ,determineRelativePosition 
                               ,timeInterpolation         
                               ,filterTimes               
                               ,detrend                   
                               ,nan_to_num                
                               ,windowData                
                               ,calculateFFT              
                               ,calculateDlm              
                               ,calculateKarr             
                               ,simulator                 
                               ,scale_karr                
                               ,detectSignals             
                               ,add_signal                
                               ,del_signal)
# from music import defineLimits, filterTimes,beamInterpolation,timeInterpolation,determineRelativePosition


import pydarn
import pydarnio


radar   = 'bpk'
sDate   = datetime.datetime(2017,1,15,1)
eDate   = datetime.datetime(2017,1,15,8)
fit_sfx = "fitacf"
data_dir = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/'
fitacf  = load_fitacf(radar,sDate,eDate,data_dir=data_dir)
dataObj = music.musicArray(fitacf,sTime=sDate,eTime=eDate,fovModel='GS')
defineLimits(dataObj,gateLimits=[27,42])

numtaps = 101
timeres = 120
# numtaps = 8
# timeres = 7

#Let's also say that we are interested in the MSTID feature between 1400 and 1600 UT.
sTime_of_interest = datetime.datetime(2017,1,15,5,20)
eTime_of_interest = datetime.datetime(2017,1,15,8,30)

#Now calculate the new start and end times...
new_times = filterTimes(sTime_of_interest, eTime_of_interest, timeres, numtaps)
defineLimits(dataObj,timeLimits=new_times)
dataObj.active.applyLimits()
beamInterpolation(dataObj)
timeInterpolation(dataObj,timeRes=timeres)
determineRelativePosition(dataObj)
filt = music.filter(dataObj, numtaps=numtaps, cutoff_low=0.0003, cutoff_high=0.0012)
dataObj.active.applyLimits()
calculateFFT(dataObj)
# import ipdb;ipdb.set_trace()
calculateDlm(dataObj)
import time
time1 = time.time()
calculateKarr(dataObj)
time2 = time.time() - time1
print("py running time in seconds:", time2)
