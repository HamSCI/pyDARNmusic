# !/usr/bin/env python
#%%
import datetime

from pydarn import time2datetime

from pyDARNmusic import load_fitacf
from pyDARNmusic import music

radar   = 'bks'
sDate   = datetime.datetime(2010,11,19,1)
eDate   = datetime.datetime(2010,11,19,23)
fit_sfx = "fitacf"
data_dir = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/'
fitacf  = load_fitacf(radar,sDate,eDate,data_dir=data_dir)

dataObj = music.musicArray(fitacf,sTime=sDate,eTime=eDate,fovModel='IS')

# import ipdb; ipdb.set_trace()

# %%
# from pyDARNmusic.musicPlot import musicRTI,musicFan
# # fig = musicRTI(dataObj,beam=7)
# plotTime = datetime.datetime(2010,11,19,14)
# fig = musicFan(dataObj,time=plotTime)

# %%
