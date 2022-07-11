# !/usr/bin/env python
#%%
import datetime

from pyDARNmusic import load_fitacf
from pyDARNmusic import music

radar   = 'bks'
sDate   = datetime.datetime(2010,11,19,12,1)
eDate   = datetime.datetime(2010,11,19,12,2)
fit_sfx = "fitacf"
data_dir = f'/home/fran/code/SuperdarnW3usr/ForGitRepo/'
fitacf  = load_fitacf(radar,sDate,eDate,data_dir=data_dir)

dataObj = music.musicArray(fitacf,sTime=sDate,eTime=eDate,fovModel='IS')

import ipdb; ipdb.set_trace()

# %%
from pyDARNmusic.musicPlot import musicRTI
fig = musicRTI(dataObj,beam=12)

# %%
