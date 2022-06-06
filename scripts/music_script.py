#!/usr/bin/env python
import datetime

from pyDARNmusic import load_fitacf
from pyDARNmusic import music

radar   = 'bks'
sDate   = datetime.datetime(2010,11,19,12)
eDate   = datetime.datetime(2010,11,19,22)

fitacf  = load_fitacf(radar,sDate,eDate)

dataObj = music.musicArray(fitacf,fovModel='IS')

import ipdb; ipdb.set_trace()
