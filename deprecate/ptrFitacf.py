#!/usr/bin/env python3
import os
import bz2
import glob
import datetime
import tqdm
import numpy as np
import pydarnio
import pydarn
import copy

import matplotlib as mpl
mpl.use('Agg')
# from matplotlib import pyplot as plt
# from pydarn import Coords
# import pandas as pd
class PtrFit():
    def __init__(self, radar = None, sTime=None,eTime=None):
        self.radar = radar
        self.sTime = sTime
        self.eTime = eTime
        self.fitacf = self.getFit()
        self.__hdw_data    = pydarn.read_hdw_file(self.radar,self.sTime)
        self.n_beams     = self.__hdw_data.beams

    def getFit(self):
        thisdate = self.sTime
        edate = self.eTime
        fitacf = None
        while thisdate < edate:
            sTime       = thisdate
            eTime       = sTime + datetime.timedelta(days = 1)

            fitacf      = self.load_fitacf(sTime,eTime,self.radar)
            thisdate += datetime.timedelta(days = 1) 
        return fitacf

    def load_fitacf(self,sTime,eTime,radar,data_dir='/home/fran/code/SuperdarnW3usr/ForGitRepo/data_vt_mcm/sd-data/2011/fitacf/mcm',fit_sfx='fitacf'):
        """
        Load FITACF data from multiple FITACF files by specifying a date range.

        This routine assumes bz2 compression.
        """
        sDate   = datetime.datetime(sTime.year,sTime.month,sTime.day)
        eDate   = datetime.datetime(eTime.year,eTime.month,eTime.day)

        # Create a list of days we need fitacf files from.
        dates   = [sDate]
        while dates[-1] < eDate:
            next_date   = dates[-1] + datetime.timedelta(days=1)
            dates.append(next_date)

        # Find the data files that fall in that date range.
        fitacf_paths_0    = []
        for date in dates:
            #radarType = 'a'
            date_str        = date.strftime('%Y%m%d')
            fpattern        = os.path.join(data_dir,'{!s}*{!s}*.{!s}.bz2'.format(date_str,radar,fit_sfx))
            fitacf_paths_0   += glob.glob(fpattern)

        # Sort the files by name.
        fitacf_paths_0.sort()

        # Get rid of any files we don't need.
        fitacf_paths = []
        for fpath in fitacf_paths_0:
            date_str    = os.path.basename(fpath)[:13]
            this_time   = datetime.datetime.strptime(date_str,'%Y%m%d.%H%M')

            if this_time <= eTime:
                fitacf_paths.append(fpath)

        # Load and append each data file.

        fitacf = []
        for fitacf_path in tqdm.tqdm(fitacf_paths,desc='Loading {!s} Files'.format(fit_sfx),dynamic_ncols=True):
            tqdm.tqdm.write(fitacf_path)
            with bz2.open(fitacf_path) as fp:
                fitacf_stream = fp.read()

            reader  = pydarnio.SDarnRead(fitacf_stream, True)
            records = reader.read_fitacf()
            fitacf += records
        return fitacf