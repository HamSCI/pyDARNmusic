#!/usr/bin/env python3
import os
import bz2
import glob
import datetime
import tqdm

import numpy as np

import pydarnio
from pydarn import time2datetime

def load_fitacf(radar,sTime,eTime=None,data_dir='/sd-data',fit_sfx='fitacf'):
    """
    Load FITACF data from multiple FITACF files by specifying a date range.

    This routine assumes bz2 compression.
    """
    
    if eTime is None:
        eTime = sTime + datetime.timedelta(days=1)

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
        date_str        = date.strftime('%Y%m%d')
        year_str        = date.strftime('%Y')
        fpattern        = os.path.join(data_dir,year_str,fit_sfx,radar,'{!s}*{!s}*.{!s}.bz2'.format(date_str,radar,fit_sfx))
        fitacf_paths_0   += glob.glob(fpattern)


    # Sort the files by name.
    fitacf_paths_0.sort()

    # Get rid of any files we don't need.
    fitacf_paths = []
    time_deltas  = []
    for finx, fpath in enumerate(fitacf_paths_0):
        date_str    = os.path.basename(fpath)[:13]
        this_time   = datetime.datetime.strptime(date_str,'%Y%m%d.%H%M')

        # Eliminate files after the time we need
        if this_time <= eTime:
            fitacf_paths.append(fpath)
            time_deltas.append((sTime-this_time).total_seconds())

    # Eliminate files before the time that we need
    tds     = np.array(time_deltas)
    min_inx = np.argmin(tds[tds>=0])
    fitacf_paths = fitacf_paths[min_inx:]

    # Load and append each data file.
    print()
    fitacf = []
    for fitacf_path in tqdm.tqdm(fitacf_paths,desc='Loading {!s} Files'.format(fit_sfx),dynamic_ncols=True):
        tqdm.tqdm.write(fitacf_path)
        with bz2.open(fitacf_path) as fp:
            fitacf_stream = fp.read()

        reader  = pydarnio.SDarnRead(fitacf_stream, True)
        records = reader.read_fitacf()
        fitacf += records

    # Remove uneeded fitacf records.
    fitacf_new = []
    for record in fitacf:
        this_time = time2datetime(record)
        if this_time >= sTime and this_time < eTime:
            fitacf_new.append(record)
    return fitacf_new

if __name__ == '__main__':
    output_dir  = 'plots'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sTime       = datetime.datetime(2010,11,1)
    eTime       = datetime.datetime(2010,11,2)
    radar       = 'bks'
    fitacf      = load_fitacf(sTime,eTime,radar)
