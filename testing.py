#%%
import datetime

# from matplotlib import pyplot as plt 
# import numpy as np
from music import musicArray,musicDataObj
# from ptrFitacf import PtrFit
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
dataObj.get_data_sets()
dataObj.DS000_originalFit.printMetadata()
dataObj.DS000_originalFit.printHistory()
dataObj.active.printHistory()
from musicPlot import musicRTI
fig = musicRTI(dataObj)
# scan1 = myPtr.readScan() 
# print(scan1)
# %%
