import numpy as np
import datetime 
import copy
import logging

from ..utils.musicUtils import applyLimits

class musicDataObj(object):
    """This class is the basic container for holding MUSIC data.

    Parameters
    ---------- 
    time : list of datetime.datetime
        list of times corresponding to data
    data : numpy.array
        3-dimensional array of data
    fov : Optional[pydarn.radar.radFov.fov]
        Radar field-of-view object.
    comment : Optional[str]
        String to be appended to the history of this object
    parent : Optional[musicArray]
        reference to parent musicArray object
    **metadata
        keywords sent to matplot lib, etc.

    Attributes
    ----------
    time : numpy.array of datetime.datetime
        numpy array of times corresponding to data
    data : numpy.array
        3-dimensional array of data
    fov : Optional[pydarn.radar.radFov.fov]
        Radar field-of-view object.
    metadata : dict
        keywords sent to matplot lib, etc.
    history : dict 

    Methods
    ---------
    copy
    setActive
    nyquistFrequency
    samplePeriod
    applyLimits
    setMetadata
    printMetadata
    appendHistory
    printHistory

    Written by Nathaniel A. Frissell, Fall 2013
    Updated by: Francis Tholley, 2022
    """

    def __init__(self, time, data, fov=None, comment=None, parent=0, **metadata):
        self.parent = parent

        self.time     = np.array(time)
        self.data     = np.array(data)
        self.fov      = fov
        self.metadata = {}
        for key in metadata: self.metadata[key] = metadata[key]

        self.history = {datetime.datetime.now():comment}

    def copy(self,newsig,comment):
        """Copy a musicDataObj object.  This deep copies data and metadata, updates the serial
        number, and logs a comment in the history.  Methods such as plot are kept as a reference.

        Parameters
        ----------
        newsig : str
            Name for the new musicDataObj object.
        comment : str
            Comment describing the new musicDataObj object.

        Returns
        -------
        newsigobj : musicDataObj
            Copy of the original musicDataObj with new name and history entry.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """

        serial = self.metadata['serial'] + 1
        newsig = '_'.join(['DS%03d' % serial,newsig])

        setattr(self.parent,newsig,copy.copy(self))
        newsigobj = getattr(self.parent,newsig)

        newsigobj.time      = copy.deepcopy(self.time)
        newsigobj.data      = copy.deepcopy(self.data)
        newsigobj.fov       = copy.deepcopy(self.fov)
        newsigobj.metadata  = copy.deepcopy(self.metadata)
        newsigobj.history   = copy.deepcopy(self.history)

        newsigobj.metadata['dataSetName'] = newsig
        newsigobj.metadata['serial']      = serial
        newsigobj.history[datetime.datetime.now()] = '['+newsig+'] '+comment
        
        return newsigobj
  
    def setActive(self):
        """Sets this signal as the currently active signal.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        self.parent.active = self

    def nyquistFrequency(self,timeVec=None):
        """Calculate the Nyquist frequency of a vt sigStruct signal.

        Parameters
        ----------
        timeVec : Optional[list of datetime.datetime]
            List of datetime.datetime to use instead of self.time.

        Returns
        -------
        nq : float
            Nyquist frequency of the signal in Hz.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """

        dt  = self.samplePeriod(timeVec=timeVec)
        nyq = float(1. / (2*dt))
        return nyq

    def samplePeriod(self,timeVec=None):
        """Calculate the sample period of a vt sigStruct signal.

        Parameters
        ----------
        timeVec : Optional[list of datetime.datetime]
            List of datetime.datetime to use instead of self.time.

        Returns
        -------
        samplePeriod : float
            samplePeriod: sample period of signal in seconds.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        
        if timeVec  == None: timeVec = self.time

        diffs       = np.diff(timeVec)
        diffs_unq   = np.unique(diffs)
        self.diffs  = diffs_unq

        if len(diffs_unq) == 1:
            samplePeriod = diffs[0].total_seconds()
        else:
            diffs_sec   = np.array([x.total_seconds() for x in diffs])
            maxDt       = np.max(diffs_sec)
            avg         = np.mean(diffs_sec)

            md          = self.metadata
            warn        = 'WARNING'
            if 'title' in md: warn = ' '.join([warn,'FOR','"'+md['title']+'"'])
            logging.warning(warn + ':')
            logging.warning('   Date time vector is not regularly sampled!')
            logging.warning('   Maximum difference in sampling rates is ' + str(maxDt) + ' sec.')
            logging.warning('   Using average sampling period of ' + str(avg) + ' sec.')
            samplePeriod = avg
            

        return samplePeriod

    def applyLimits(self,rangeLimits=None,gateLimits=None,timeLimits=None,newDataSetName='limitsApplied',comment='Limits Applied'):
        """Removes data outside of the rangeLimits, gateLimits, and timeLimits boundaries.

        Parameters
        ----------
        rangeLimits : Optional[interable]
            Two-element array defining the maximum and minumum slant ranges to use. [km]
        gateLimits : Optional[iterable]
            Two-element array defining the maximum and minumum gates to use.
        timeLimits :  Optional[]

        newDataSetName : Optional[str]
            Name of the new musicDataObj to be created in the current musicArray object as a result of this processing.
        comment : Optional[str]
            String to be appended to the history of this object.

        Returns
        -------
        newMusicDataObj : musicDataObj
            New musicDataObj.  The musicDataObj is also stored in it's parent musicArray object.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        return applyLimits(self.parent,self.metadata['dataSetName'],rangeLimits=rangeLimits,gateLimits=gateLimits,timeLimits=timeLimits,newDataSetName=newDataSetName,comment=comment)

    def setMetadata(self,**metadata):
        """Adds information to the current musicDataObj's metadata dictionary.
        Metadata affects various plotting parameters and signal processing routinges.

        Parameters
        ----------
        **metadata :
            keywords sent to matplot lib, etc.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        self.metadata = dict(list(self.metadata.items()) + list(metadata.items()))

    def printMetadata(self):
        """Nicely print all of the metadata associated with the current musicDataObj object.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        keys = list(self.metadata.keys())
        keys.sort()
        for key in keys:
            print(key+':',self.metadata[key])

    def appendHistory(self,comment):
        """Add an entry to the processing history dictionary of the current musicDataObj object.

        Parameters
        ----------
        comment : string
            Infomation to add to history dictionary.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        self.history[datetime.datetime.now()] = '['+self.metadata['dataSetName']+'] '+comment

    def printHistory(self):
        """Nicely print all of the processing history associated with the current musicDataObj object.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        keys = list(self.history.keys())
        keys.sort()
        for key in keys:
            print(key,self.history[key])
