import numpy as np
import copy
import logging

class filter(object):
    """Filter a VT sig/sigStruct object and define a FIR filter object.
    If only cutoff_low is defined, this is a high pass filter.
    If only cutoff_high is defined, this is a low pass filter.
    If both cutoff_low and cutoff_high is defined, this is a band pass filter.

    Uses scipy.signal.firwin()
    High pass and band pass filters inspired by Matti Pastell's page:
      http://mpastell.com/2010/01/18/fir-with-scipy/

    Metadata keys:
      'filter_cutoff_low'   --> cutoff_low
      'filter_cutoff_high'  --> cutoff_high
      'filter_numtaps'      --> cutoff_numtaps

    Parameters
    ----------
    dataObj : musicArray
        musicArray object
    dataSet : Optional[str]
        which dataSet in the musicArray object to process
    numtaps : Optional[int]
        Length of the filter (number of coefficients, i.e. the filter
        order + 1).  `numtaps` must be even if a passband includes the
        Nyquist frequency.

        If dataObj.dataSet.metadata['filter_numptaps'] is set and this keyword is None,
        the metadata value will be used.

    cutoff_low : Optional[float, 1D array_like or None]
        High pass cutoff frequency of filter (expressed in the same units as `nyq`)
        OR an array of cutoff frequencies (that is, band edges). In the
        latter case, the frequencies in `cutoff` should be positive and
        monotonically increasing between 0 and `nyq`.  The values 0 and
        `nyq` must not be included in `cutoff`. If None, a low-pass filter will not
        be applied.

        If dataObj.dataSet.metadata['filter_cutoff_low'] is set and this keyword is None,
        the metadata value will be used.

    cutoff_high : Optional[float, 1D array_like, or None]
        Like cutoff_low, but this is the low pass cutoff frequency of the filter.

        If dataObj.dataSet.metadata['filter_cutoff_high'] is set and this keyword is None,
        the metadata value will be used.

    width : Optional[float]
        If `width` is not None, then assume it is the approximate width
        of the transition region (expressed in the same units as `nyq`)
        for use in Kaiser FIR filter design.  In this case, the `window`
        argument is ignored.

    window : Optional[string or tuple of string and parameter values]
        Desired window to use. See `scipy.signal.get_window` for a list
        of windows and required parameters.

    pass_zero : Optional[bool]
        If True, the gain at the frequency 0 (i.e. the "DC gain") is 1.
        Otherwise the DC gain is 0.

    scale : Optional[bool]
        Set to True to scale the coefficients so that the frequency
        response is exactly unity at a certain frequency.
        That frequency is either:
            0 (DC) if the first passband starts at 0 (i.e. pass_zero is True);
            nyq` (the Nyquist rate) if the first passband ends at
            `nyq` (i.e the filter is a single band highpass filter);
            center of first passband otherwise.

    Attributes
    ----------
    comment : str

    cutoff_low : float, 1D array_like or None
        High pass cutoff frequency of filter (expressed in the same units as `nyq`)
        OR an array of cutoff frequencies (that is, band edges).
    cutoff_high : float, 1D array_like, or None
        Like cutoff_low, but this is the low pass cutoff frequency of the filter.
    nyq : float
        the Nyquist rate
    ir : 

    Methods
    -------
    plotTransferFunction
    plotImpulseResponse
    filter

    Written by Nathaniel A. Frissell, Fall 2013
    Updated by: Francis Tholley, 2022
    """
    def __init__(self, dataObj, dataSet='active', numtaps=None, cutoff_low=None, cutoff_high=None, width=None, window='blackman', pass_zero=True, scale=True,newDataSetName='filtered'):
        from scipy import signal
        sigObj = getattr(dataObj,dataSet)
        nyq = sigObj.nyquistFrequency()

        #Get metadata for cutoffs and numtaps.
        md = sigObj.metadata
        if cutoff_high == None:
            if 'filter_cutoff_high' in md:
                cutoff_high = md['filter_cutoff_high']

        if cutoff_low == None:
            if 'filter_cutoff_low' in md:
                cutoff_low = md['filter_cutoff_low']

        if numtaps == None:
            if 'filter_numtaps' in md:
                numtaps = md['filter_numtaps']
            else:
                logging.warning('You must provide numtaps.')
                return


        if   cutoff_high != None:    #Low pass
            lp = signal.firwin(numtaps=numtaps, cutoff=cutoff_high, width=width, window=window, pass_zero=pass_zero, scale=scale, nyq=nyq)
            d = lp

        if   cutoff_low != None:    #High pass
            hp = -signal.firwin(numtaps=numtaps, cutoff=cutoff_low, width=width, window=window, pass_zero=pass_zero, scale=scale, nyq=nyq)
            hp[numtaps//2] = hp[numtaps//2] + 1
            d = hp

        if cutoff_high != None and cutoff_low != None:
            d = -(lp+hp)
            d[numtaps//2] = d[numtaps//2] + 1
            d = -1.*d #Needed to correct 180 deg phase shift.

        if cutoff_high == None and cutoff_low == None:
            logging.warning("You must define cutoff frequencies!")
            return
    
        self.comment = ' '.join(['Filter:',window+',','Nyquist:',str(nyq),'Hz,','Cuttoff:','['+str(cutoff_low)+', '+str(cutoff_high)+']','Hz,','Numtaps:',str(numtaps)])
        self.cutoff_low   = cutoff_low
        self.cutoff_high  = cutoff_high
        self.nyq = nyq
        self.ir = d

        self.filter(dataObj,dataSet=dataSet,newDataSetName=newDataSetName)


    def __str__(self):
        return self.comment

    def plotTransferFunction(self,xmin=0,xmax=None,ymin_mag=-150,ymax_mag=5,ymin_phase=None,ymax_phase=None,worN=None,fig=None):
        from scipy import signal
        """Plot the frequency and phase response of the filter object.

        Parameters
        ----------
        xmin : Optional[float]
            Minimum value for x-axis.
        xmax : Optional[float]
            Maximum value for x-axis.
        ymin_mag : Optional[float]
            Minimum value for y-axis for the frequency response plot.
        ymax_mag : Optional[float]
            Maximum value for y-axis for the frequency response plot.
        ymin_phase : Optional[float]
            Minimum value for y-axis for the phase response plot.
        ymax_phase : Optional[float]
            Maximum value for y-axis for the phase response plot.
        worN : Optional[int]
            passed to scipy.signal.freqz()
            If None, then compute at 512 frequencies around the unit circle.
            If the len(filter) > 512, then compute at len(filter) frequencies around the unit circle.
            If a single integer, the compute at that many frequencies.
            Otherwise, compute the response at frequencies given in worN
        fig : Optional[matplotlib.Figure]
            Figure object on which to plot.  If None, a figure will be created.

        Returns
        -------
        fig : matplotlib.Figure
            Figure object containing the plot.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """

        if fig == None:
            from matplotlib import pyplot as plt
            fig   = plt.figure(figsize=(20,10))

        if worN == None:
            if len(self.ir) > 512: worN = len(self.ir)
            else: worN = None
        else: pass

        w,h = signal.freqz(self.ir,1,worN=worN)
        h_dB = 20 * np.log10(abs(h))
        axis = fig.add_subplot(211)

        #Compute frequency vector.
        w = w/max(w) * self.nyq
        axis.plot(w,h_dB,'.-')
        #mp.axvline(x=self.fMax,color='r',ls='--',lw=2)

        if xmin is not None: axis.set_xlim(xmin=xmin)
        if xmax is not None: axis.set_xlim(xmax=xmax)
        if ymin_mag is not None: axis.set_ylim(ymin=ymin_mag)
        if ymax_mag is not None: axis.set_ylim(ymax=ymax_mag)

        axis.set_xlabel(r'Frequency (Hz)')
        axis.set_ylabel('Magnitude (db)')

        axis.set_title(r'Frequency response')

        axis = fig.add_subplot(212)
        h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
        axis.plot(w,h_Phase,'.-')

        if xmin is not None: axis.set_xlim(xmin=xmin)
        if xmax is not None: axis.set_xlim(xmax=xmax)
        if ymin_phase is not None: axis.set_ylim(ymin=ymin_phase)
        if ymax_phase is not None: axis.set_ylim(ymax=ymax_phase)

        axis.set_ylabel('Phase (radians)')
        axis.set_xlabel(r'Frequency (Hz)')
        axis.set_title(r'Phase response')
        fig.suptitle(self.comment)
        fig.subplots_adjust(hspace=0.5)

        return fig

    def plotImpulseResponse(self,xmin=None,xmax=None,ymin_imp=None,ymax_imp=None,ymin_step=None,ymax_step=None,fig=None):
        from scipy import signal
        """Plot the frequency and phase response of the filter object.

        Parameters
        ----------
        xmin : Optional[float]
            Minimum value for x-axis.
        xmax : Optional[float]
            Maximum value for x-axis.
        ymin_imp : Optional[float]
            Minimum value for y-axis for the impulse response plot.
        ymax_imp : Optional[float]
            Maximum value for y-axis for the impulse response plot.
        ymin_step : Optional[float]
            Minimum value for y-axis for the step response plot.
        ymax_step : Optional[float]
            Maximum value for y-axis for the step response plot.
        fig : Optional[matplotlib.Figure]
            Figure object on which to plot.  If None, a figure will be created.

        Returns
        -------
        fig : matplotlib.Figure
            Figure object containing the plot.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """

        if fig == None:
            from matplotlib import pyplot as plt
            fig   = plt.figure(figsize=(20,10))

        l = len(self.ir)
        impulse = np.repeat(0.,l); impulse[0] =1.
        x = np.arange(0,l)
        response = signal.lfilter(self.ir,1,impulse)
        axis = fig.add_subplot(211)
        axis.stem(x, response)
        axis.set_ylabel('Amplitude')
        axis.set_xlabel(r'n (samples)')
        axis.set_title(r'Impulse response')

        axis = fig.add_subplot(212)
        step = np.cumsum(response)
        axis.stem(x, step)
        axis.set_ylabel('Amplitude')
        axis.set_xlabel(r'n (samples)')
        axis.set_title(r'Step response')
        fig.suptitle(self.comment)
        fig.subplots_adjust(hspace=0.5)

        return fig

    def filter(self,dataObj,dataSet='active',newDataSetName='filtered'):
        """Apply the filter to a vtsig object.

        Parameters
        ----------
        dataObj : musicArray
            musicArray object
        dataSet : Optional[str]
            which dataSet in the musicArray object to process
        newDataSetName : Optional[str]
            Name of the new musicDataObj to be created in the current musicArray object as a result of this processing.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """
        from scipy import signal

        sigobj = getattr(dataObj,dataSet)
        vtsig  = sigobj.parent

        nrTimes,nrBeams,nrGates = np.shape(sigobj.data)

        #Filter causes a delay in the signal and also doesn't get the tail end of the signal...  Shift signal around, provide info about where the signal is valid.
        shift = np.int32(-np.floor(len(self.ir)/2.))

        start_line    = np.zeros(nrTimes)
        start_line[0] = 1
        start_line    = np.roll(start_line,shift)

        tinx0 = abs(shift)
        tinx1 = np.where(start_line == 1)[0][0]

        val_tm0 = sigobj.time[tinx0]
        val_tm1 = sigobj.time[tinx1]

        filteredData = np.zeros_like(sigobj.data)

        #Apply filter
        for bm in range(nrBeams):
            for rg in range(nrGates):
                tmp = signal.lfilter(self.ir,[1.0],sigobj.data[:,bm,rg])
                tmp = np.roll(tmp,shift)
                filteredData[:,bm,rg] = tmp[:]

        #Create new signal object.
        newsigobj = sigobj.copy(newDataSetName,self.comment)
        #Put in the filtered data.
        newsigobj.data = copy.copy(filteredData)
        newsigobj.time = copy.copy(sigobj.time)

        #Clear out ymin and ymax from metadata; make sure meta data block exists.
        #If not, create it.

        if hasattr(newsigobj,'metadata'):
            delMeta = ['ymin','ymax','ylim']
            for key in delMeta:
                if key in newsigobj.metadata:
                    del newsigobj.metadata[key]
        else:
            newsigobj.metadata = {}

        newsigobj.metadata['timeLimits'] = (val_tm0,val_tm1)

        key = 'title'
        if key in newsigobj.metadata:
            newsigobj.metadata[key] = ' '.join(['Filtered',newsigobj.metadata[key]])
        else:
            newsigobj.metadata[key] = 'Filtered'

        newsigobj.metadata['fir_filter'] = (self.cutoff_low,self.cutoff_high)
        newsigobj.setActive()
