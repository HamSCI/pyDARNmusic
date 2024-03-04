import datetime
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np

from pydarn import (Re, time2datetime, Coords, SuperDARNRadars,RangeEstimation)
from ..utils.geoPack import (greatCircleDist,greatCircleMove,greatCircleAzm)
from .music_data_object import musicDataObj

class musicArray(object):
    """This class is the basic container for holding MUSIC data.

    Parameters
    ----------
    myPtr : pydarn.sdio.radDataTypes.radDataPtr
        contains the pipeline to the data we are after
    sTime : Optional[datetime.datetime]
        start time UT (if None myPtr.sTime is used)
    eTime : Optional[datetime.datetime]
        end time UT (if None myPtr.eTime is used)
    param : Optional[str]
        Radar FIT parameter to load and process.  Any appropriate attribute of the
        FIT data structure is allowed.
    gscat : Optional[int]
        Ground scatter flag.
            0: all backscatter data 
            1: ground backscatter only
            2: ionospheric backscatter only
            3: all backscatter data with a ground backscatter flag.
    fovElevation : Optional[float]
        Passed directly to pydarn.radar.radFov.fov()
    fovModel : Optional[str]
        Scatter mapping model.
        GS : Ground Scatter Mapping Model.  See Bristow et al. [1994] (default)
        IS : Standard SuperDARN scatter mapping model.
        S  : Standard projection model
        E1 : for Chisham E-region 1/2-hop ionospheric projection model
        F1 : for Chisham F-region 1/2-hop ionospheric projection model
        F3 : for Chisham F-region 1 1/2-hop ionospheric projection model
        C  : Chisham projection model
        None : if you trust your elevation or altitude values
    fovCoords : Optional[str]
        Map coordinate system. WARNING: 'geo' is curently only tested coordinate system.
    full_array : Optional[bool]
        If True, make the data array the full beam, gate dimensions listed in the hdw.dat file.
        If False, truncate the array to the maximum dimensions that there is actually data.
        False will save space without throwing out any data, but sometimes it is easier to work
        with the full-size array.

    Attributes
    ----------
    messages : list

    prm : 

    Methods
    -------
    get_data_sets

    Example
    -------
        #Set basic event parameters.
        rad         ='wal'
        sTime       = datetime.datetime(2011,5,9,8,0)
        eTime       = datetime.datetime(2011,5,9,19,0)
        #Connect to a SuperDARN data source.
        myPtr       = pydarn.sdio.radDataOpen(sTime,rad,eTime=eTime)
        #Create the musicArray Object.
        dataObj     = music.musicArray(myPtr,fovModel='GS')

    References
    ----------
    Bristow, W. A., R. A. Greenwald, and J. C. Samson (1994), Identification of high-latitude acoustic gravity wave sources
        using the Goose Bay HF Radar, J. Geophys. Res., 99(A1), 319-331, doi:10.1029/93JA01470.

    Written by Nathaniel A. Frissell, Fall 2013
    Updated by: Francis Tholley, 2022
    """

    def __init__(self,fitacf,sTime=None,eTime=None,param='p_l',gscat=1,
            fovElevation=None,fovModel='GS',fovCoords='geo',full_array=True, scan_index= 1, channel = 'all',file_type='fitacf'):    
        self.messages   = []

        no_data_message = 'No data for this time period.'
        # If no data, report and return.
        if fitacf is None:
            self.messages.append(no_data_message)
            return

        if len(fitacf) == 0:
            self.messages.append(no_data_message)
            return

        fitacf_times = []
        for record in fitacf:
            this_time = datetime.datetime(record['time.yr'],record['time.mo'],record['time.dy'],
                            record['time.hr'],record['time.mt'],record['time.sc'],record['time.us'])
            fitacf_times.append(this_time)

        if sTime == None: sTime = min(fitacf_times)
        if eTime == None: eTime = max(fitacf_times)

        scanTimeList    = []
        dataList        = []
        cpidList        = []
        #Subscripts of columns in the dataList/dataArray
        scanInx = 0
        dateInx = 1
        beamInx = 2
        gateInx = 3
        dataInx = 4

        beamTime    = sTime
        fov         = None

        # Get scan numbers for each record

        # Create a place to store the prm data.
        prm                 = {}
        prm['time']         = []
        prm['mplgs']        = []
        prm['nave']         = []
        prm['noisesearch']  = []
        prm['scan']         = []
        prm['smsep']        = []
        prm['mplgexs']      = []
        prm['xcf']          = []
        prm['noisesky']     = []
        prm['rsep']         = []
        prm['mppul']        = []
        prm['inttsc']       = []
        prm['frang']        = []
        prm['bmazm']        = []
        prm['lagfr']        = []
        prm['ifmode']       = []
        prm['noisemean']    = []
        prm['tfreq']        = []
        prm['inttus']       = []
        prm['rxrise']       = []
        prm['mpinc']        = []
        prm['nrang']        = []

        stid = None
        radCode = None
        cp = None
        
        # Get scan numbers for each record
        # beam_scan = build_scan(fitacf)
        # scan_ids = np.unique(beam_scan) # get the total number of scans by getting only the unique num of scan
        
        # loop through all of the scan numbers for the data
        # for each scan number get all of the data associated with that scan
        scan_ids = []
        scan_id  = 0
        for ftf in fitacf:
#            beamTime = time2datetime(ftf)
            # print(ftf['bmnum'],beamTime,ftf['cp'],ftf['scan'])

            # These are unsupported control program IDs.
            # CPID 8600 seems to give beamnumbers of -1
            # See https://superdarn.thayer.dartmouth.edu/wg-scd.html for the SuperDARN CPID database.
            if np.abs(ftf['cp']) >= 20000 or int(np.abs(ftf['cp'])) in [8600]:
                scan_ids.append(-1)
                self.messages.append("Invalid control program ID: {!s}".format(ftf['cp']))
                continue
            
            if np.abs(ftf['scan']) == 1:
                scan_id += 1
                
            scan_ids.append(scan_id)
        scan_ids        = np.array(scan_ids)
        unq_scan_ids    = np.unique(scan_ids)
                
        for scanNum in unq_scan_ids:
            if scanNum == -1:
                continue
            
            myScan  = []
            # get data for scan number
            fitacf_inxs = np.where(scan_ids== scanNum)[0]
            
            for fit_inx in fitacf_inxs:
                myScan.append(fitacf[fit_inx])
            
            goodScan = False # This flag turns to True as soon as good data is found for the scan.
            # loop through each scan and find good data 
            for myBeam in myScan:
                beamTime = time2datetime(myBeam)
                if beamTime < eTime:
                    if stid is None or radCode is None or cp is None:
                        stid = myBeam["stid"]
                                            
                    bmnum    = myBeam["bmnum"]

                    #Calculate the field of view if it has not yet been calculated.
                    if fov == None:
                        hdw_info    = SuperDARNRadars.radars[stid].hardware_info
                        radCode     = hdw_info.abbrev
                        coords = Coords.GEOGRAPHIC
                        ranges = [0, fitacf[0]['nrang']]

                        if fovModel == "GS":
                            range_estimation = RangeEstimation.GSMR
                        elif fovModel == "HALF_SLANT":
                            range_estimation = RangeEstimation.HALF_SLANT
                        else:
                            range_estimation = RangeEstimation.SLANT_RANGE

                        beam_corners_lats, beam_corners_lons =\
                            coords(stid=myBeam['stid'],
                                rsep=myBeam["rsep"], frang=myBeam["frang"],
                                gates=ranges, date=beamTime,range_estimation=range_estimation
                                )

                        if fovModel == 'GS':
                            # The Ground Scatter Mapped Range Equation used in Bristow et al. (1994)
                            # goes imaginary at close-in ranges and therefore should return NaNs at those ranges.
                            # However, pyDARN just truncates the FOV array rather than returning NaNs.
                            # This behavior throws off keeping track of the range indices.
                            # To fix this, we add the NaNs back in at close range gates to fill out the FOV.
                            if beam_corners_lats.shape[0] != (ranges[1]+1):
                                full_shape  = (ranges[1]+1, beam_corners_lats.shape[1])
                                sInx        = (ranges[1]+1) - beam_corners_lats.shape[0]

                                tmp                 = np.zeros(full_shape)*np.nan
                                tmp[sInx:,:]        = beam_corners_lats
                                beam_corners_lats   = tmp

                                tmp                 = np.zeros(full_shape)*np.nan
                                tmp[sInx:,:]        = beam_corners_lons
                                beam_corners_lons   = tmp

                        fov = {}
                        fov["latFull"] = beam_corners_lats
                        fov["lonFull"] = beam_corners_lons

                        fov["nr_beams"] = beam_corners_lats.shape[1]-1
                        fov["nr_gates"] = beam_corners_lats.shape[0]-1

                    #Get information from each beam in the scan.
                        
                    # Save all of the radar operational parameters.
                    prm['time'].append(beamTime)
                    prm['mplgs'].append(myBeam.get("mplgs"))
                    prm['nave'].append(myBeam.get("nave"))
                    prm['noisesearch'].append(myBeam.get("noise.search"))
                    prm['scan'].append(myBeam.get("scan"))
                    prm['smsep'].append(myBeam.get("smsep"))
                    prm['mplgexs'].append(myBeam.get("mplgexs"))
                    prm['xcf'].append(myBeam.get("xcf"))
                    prm['noisesky'].append(myBeam.get("noise.sky"))
                    prm['rsep'].append(myBeam.get("rsep"))
                    prm['mppul'].append(myBeam.get("mppul"))
                    prm['inttsc'].append(myBeam.get("intt.sc"))
                    prm['frang'].append(myBeam.get("frang"))
                    prm['bmazm'].append(myBeam.get("bmazm"))
                    prm['lagfr'].append(myBeam.get("lagfr"))
                    prm['ifmode'].append(myBeam.get("ifmode"))
                    prm['noisemean'].append(myBeam.get("noise.mean"))
                    prm['tfreq'].append(myBeam.get("tfreq"))
                    prm['inttus'].append(myBeam.get("intt.us"))
                    prm['rxrise'].append(myBeam.get("rxrise"))
                    prm['mpinc'].append(myBeam.get("mpinc"))
                    prm['nrang'].append(myBeam.get("nrang"))

                    if "p_l" not in myBeam.keys():
                        continue
                        
                    fitDataList = myBeam['p_l']
                    slist       = myBeam["slist"]
                    gflag       = myBeam["gflg"]
                        
                    if len(slist) > 1:
                        for (gate,data,flag) in zip(slist,fitDataList,gflag):
                            #Get information from each gate in scan.  Skip record if the chosen ground scatter option is not met.
                            if (gscat == 1) and (flag == 0): continue
                            if (gscat == 2) and (flag == 1): continue
                            tmp = (scanNum,beamTime,bmnum,gate,data)
                            dataList.append(tmp)
                            goodScan = True
                    elif len(slist) == 1:
                        gate,data,flag = (slist[0],fitDataList[0],gflag[0])
                        #Get information from each gate in scan.  Skip record if the chosen ground scatter option is not met.
                        if (gscat == 1) and (flag == 0): continue
                        if (gscat == 2) and (flag == 1): continue
                        tmp = (scanNum,beamTime,bmnum,gate,data)
                        dataList.append(tmp)
                        goodScan = True
                    else:
                        continue
                    
            if goodScan:
                #Determine the start time for each scan and save to list.
                scanTimeList.append( (scanNum, time2datetime(fitacf[fitacf_inxs[0]])) )
         
        #Convert lists to numpy arrays.
        dataListArray   = np.array(dataList)
        
        # If no data, report and return.
        if dataListArray.size == 0:
            self.messages.append(no_data_message)
            return

        #Figure out what size arrays we need and initialize the arrays...
        nrTimes = int(np.max(dataListArray[:,scanInx]) + 1)
        nrBeams = int(np.max(dataListArray[:,beamInx]) + 1)
        nrGates = int(np.max(dataListArray[:,gateInx]) + 1)

        # Get location of radar
        radar_lat = []
        radar_lon = []
        if stid:
            radar_lat = SuperDARNRadars.radars[stid].hardware_info.geographic.lat
            radar_lon = SuperDARNRadars.radars[stid].hardware_info.geographic.lon

        if fov:
            # calculates slantRFull
            fov["slantRFull"] = np.empty((fov['nr_gates']+1,fov['nr_beams']+1))
            fov["azmFull"] = np.empty((fov['nr_gates']+1,fov['nr_beams']+1))
            for beam in range(fov['nr_beams']+1):
                for gate in range(fov['nr_gates']+1):
                    beam_lat = fov['latFull'][gate][beam]
                    beam_lon = fov['lonFull'][gate][beam]
                    fov["slantRFull"][gate][beam] = greatCircleDist(radar_lat,radar_lon,beam_lat,beam_lon) * Re #in km
                    fov['azmFull'][gate][beam] = greatCircleAzm(radar_lat,radar_lon,beam_lat,beam_lon)
                        
            #calculate azm center        
            fov["azmRCenter"] = np.empty((fov['nr_gates'],fov['nr_beams'])) 
            for gate in range(fov['nr_gates']):        
                for beam in range(fov['nr_beams']):
                    a1 = fov['azmFull'][gate][beam]
                    a2 = fov['azmFull'][gate][beam] #loop through horizontally
                    fov['azmRCenter'][gate][beam] = ((a2-a1) / 2) + a1
                         
            # calculates slantRCenter
            fov["slantRCenter"] = np.empty((fov['nr_gates'],fov['nr_beams']))      
            for beam in range(fov['nr_beams']):
                for gate in range(fov['nr_gates']):
                    l1 = fov['slantRFull'][gate][beam]
                    l2 = fov['slantRFull'][gate][beam] #loop through vertically
                    fov['slantRCenter'][gate][beam] = ((l2-l1) / 2) + l1
                    
            # calculates great circle move
            fov["latCenter"] = np.empty((fov['nr_gates'],fov['nr_beams'])) 
            fov["lonCenter"] = np.empty((fov['nr_gates'],fov['nr_beams'])) 
            for beam in range(fov['nr_beams']):
                for gate in range(fov['nr_gates']): 
                    lc = fov['slantRCenter'][gate][beam]
                    ac = fov['azmRCenter'][gate][beam]     
                    fov['latCenter'][gate][beam],fov['lonCenter'][gate][beam] = greatCircleMove(radar_lat,radar_lon,lc,ac)
                     

        # transpose the values so that the 2d arrays are arranged as (beam,gate)
        fov["latCenter"]     = np.transpose(fov["latCenter"])
        fov["lonCenter"]     = np.transpose(fov["lonCenter"])
        fov["slantRCenter"]  = np.transpose(fov["slantRCenter"])
        fov["latFull"]       = np.transpose(fov["latFull"])
        fov["lonFull"]       = np.transpose(fov["lonFull"])
        fov["slantRFull"]    = np.transpose(fov["slantRFull"])
        
        fov["beams"] = np.arange(0,fov["nr_beams"]) # get the range of beams
        #Make sure the FOV is the same size as the data array.
        if fov['nr_beams'] != nrBeams:
          fov["beams"]         = fov["beams"][0:nrBeams]
          fov["latCenter"]     = fov["latCenter"][0:nrBeams,:]
          fov["lonCenter"]     = fov["lonCenter"][0:nrBeams,:]
          fov["slantRCenter"]  = fov["slantRCenter"][0:nrBeams,:]
          fov["latFull"]       = fov["latFull"][0:nrBeams+1,:]
          fov["lonFull"]       = fov["lonFull"][0:nrBeams+1,:]
          fov["slantRFull"]    = fov["slantRFull"][0:nrBeams+1,:]
        
        fov["gates"] = np.arange(0,fov["nr_gates"]) # get the range of gates
        if fov['nr_gates'] != nrGates:
          fov["gates"]         = fov["gates"][0:nrGates]
          fov["latCenter"]     = fov["latCenter"][:,0:nrGates] #lats cal using great circle move
          fov["lonCenter"]     = fov["lonCenter"][:,0:nrGates] #lons cal using great circle move
          fov["slantRCenter"]  = fov["slantRCenter"][:,0:nrGates] # lc values
          fov["latFull"]       = fov["latFull"][:,0:nrGates+1]
          fov["lonFull"]       = fov["lonFull"][:,0:nrGates+1]
          fov["slantRFull"]    = fov["slantRFull"][:,0:nrGates+1]  #great circle dist then multiple by earth radius 
        

        fov["coords"] = fovCoords
        #Convert the dataListArray into a 3 dimensional array.
        dataArray     = np.ndarray([nrTimes,nrBeams,nrGates])
        dataArray[:]  = np.nan
        
        for inx in range(len(dataListArray)):
            dataArray[int(dataListArray[inx,scanInx]),int(dataListArray[inx,beamInx]),int(dataListArray[inx,gateInx])] = dataListArray[inx,dataInx]
        # Convert scanTimeList into a numpy array.
        # Account for the possibility that scanNums may not be continuous.
        timeListArray = [0]*nrTimes
        for scanNum, time in scanTimeList:
            timeListArray[scanNum] = time
        timeArray = np.array(timeListArray)

        # Remove skipped scans (where timestamp == 0)
        tf          = timeArray != 0
        timeArray   = timeArray[tf]
        dataArray   = dataArray[tf,:,:]

        #Make metadata block to hold information about the processing.
        metadata = {}
        metadata['dType']     = "dmap"
        metadata['stid']      = stid
        metadata['name']      = ' ' + SuperDARNRadars.radars[stid]\
                    .name
        metadata['code']      = ' ' + radCode
        metadata['fType']     = file_type
        metadata['cp']        = cp
        metadata['channel']   = channel
        metadata['sTime']     = sTime
        metadata['eTime']     = eTime
        metadata['param']     = param
        metadata['gscat']     = gscat
        # metadata['slist']     = slist # added for easy use
        metadata['elevation'] = fovElevation
        metadata['model']     = fovModel
        metadata['coords']    = fovCoords
        dataSet = 'DS000_originalFit'
        metadata['dataSetName'] = dataSet
        metadata['serial']      = 0
        comment = '['+dataSet+'] '+ 'Original Fit Data'
        #Save data to be returned as self.variables
        
        setattr(self,dataSet,musicDataObj(timeArray,dataArray,fov=fov,parent=self,comment=comment))
        newSigObj = getattr(self,dataSet)
        setattr(newSigObj,'metadata',metadata)

        #Set the new data active.
        newSigObj.setActive()

        #Make prm data part of the object.
        self.prm = prm

    def get_data_sets(self):
        """Return a sorted list of musicDataObj's contained in this musicArray.

        Returns
        -------
        dataSets : list of str
            Names of musicDataObj's contained in this musicArray.

        Written by Nathaniel A. Frissell, Fall 2013
        Updated by: Francis Tholley, 2022
        """

        attrs = dir(self)

        dataSets = []
        for item in attrs:
            if item.startswith('DS'):
                dataSets.append(item)
        dataSets.sort()
        return dataSets
