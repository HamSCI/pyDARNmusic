# Copyright (C) 2012  VT SuperDARN Lab
# Full license can be found in LICENSE.txt
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
.. module:: radDataTypes
   :synopsis: the classes needed for reading, writing, and storing fundamental
              radar data (iq,raw,fit)
.. moduleauthor:: AJ, 20130108

pydarn.sdio.radDataTypes
-------------------------

Classes
--------
radDataPtr
radBaseData
scanData
beamData
prmData
fitData
rawData
iqData
"""

# import davitpy
import logging
# from davitpy.utils import twoWayDict
alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q',
         'r','s','t','u','v','w','x','y','z']

class radDataPtr():
    """A class which contains a pipeline to a data source

    Public Attributes
    ------------------
    sTime : (datetime)
        start time of the request
    eTime : (datetime)
        end time of the request
    stid : (int)
        station id of the request
    channel : (str/NoneType)
        The 1-letter code to specify the UAF channel (not stereo),
        e.g. 'a','b',... If 'all', ALL channels were obtained.
        (default=None, meaning don't check for UAF named data files)
    bmnum : (int)
        beam number of the request
    cp : (int)
        control prog id of the request
    fType : (str)
        the file type, 'fitacf', 'fitacf3', 'rawacf', 'iqdat', 'fitex',
        'lmfit'
    fBeam : (pydarn.sdio.radDataTypes.beamData)
        the first beam of the next scan, useful for when reading into scan
        objects
    recordIndex : (dict)
        look up dictionary for file offsets for all records 
    scanStartIndex : (dict)
        look up dictionary for file offsets for scan start records

    Private Attributes
    --------------------
    ptr : (file or mongodb query object)
        the data pointer (different depending on mongodo or dmap)
    fd : (int)
        the file descriptor
    filtered : (bool)
        use Filtered datafile
    nocache : (bool)
        do not use cached files, regenerate tmp files
    src : (str)
        local or sftp

    Methods
    ----------
    open
    close
        Close file pointer
    createIndex
        Index the offsets for all records and scan boundaries
    offsetSeek
        Seek file to requested byte offset, checking to make sure it in the
        record index
    offsetTell
        Current byte offset
    rewind
        rewind file back to the beginning
    readRec
        read record at current file offset
    readScan
        read scan associated with current record
    readAll
        read all records

    Written by AJ 20130108
    """
    def __init__(self, sTime=None, radcode=None, eTime=None, stid=None,
                 channel=None, bmnum=None, cp=None, fileType=None,
                 filtered=False, src=None, fileName=None, noCache=False,
                 local_dirfmt=None, local_fnamefmt=None, local_dict=None,
                 tmpdir=None, remove=False, try_file_types=True):
        import datetime as dt
        import os,glob,string
        # from davitpy.pydarn.radar import network
        from utils.timeUtils import datetimeToEpoch
        from utils.fetchUtils import fetch_local_files
        
        import time


        self.sTime = sTime
        self.eTime = eTime
        self.stid = stid
        self.channel = channel
        self.bmnum = bmnum
        self.cp = cp
        self.fType = fileType
        self.dType = None
        self.fBeam = None
        self.recordIndex = None
        self.scanStartIndex = None
        self.__filename = fileName
        self.__filtered = filtered
        self.__nocache = noCache
        self.__src = src
        self.__fd = None
        self.__ptr =  None
        self.__records = []
        self.__read_one_rec = None
        self.__seven_tracker = False

        # check inputs
        estr = "fileType must be one of: rawacf, fitacf, fitacf3, fitex,"
        estr += " lmfit, iqdat"
        assert isinstance(self.sTime,dt.datetime), \
            logging.error('sTime must be datetime object')
        assert self.eTime == None or isinstance(self.eTime, dt.datetime), \
            logging.error('eTime must be datetime object or None')
        assert(self.channel == None or self.channel == 'all' or
               (isinstance(self.channel,str) and len(self.channel) == 1)), \
            logging.error('channel must be None or a 1-letter string')
        assert bmnum == None or isinstance(bmnum,int), \
            logging.error('bmnum must be an int or None')
        assert cp == None or isinstance(cp, int), \
            logging.error('cp must be an int or None')
        assert(fileType == 'rawacf' or fileType == 'fitacf' or
               fileType == 'fitacf3' or fileType == 'fitex' or
               fileType == 'lmfit' or fileType == 'iqdat'), \
               logging.error(estr)
        assert fileName == None or isinstance(fileName,str), \
            logging.error('fileName must be None or a string')
        assert isinstance(filtered, bool), \
            logging.error('filtered must be True of False')
        assert src == None or src == 'local' or src == 'sftp', \
            logging.error('src must be one of: None, local, sftp')

        # If channel is all, then make the channel a wildcard, then it will pull
        # in all UAF channels
        if self.channel=='all':
            channel = '.'

        if(self.eTime == None):
            self.eTime = self.sTime + dt.timedelta(days=1)

        filelist = []
        arr = [fileType]

        if try_file_types:
            all_file_types = ['fitex', 'fitacf', 'fitacf3', 'lmfit']
            try:
                all_file_types.pop(all_file_types.index(fileType))
                arr.extend(all_file_types)
            except:
                pass
        # import pdb; pdb.set_trace()
        # a temporary directory to store a temporary file
        if tmpdir is None:
            tmpdir = "/tmp/sd/"

        d = os.path.dirname(tmpdir)
        if not os.path.exists(d):
            os.makedirs(d)

        cached = False

        # FIRST, check if a specific filename was given
        if fileName != None:
            try:
                if(not os.path.isfile(fileName)):
                    estr = 'problem reading {:s} :file does '.format(fileName)
                    logging.error("{:s}not exist".format(estr))
                    return None
                outname = tmpdir + \
                          str(int(datetimeToEpoch(dt.datetime.now())))
                if(fileName.find('.bz2') != -1):
                    outname = fileName.replace('.bz2','')
                    logging.debug('bunzip2 -c '+fileName+' > '+outname+'\n')
                    os.system('bunzip2 -c '+fileName+' > '+outname)
                elif(fileName.find('.gz') != -1):
                    outname = fileName.replace('.gz','')
                    logging.debug('gunzip -c '+fileName+' > '+outname+'\n')
                    os.system('gunzip -c '+fileName+' > '+outname)
                else:
                    os.system('cp '+fileName+' '+outname)
                    logging.debug('cp '+fileName+' '+outname)
                filelist.append(outname)
                self.dType = 'dmap'
            except Exception as e:
                logging.exception(e)
                logging.exception('problem reading file', fileName)
                return None

        # Next, check for a cached file
        if fileName == None and not noCache:
            try:
                if self.channel is None:
                    gl = glob.glob("%s????????.??????.????????.??????.%s.%s" %
                                   (tmpdir, radcode, fileType))
                    for f in gl:
                        try:
                            ff = f.replace(tmpdir, '')
                            # check time span of file
                            t1 = dt.datetime(int(ff[0:4]), int(ff[4:6]),
                                             int(ff[6:8]), int(ff[9:11]),
                                             int(ff[11:13]), int(ff[13:15]))
                            t2 = dt.datetime(int(ff[16:20]), int(ff[20:22]),
                                             int(ff[22:24]), int(ff[25:27]),
                                             int(ff[27:29]), int(ff[29:31]))
                            #check if file covers our timespan
                            if t1 <= self.sTime and t2 >= self.eTime:
                                cached = True
                                filelist.append(f)
                                logging.info('Found cached file: %s' % f)
                                break
                        except Exception as e:
                            logging.exception(e)
                else:
                    gl = glob.glob("%s????????.??????.????????.??????.%s.%s.%s"
                                   % (tmpdir, radcode, self.channel, fileType))
                    for f in gl:
                        try:
                            ff = f.replace(tmpdir,'')
                            # check time span of file
                            t1 = dt.datetime(int(ff[0:4]), int(ff[4:6]),
                                             int(ff[6:8]), int(ff[9:11]),
                                             int(ff[11:13]), int(ff[13:15]))
                            t2 = dt.datetime(int(ff[16:20]), int(ff[20:22]),
                                             int(ff[22:24]), int(ff[25:27]),
                                             int(ff[27:29]), int(ff[29:31]))
                            # check if file covers our timespan
                            if t1 <= self.sTime and t2 >= self.eTime:
                                cached = True
                                filelist.append(f)
                                logging.info('Found cached file: %s' % f)
                                break
                        except Exception as e:
                            logging.exception(e)
            except Exception as e:
                logging.exception(e)
        # import pdb; pdb.set_trace()
        # Next, LOOK LOCALLY FOR FILES
        if not cached and (src == None or src == 'local') and fileName == None:
            try:
                for ftype in arr:
                    estr = "\nLooking locally for {:} files with".format(ftype)
                    estr = "{:} radcode: {:} channel: {:}".format(estr, radcode,
                                                               self.channel)
                    logging.info(estr)

                    # If the following aren't already, in the near future they
                    # will be assigned by a configuration dictionary much like
                    # matplotlib's rcsetup.py (matplotlibrc)
                    if local_dirfmt is None:
                        local_dirfmt = 'sd-data/{year}/{ftype}/{radar}/'
                        estr = 'Config entry DAVIT_LOCAL_DIRFORMAT not set,'
                        estr = '{:s} using default: '.format(estr)
                        logging.exception("{:s}{:}".format(estr,
                                                           local_dirfmt))

                    if local_dict is None:
                        local_dict = {'radar':radcode, 'ftype':ftype,
                                      'channel':channel}
                    if 'ftype' in list(local_dict.keys()):
                        local_dict['ftype'] = ftype

                    if local_fnamefmt is None:
                        local_fnamefmt = \
                          ['{date}.{hour}......{radar}.{ftype}',
                           '{date}.{hour}......{radar}.{channel}.{ftype}']
                        estr = 'Config entry DAVIT_LOCAL_FNAMEFMT not set, '
                        estr = '{:s} using default: '.format(estr)
                        logging.exception("{:s}{:}".format(estr,
                                                           local_fnamefmt))

                    outdir = tmpdir

                    # check to see if channel was specified and only use
                    # fnamefmts with channel in them
                    for f,fname in enumerate(local_fnamefmt):
                        if channel is not None and 'channel' not in fname:
                            local_fnamefmt.pop(f)
                    if len(local_fnamefmt) == 0:
                        estr = 'No file name formats containing channel exists!'
                        logging.error(estr)
                        break

                    # fetch the local files
                    temp = fetch_local_files(self.sTime, self.eTime,
                                                    local_dirfmt, local_dict,
                                                    outdir, local_fnamefmt,
                                                    remove=remove)
                    t1 = time.time()

                    # check to see if the files actually have data between stime
                    # and etime
                    valid = self.__validate_fetched(temp, self.sTime,
                                                    self.eTime)
                    t2 = time.time()
                    print(f"Runtime: {(t2-t1)}")
                    filelist = [x[0] for x in zip(temp,valid) if x[1]]
                    invalid_files = [x[0] for x in zip(temp,valid) if not x[1]]

                    if len(invalid_files) > 0:
                        for f in invalid_files:
                            logging.debug('removing invalid file: ' + f)
                            os.system('rm ' + f)

                    # If we have valid files then continue
                    if len(filelist) > 0:
                        logging.info('found ' + ftype + ' data in local files')
                        self.fType = ftype
                        self.dType = 'dmap'
                        fileType = ftype
                        break
                    else:
                        estr = "couldn't find local [{}] data".format(ftype)
                        logging.info(estr)
            except Exception as e:
                logging.exception(e)
                estr = "Unable to read local data, possible problem with "
                estr = "{:s}local_dirfmt input or rcParameter ".format(estr)
                estr = "{:s}DAVIT_LOCAL_DIRFORMAT\nWill attempt to".format(estr)
                estr = "{:s} fetch data from remote.".format(estr)
                logging.exception(estr)
                src = None

        # check if we have found files
        if len(filelist) != 0:
            # concatenate the files into a single file
            if not cached:
                logging.info('Concatenating all the files in to one')
                # choose a temp file name with time span info for cacheing
                if (self.channel is None):
                    tmpName = '%s%s.%s.%s.%s.%s.%s' % \
                              (tmpdir, self.sTime.strftime("%Y%m%d"),
                               self.sTime.strftime("%H%M%S"),
                               self.eTime.strftime("%Y%m%d"),
                               self.eTime.strftime("%H%M%S"), radcode, fileType)
                else:
                    tmpName = '%s%s.%s.%s.%s.%s.%s.%s' % \
                              (tmpdir, self.sTime.strftime("%Y%m%d"),
                               self.sTime.strftime("%H%M%S"),
                               self.eTime.strftime("%Y%m%d"),
                               self.eTime.strftime("%H%M%S"),
                               radcode, self.channel, fileType)
                logging.debug('cat ' + "".join(filelist) + ' > ' + tmpName)
                os.system('cat ' + "".join(filelist) + ' > ' + tmpName)
                for filename in filelist:
                    logging.debug('rm ' + filename)
                    os.system('rm ' + filename)
            else:
                tmpName = filelist[0]
                self.fType = fileType
                self.dType = 'dmap'

            # filter(if desired) and open the file
            if not filtered:
                self.__filename=tmpName
                self.open()
            else:
                if not fileType+'f' in tmpName:
                    try:
                        fTmpName = tmpName + 'f'
                        command = 'fitexfilter ' + tmpName + ' > ' + fTmpName
                        logging.debug("performing: {:s}".format(command))
                        os.system(command)
                    except Exception as e:
                        estr = 'problem filtering file, using unfiltered'
                        logging.warning(estr)
                        fTmpName = tmpName
                else:
                    fTmpName = tmpName
                try:
                    self.__filename=fTmpName
                    self.open()
                except Exception as e:
                    logging.exception('problem opening file')
                    logging.exception(e)

        if(self.__ptr != None):
            if(self.dType == None): self.dType = 'dmap'
        else:
            logging.error('Sorry, we could not find any data for you :(')

        os.system('rm -r /tmp/sd')


    def __repr__(self):
        myStr = 'radDataPtr: \n'
        for key,var in self.__dict__.items():
            if(isinstance(var, radDataPtr) or
               isinstance(var, type({}))):
                myStr += '%s = %s \n' % (key,'object')
            else:
                myStr += '%s = %s \n' % (key,var)
        return myStr

    def __del__(self):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        beam = self.readRec()
        if beam is None:
            raise StopIteration
        else:
            return beam

    def open(self):
        """open the associated dmap filename."""
        import os
        self.__fd = os.open(self.__filename,os.O_RDONLY)
        self.__ptr = os.fdopen(self.__fd)

    def close(self):
        """close associated dmap file."""
        import os

        if self.__ptr is not None:
            self.__ptr.close()
            self.__fd = None


    def __validate_fetched(self,filelist,stime,etime):
        """ This function checks if the files in filelist contain data
        for the start and end times (stime,etime) requested by a user.

        Parameters
        -------------
        filelist : (list)
            List of filenames 
        stime : (datetime.datetime)
            Starting time for list of filenames
        etime : (datetime.datetime)
            Ending time for list of filenames

        Returns:
        List of booleans. True if a file contains data in the time
        range (stime, etime)
        """
        # This method will need some modification for it to work with
        # file formats that are NOT DMAP (i.e. HDF5). Namely, the dmapio
        # specific code will need to be modified (readDmapRec).
        import os
        import datetime as dt
        import numpy as np
        from pydarn.io.superdarn_io import SuperDARNRead
        import pydarnio
        
        valid = []
        for f in filelist:
            logging.debug('Checking file: ' + f)
            stimes = []
            etimes = []

            # Open the file and create a file pointer
            self.__filename = f
            self.open()
            import time

            reader = pydarnio.SDarnRead(f)

            tt1 = time.time()
            records  = reader.read_fitacf()
            tt2 = time.time()
            print(f"RunTime: {tt2-tt1}")
            # Iterate through the file and grab the start time for beam
            # integration and calculate the end time from intt.sc and intt.us
            recCount = 0;
            while(recCount < len(records)):
                dfile = records[recCount]
                # read the next record from the dmap file

                if(dfile is None):
                    break
                else:
                    year = dfile["time.yr"]
                    month = dfile["time.mo"]
                    day = dfile["time.dy"]
                    hour = dfile["time.hr"]
                    minute = dfile["time.mt"]
                    second = dfile["time.sc"]
                    microsecond = dfile["time.us"]
                    temp = dt.datetime(year,month,day,hour,minute,second,microsecond)
                    stimes.append(temp)
                    sec = dfile['intt.sc'] + dfile['intt.us'] / (10. ** 6)
                    etimes.append(temp + dt.timedelta(seconds=sec))
                    
                recCount += 1
            self.__records += records
            # Close the file and clean up
            self.close()
            self.__ptr = None

            inds = np.where((np.array(stimes) >= stime) &
                            (np.array(stimes) <= etime))
            inde = np.where((np.array(etimes) >= stime) &
                            (np.array(etimes) <= etime))
            if (np.size(inds) > 0) or (np.size(inde) > 0):
                valid.append(True)
            else:
                valid.append(False)
        
        from radDataTypes import readOneRec
        # Initialize reading one rec class
        self.__read_one_rec = readOneRec(self.__records)
        return valid

    def readScan(self, firstBeam=None, useEvery=None, warnNonStandard=True,
                 showBeams=False):
        """A function to read a full scan of data from a
        :class:`pydarn.sdio.radDataTypes.radDataPtr` object.
        This function is capable of reading standard scans and extracting
        standard scans from patterned or interleaved scans (see Notes).

        Parameters
        ----------
        firstBeam : (int/NoneType)
            If manually specifying a scan pattern, will start picking beams at
            this index in the scan. Requires useEvery to also be specified.
            (default=None)
        useEvery : (int/NoneType)
            If manually specifying a scan pattern, will pick every `useEvery`
            beam. Requires firstBeam to also be specified. (default=None)
        warnNonStandard : (bool)
            If True, display a warning when auto-detecting a non-standard scan
            pattern (``firstBeam != 0`` or ``useEvery != 1``). (default=True)
        showBeams : (bool)
            `showBeams` will print the collected scan numbers. Useful for
            debugging or if you manually want to find the correct combination
            of `firstBeam` and `useEvery`. (default=False)

        Returns
        -------
        myScan : :class:`~pydarn.sdio.radDataTypes.scanData` or None
            A sequence of beams (``None`` when no more data are available)

        Notes
        -----
        For patterned scans (e.g. if beam numbers are [5, 0, 5, 1, 5, 2, ...])
        the function will try to find  a subset of the beams where beam numbers
        are increasing/decreasing by 1 throughout the scan.  Alternatively you
        can specify the pattern manually by using `firstBeam` and `useEvery`.
        You will then get the subset of the beams starting at `firstBeam`
        (which is the beam's index in the list of beams in the scan, not the
        beam number) and only including every `useEvery` beam.

        This will ignore any bmnum request in
        :func:`~pydarn.sdio.radDataRead.radDataOpen`.
        Also, if no channel was specified, it will only read channel 'a'.
        """
        from radDataTypes import scanData
        # Initialize reading one rec class
        if None in [firstBeam, useEvery] and firstBeam is not useEvery:
            estr = 'firstBeam and useEvery must both either be None or '
            raise ValueError('{:s}specified'.format(estr))

        # Save the radDataPtr's bmnum setting temporarily and set it to None
        orig_beam = self.bmnum
        self.bmnum = None

        if self.__ptr is None:
            estr = 'Self.__ptr is None.  There is probably no data available '
            logging.error('{:s}for your selected time.'.format(estr))
            self.bmnum = orig_beam
            return None

        if self.__ptr.closed:
            logging.error('Your file pointer is closed')
            self.bmnum = orig_beam
            return None

        myScan = scanData()
        recs = self.__records

        # get first beam in the scan
        myBeam = self.readRec()
        
        if myBeam is None:  # no more data
            self.bmnum = orig_beam
            return None
        while not myBeam.prm.scan:
            # continue to read until we encounter a set scan flag
            myBeam = self.readRec()
            if myBeam is None:
                # no more data
                self.bmnum = orig_beam
                return None

        # myBeam is now the first beam we encountered where scan flag is set

        myScan.append(myBeam)
        firstBeamNum = myBeam.bmnum
        # import pdb; pdb.set_trace()
        # get the rest of the beams in the scan
        while True:
            # get current offset (in case we have to revert) and next beam
            offset = myBeam.offset
            # import pdb; pdb.set_trace()
            # print("F")
            myBeam = self.readRec()
            # if(myBeam.bmnum == 7):
            #     import pdb; pdb.set_trace()

            if myBeam is None:
                # no more data
                break

            # Scan detection algorithm: We have a new scan if scan flag
            # is set AND beam number is the same as the first beam number
            # in the previous scan.
            # The latter condition is important since we can encounter
            # patterned scans with
            #   scan flags:    [1, 1, 0, 0, 0, 0, ...]
            #   beam numbers:  [5, 0, 5, 1, 5, 2, ...]
            # and we don't want to break out on the 2nd beam in this case.
            if myBeam.prm.scan and myBeam.bmnum == firstBeamNum:
                self.__read_one_rec.set_current_index()
                # print(f"BEAMNUMBER: {myBeam.bmnum} FIRSTBEAMNUM: {firstBeamNum} SCAN: {myBeam.prm.scan}")
                break
                # if start of (next) scan revert offset to start of scan and
                # break out of loop
                # pydarn.dmapio.setDmapOffset(self.__fd, offset)
            else:
                # append beam to current scan
                myScan.append(myBeam)
           
        # import pdb; pdb.set_trace()
        self.bmnum = orig_beam

        # use scan pattern from parameters if given
        if None not in [firstBeam, useEvery]:
            if showBeams:
                estr = 'Beam numbers in scan pattern for firstBeam='
                estr = '{:s}{}, useEvery={}: '.format(estr, firstBeam, useEvery)
                estr = '{:s}{}'.format(estr, [beam.bmnum for beam in myScan])
                logging.info(estr)
            # return None if scan is empty
            return myScan[firstBeam::useEvery] or None
        # import pdb; pdb.set_trace()
        # try to find the scan pattern automatically
        # print(f"AFIRSTBEAM: {firstBeam}")
        import itertools
        import numpy as np
        count = 0
        for firstBeam, useEvery in itertools.product(list(range(24)), list(range(1, 24))):
            # import pdb; pdb.set_trace()
            scan = myScan[firstBeam::useEvery]

            bmnums = [beam.bmnum for beam in scan]
            # assume correct pattern if beam numbers are increasing/decreasing
            # by one throughout the scan
            # count+=1
            # if(count == 2):
                # import pdb; pdb.set_trace()

            if np.all(np.diff(bmnums) == 1) or np.all(np.diff(bmnums) == -1):
                if showBeams or (warnNonStandard and (firstBeam != 0 or
                                                      useEvery != 1)):
                    estr = 'Auto-detected scan pattern with firstBeam='
                    estr = '{:s}{}, useEvery='.format(estr, firstBeam)
                    estr = '{:s}{} beam numbers are '.format(estr, useEvery)
                    estr = '{:s}{}'.format(estr, [beam.bmnum for beam in scan])
                    logging.info(estr)

                # return None if scan is empty
                # import pdb; pdb.set_trace()
                return scan or None
        # the only reason for not having returned yet is that the automatic
        # detection failed
        estr = 'Auto-detection of scan pattern failed, set pattern manually '
        estr = '{:s}using the firstBeam and useEvery parameters'.format(estr)
        raise ValueError(estr)

    def readRec(self):
        """A function to read a single record of radar data from a
        :class:`pydarn.sdio.radDataTypes.radDataPtr` object

        Returns
        ---------
        myBeam : (:class:`pydarn.sdio.radDataTypes.beamData`/NoneType)
        an object filled with the data we are after.  Will return None when
        finished reading.
        """
        from radDataTypes import radDataPtr, beamData, \
            fitData, prmData
        import datetime as dt

        # check input
        if(self.__ptr == None):
            logging.error('Your pointer does not point to any data')
            return None
        if self.__ptr.closed:
            logging.error('Your file pointer is closed')
            return None
        myBeam = beamData()
        # do this until we reach the requested start time
        # and have a parameter match
        records = self.__read_one_rec

        # import pdb;pdb.set_trace()

        while(1):
            dfile = records.read_one_record()
            offset = records.get_previous()
            # import pdb; pdb.set_trace()
            # index += 1
            # Handles gbr radar when scan 7 is skipped in the data retrieved
            # For example gbr data returns [7,5,7,6,7,8,7] instead of [7,5,7,6,7,7,7,8,7]
            if(dfile["bmnum"] == 8 and self.__records[offset-2]["bmnum"] != 7 and self.__seven_tracker is False):
                # import pdb;pdb.set_trace()
                dfile = records.get_record_by_index(offset-1)
                records.set_current_index(offset-1)
                self.__seven_tracker = True

            year = dfile["time.yr"]
            month = dfile["time.mo"]
            day = dfile["time.dy"]
            hour = dfile["time.hr"]
            minute = dfile["time.mt"]
            second = dfile["time.sc"]
            microsecond = dfile["time.us"]
            dfiletime = dt.datetime(year,month,day,hour,minute,second,microsecond)
            # check for valid data
            if(dfile == None or dfiletime > self.eTime):
                # if we dont have valid data, clean up, get out
                logging.info('reached end of data')
                #self.close()
                return None
            # check that we're in the time window, and that we have a
            # match for the desired params
            # if dfile['channel'] < 2: channel = 'a'  THIS CHECK IS BAD.
            # 'channel' in a dmap file specifies STEREO operation or not.
            #else: channel = alpha[dfile['channel']-1]
            if(dfiletime >= self.sTime and dfiletime <= self.eTime and
               (self.stid == None or self.stid == dfile['stid']) and
               #(self.channel == None or self.channel == channel) and
               # ASR removed because of bad check as above.
               (self.bmnum == None or self.bmnum == dfile['bmnum']) and
               (self.cp == None or self.cp == dfile['cp'])):
                # fill the beamdata object
                myBeam.updateValsFromDict(dfile)
                myBeam.recordDict = dfile
                myBeam.fType = self.fType
                myBeam.fPtr = self
                myBeam.offset = offset
                myBeam.time = dfiletime
                # file prm object
                myBeam.prm.updateValsFromDict(dfile)
                if myBeam.fType == "rawacf":
                    myBeam.rawacf.updateValsFromDict(dfile)
                if myBeam.fType == "iqdat":
                    myBeam.iqdat.updateValsFromDict(dfile)
                if(myBeam.fType == 'fitacf' or myBeam.fType == 'fitacf3' or
                   myBeam.fType == 'fitex' or myBeam.fType == 'lmfit' ):
                    myBeam.fit.updateValsFromDict(dfile)
                if type(myBeam.fit.slist) == None:
                    myBeam.fit.slist = []
                # import pdb;pdb.set_trace()
                return myBeam
                
class readOneRec():
    def __init__(self,records = None):
        self.__records = records
        self.__recLength = len(self.__records)
        self.__index = 0
        # self.currentIndex = None
        self.__previousIndex = 0
    
    def read_one_record(self):
        if(self.__records is None or self.__recLength==0 or self.__recLength == self.__index):
            return None
        if(self.__index == 0):
            record = self.__records[self.__index]
            self.__previousIndex = self.__index
            self.__index += 1
            return record
        if(self.__index > 0 and self.__recLength != self.__index):
            record = self.__records[self.__index]
            self.__previousIndex = self.__index
            self.__index += 1
            return record
    
    def get_record_by_index(self, index = 0):
        if(self.__recLength == index): return None
        return self.__records[index]
    def get_previous(self):
        return self.__previousIndex
    def set_current_index(self,index=0):
        if index == 0:
            self.__index = self.__previousIndex
        else:
            self.__index = index

class radBaseData():
    """a base class for the radar data types.  This allows for single
    definition of common routines

    Parameters
    -----------
    None

    Methods
    --------
    copyData : (func)
        Recursively copy contents into a new object
    updateValsFromDict : (func)
        converts a dict from a dmap file to radBaseData

    Written by AJ 20130108
    """

    def copyData(self,obj):
        """This method is used to recursively copy all of the contents from
        input object to self

        Parameters
        -----------
        obj : (:class:`pydarn.sdio.radDataTypes.radBaseData`)
            the object to be copied

        Returns
        --------
        Void

        Example
        ::

        myradBaseData.copyData(radBaseDataObj)

        Note
        -----
        In general, users will not need to use this.

        written by AJ, 20130402
        """
        for key, val in obj.__dict__.items():
            if isinstance(val, radBaseData):
                try:
                    getattr(self, key).copyData(val)
                except:
                    pass
            else:
                setattr(self, key, val)

    def updateValsFromDict(self, aDict):
        """A function to to fill a radar params structure with the data in a
        dictionary that is returned from the reading of a dmap file

        Parameters
        ------------
        aDict : (dict)
            The dictionary containing the radar data

        Returns
        --------
        Void

        Note
        ------
        In general, users will not need to us this.

        Written by AJ 20121130
        """
        import datetime as dt

        # iterate through prmData's attributes
        # REMOVED BY ASR on 11 SEP 2014
        # the channel attribute in fitted files (fitacf, lmfit, fitex) specifies
        # if the data came from a STEREO radar, so we shouldn't clobber the
        # value from the dmap file.
        #    elif(attr == 'channel'):
        #      if(aDict.has_key('channel')):
        #        if(isinstance(aDict.has_key('channel'), int)):
        #          if(aDict['channel'] < 2): self.channel = 'a'
        #          else: self.channel = alpha[aDict['channel']-1]
        #        else: self.channel = aDict['channel']
        #      else: self.channel = 'a'
        #      continue

        for attr, value in self.__dict__.items():
            #check for special params
            if attr == 'time':
                #convert from epoch to datetime
                if attr in aDict and isinstance(aDict[attr], float):
                    setattr(self, attr,
                            dt.datetime.utcfromtimestamp(aDict[attr]))
                continue
            elif attr == 'channel':
                if 'channel' in aDict:
                    self.channel = aDict['channel']
                continue
            elif attr == 'inttus':
                if 'intt.us' in aDict:
                    self.inttus = aDict['intt.us']
                continue
            elif attr == 'inttsc':
                if 'intt.sc' in aDict:
                    self.inttsc = aDict['intt.sc']
                continue
            elif attr == 'noisesky':
                if 'noise.sky' in aDict:
                    self.noisesky = aDict['noise.sky']
                continue
            elif attr == 'noisesearch':
                if 'noise.search' in aDict:
                    self.noisesearch = aDict['noise.search']
                continue
            elif attr == 'noisemean':
                if 'noise.mean' in aDict:
                    self.noisemean = aDict['noise.mean']
                continue
            elif attr == 'acfd' or attr == 'xcfd':
                if attr in aDict:
                    setattr(self, attr, [])
                    for i in range(self.parent.prm.nrang):
                        rec = []
                        for j in range(self.parent.prm.mplgs):
                            samp = []
                            for k in range(2):
                                aa = (i * self.parent.prm.mplgs + j) * 2 + k
                                samp.append(aDict[attr][aa])
                            rec.append(samp)
                        getattr(self, attr).append(rec)
                else:
                    setattr(self, attr, [])
                continue
            elif attr == 'mainData':
                if 'data' in aDict:
                    if(len(aDict['data']) ==
                       aDict['smpnum'] * aDict['seqnum'] * 2 * 2):
                        fac = 2
                    else:
                        fac = 1
                    setattr(self, attr, [])
                    for i in range(aDict['seqnum']):
                        rec = []
                        for j in range(aDict['smpnum']):
                            samp = []
                            for k in range(2):
                                aa = (i * fac * aDict['smpnum'] + j) * 2 + k
                                samp.append(aDict['data'][aa])
                            rec.append(samp)
                        getattr(self, attr).append(rec)
                else:
                    setattr(self, attr, [])
                continue
            elif attr == 'intData':
                if 'data' in aDict:
                    if(len(aDict['data']) ==
                       aDict['smpnum'] * aDict['seqnum'] * 2 * 2):
                        fac = 2
                    else:
                        continue
                    setattr(self, attr, [])
                    for i in range(aDict['seqnum']):
                        rec = []
                        for j in range(aDict['smpnum']):
                            samp = []
                            for k in range(2):
                                aa = ((i * fac + 1) * aDict['smpnum']
                                      + j) * 2 + k
                                samp.append(aDict['data'][aa])
                            rec.append(samp)
                        getattr(self, attr).append(rec)
                else:
                    setattr(self, attr, [])
                continue
            try:
                setattr(self, attr, aDict[attr])
            except:
                #put in a default value if not another object
                if(not isinstance(getattr(self, attr), radBaseData)):
                    setattr(self, attr, None)

  #def __repr__(self):
    #myStr = ''
    #for key,var in self.__dict__.iteritems():
      #if(isinstance(var,radBaseData) and key != 'parent'):
        #print key
        #myStr += key+'\n'
        #myStr += str(var)
      #else:
        #myStr += key+' = '+str(var)+'\n'
    #return myStr

class scanData(list):
    """a class to contain a radar scan.  Extends list.
    Just a list of :class:`pydarn.sdio.radDataTypes.beamData` objects

    Attributes
    ----------
    None

    Example
    --------
    ::

    myBeam = pydarn.sdio.scanData()

    Written by AJ 20121130
    """

    def __init__(self):
        pass

class beamData(radBaseData):
    """a class to contain the data from a radar beam sounding,
    extends class :class:`pydarn.sdio.radDataTypes.radBaseData`

    Attributes
    -----------
    cp : (int)
        radar control program id number
    stid : (int)
        radar station id number
    time : (datetime)
        timestamp of beam sounding
    channel (int)
        radar operating channel defined by STEREO operations, eg 0, 1, 2.
        Zero is for non-stereo operations and 1 & 2 are for STEREO operations
        of A & B channels
    bmnum : (int)
        beam number
    prm : (pydarn.sdio.radDataTypes.prmData)
        operating params
    fit : (pydarn.sdio.radDataTypes.fitData)
        fitted params
    rawacf : (pydarn.sdio.radDataTypes.rawData)
        rawacf data
    iqdat : (pydarn.sdio.radDataTypes.iqData)
        iqdat data
    fType : (str)
        the file type, 'fitacf', 'fitacf3' 'rawacf',
        'iqdat', 'fitex', 'lmfit'

    Example
    --------
    ::

    myBeam = pydarn.sdio.radBeam()

    Written by AJ 20121130
    """
    def __init__(self, beamDict=None, myBeam=None, proctype=None):
        #initialize the attr values
        self.cp = None
        self.stid = None
        self.time = None
        self.bmnum = None
        self.channel = None
        self.exflg = None
        self.lmflg = None
        self.acflg = None
        self.rawflg = None
        self.iqflg = None
        self.fitex = None
        self.fitacf = None
        self.lmfit= None
        self.fit = fitData()
        self.rawacf = None
        self.prm = prmData()
        self.iqdat = None
        self.recordDict = None
        self.fType = None
        self.offset = None
        self.fPtr = None
        #if we are intializing from an object, do that
        if(beamDict != None):
            self.updateValsFromDict(beamDict)

    def __repr__(self):
        import datetime as dt
        myStr = 'Beam record FROM: ' + str(self.time) + '\n'
        for key,var in self.__dict__.items():
            if(isinstance(var, radBaseData) or isinstance(var, radDataPtr) or
               isinstance(var, type({}))):
                myStr += '%s  = %s \n' % (key, 'object')
            else:
                myStr += '%s  = %s \n' % (key, var)
        return myStr


class prmData(radBaseData):
    """A class to represent radar operating parameters, extends
    :class:`pydarn.sdio.radDataTypes.radBaseData `

    Attributes
    -----------
    nave : (int)
        number of averages
    lagfr : (int)
        lag to first range in us
    smsep : (int)
        sample separation in us
    bmazm : (float)
        beam azimuth
    scan : (int)
        new scan flag
    rxrise : (int)
        receiver rise time
    inttsc : (int)
        integeration time (sec)
    inttus : (int)
        integration time (us)
    mpinc : (int)
        multi pulse increment (tau, basic lag time) in us
    mppul : (int)
        number of pulses
    mplgs : (int)
        number of lags
    mplgexs : (int)
        number of lags (tauscan)
    nrang : (int)
        number of range gates
    frang : (int)
        first range gate (km)
    rsep : (int)
        range gate separation in km
    xcf : (int)
        xcf flag
    tfreq : (int)
        transmit freq in kHz
    txpl : (int)
        transmit pulse length in us
    ifmode : (int)
        if mode flag
    ptab : (mppul length list)
        pulse table
    ltab : (mplgs x 2 length list)
        lag table
    noisemean : (float)
        mean noise level
    noisesky : (float)
        sky noise level
    noisesearch : (float)
        freq search noise level

    Written by AJ 20121130
    """

    # initialize the struct
    def __init__(self, prmDict=None, myPrm=None):
        # set default values
        self.nave = None        #number of averages
        self.lagfr = None       #lag to first range in us
        self.smsep = None       #sample separation in us
        self.bmazm = None       #beam azimuth
        self.scan = None        #new scan flag
        self.rxrise = None      #receiver rise time
        self.inttsc = None      #integeration time (sec)
        self.inttus = None      #integration time (us - microsec)
        self.mpinc = None       #multipulse increment (ms) (tau, basic lag time)
        self.mppul = None       #number of pulses
        self.mplgs = None       #number of lags
        self.mplgexs = None     #number of lags (tauscan)
        self.nrang = None       #number of range gates
        self.frang = None       #first range gate (km)
        self.rsep = None        #range gate separation in km
        self.xcf = None         #xcf flag
        self.tfreq = None       #transmit freq in kHz
        self.txpl = None       #transmit freq in kHz
        self.ifmode = None      #if mode flag
        self.ptab = None        #pulse table
        self.ltab = None        #lag table
        self.noisemean = None   #mean noise level
        self.noisesky = None    #sky noise level
        self.noisesearch = None #freq search noise level

        #if we are copying a structure, do that
        if(prmDict != None):
            self.updateValsFromDict(prmDict)

    def __repr__(self):
        import datetime as dt
        myStr = 'Prm data: \n'
        for key,var in self.__dict__.items():
            myStr += '%s  = %s \n' % (key, var)
        return myStr

class fitData(radBaseData):
    """a class to contain the fitted params of a radar beam sounding,
    extends :class:`pydarn.sdio.radDataTypes.radBaseData`

    Attributes
    ------------
    pwr0 : (prm.nrang length list)
        lag 0 power
    slist : (npnts length list)
        list of range gates with backscatter
    npnts (int)
        number of range gates with scatter
    nlag : (npnts length list)
        number of good lags
    qflg : (npnts length list)
        quality flag
    gflg : (npnts length list)
        ground scatter flag
    p_l : (npnts length list)
        lambda power
    p_l_e : (npnts length list)
        lambda power error
    p_s : (npnts length list)
        sigma power
    p_s_e : (npnts length list)
        sigma power error
    v : (npnts length list)
        velocity
    v_e : (npnts length list)
        velocity error
    w_l : (npnts length list)
        lambda spectral width
    w_l_e : (npnts length list)
        lambda width error
    w_s : (npnts length list)
        sigma spectral width
    w_s_e : (npnts length list)
        sigma width error
    phi0 : (npnts length list)
        phi 0
    phi0_e : (npnts length list)
        phi 0 error
    elv : (npnts length list)
        elevation angle

    Example
    ---------
    ::

    myFit = pydarn.sdio.fitData()

    Written by AJ 20121130
    """
    # initialize the struct
    def __init__(self, fitDict=None, myFit=None):
        self.pwr0 = None      #lag 0 power
        self.slist = None     # list of range gates with backscatter
        self.npnts = None     #number of range gates with scatter
        self.nlag = None      #number of good lags
        self.qflg = None      #quality flag
        self.gflg = None      #ground scatter flag
        self.p_l = None       #lambda power
        self.p_l_e = None     #lambda power error
        self.p_s = None       #sigma power
        self.p_s_e = None     #sigma power error
        self.v = None         #velocity
        self.v_e = None       #velocity error
        self.w_l = None       #lambda spectral width
        self.w_l_e = None     #lambda width error
        self.w_s = None       #sigma spectral width
        self.w_s_e = None     #sigma width error
        self.phi0 = None      #phi 0
        self.phi0_e = None    #phi 0 error
        self.elv = None       #elevation angle

        if(fitDict != None):
            self.updateValsFromDict(fitDict)

    def __repr__(self):
        import datetime as dt
        myStr = 'Fit data: \n'
        for key,var in self.__dict__.items():
            myStr += '%s = %s \n' % (key, var)
        return myStr