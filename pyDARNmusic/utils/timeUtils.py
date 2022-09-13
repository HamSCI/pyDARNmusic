# -*- coding: utf-8 -*-
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

"""timeUtils module

A module for manipulating time representations

Functions
-----------------------------------
JulianDayFromDate    creates a Julian Day from a 'datetime-like' object Adapted from mpl_toolkits.basemap
epem                 meridian and the meridian containing the subsolar point Adapted from mpl_toolkits.basemap
daynight_terminator  Calculates the latitude, Greenwich Hour Angle, and solar declination from a given latitude and longitude
dateToDecYear        datetime to decimal
dateToYyyymmdd       datetime to string
datetimeToEpoch      datetime to epoch
julToDatetime        julian to datetime
parseDate            parse yr, mn, day
parseTime            parse hr, min, sec
timeYrsecToDate      time to datetime
yyyymmddToDate       string to date
-----------------------------------

"""
import logging
import numpy as np

def JulianDayFromDate(date,calendar='standard'):
    """
    creates a Julian Day from a 'datetime-like' object.  Returns the fractional
    Julian Day (resolution 1 second).
    if calendar='standard' or 'gregorian' (default), Julian day follows Julian 
    Calendar on and before 1582-10-5, Gregorian calendar after 1582-10-15.
    if calendar='proleptic_gregorian', Julian Day follows gregorian calendar.
    if calendar='julian', Julian Day follows julian calendar.
    Algorithm:
    Meeus, Jean (1998) Astronomical Algorithms (2nd Edition). Willmann-Bell,
    Virginia. p. 63
    Adapted from mpl_toolkits.basemap. by Francis Tholley, 2022
    """
    # based on redate.py by David Finlayson.
    year=date.year; month=date.month; day=date.day
    hour=date.hour; minute=date.minute; second=date.second
    # Convert time to fractions of a day
    day = day + hour/24.0 + minute/1440.0 + second/86400.0
    # Start Meeus algorithm (variables are in his notation)
    if (month < 3):
        month = month + 12
        year = year - 1
    A = int(year/100)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + \
         day - 1524.5
    # optionally adjust the jd for the switch from 
    # the Julian to Gregorian Calendar
    # here assumed to have occurred the day after 1582 October 4
    if calendar in ['standard','gregorian']:
        if jd >= 2299170.5:
            # 1582 October 15 (Gregorian Calendar)
            B = 2 - A + int(A/4)
        elif jd < 2299160.5:
            # 1582 October 5 (Julian Calendar)
            B = 0
        else:
            raise ValueError('impossible date (falls in gap between end of Julian calendar and beginning of Gregorian calendar')
    elif calendar == 'proleptic_gregorian':
        B = 2 - A + int(A/4)
    elif calendar == 'julian':
        B = 0
    else:
        raise ValueError('unknown calendar, must be one of julian,standard,gregorian,proleptic_gregorian, got %s' % calendar)
    # adjust for Julian calendar if necessary
    jd = jd + B
    return jd 

def epem(date):
    """
    input: date - datetime object (assumed UTC)
    ouput: gha - Greenwich hour angle, the angle between the Greenwich
           meridian and the meridian containing the subsolar point.
           dec - solar declination.
    Adapted from mpl_toolkits.basemap. by Francis Tholley, 2022
    """
    dg2rad = np.pi/180.
    rad2dg = 1./dg2rad
    # compute julian day from UTC datetime object.
    # datetime objects use proleptic gregorian calendar.
    jday = JulianDayFromDate(date,calendar='proleptic_gregorian')
    jd = np.floor(jday) # truncate to integer.
    # utc hour.
    ut = date.hour + date.minute/60. + date.second/3600.
    # calculate number of centuries from J2000
    t = (jd + (ut/24.) - 2451545.0) / 36525.
    # mean longitude corrected for aberration
    l = (280.460 + 36000.770 * t) % 360
    # mean anomaly
    g = 357.528 + 35999.050 * t
    # ecliptic longitude
    lm = l + 1.915 * np.sin(g*dg2rad) + 0.020 * np.sin(2*g*dg2rad)
    # obliquity of the ecliptic
    ep = 23.4393 - 0.01300 * t
    # equation of time
    eqtime = -1.915*np.sin(g*dg2rad) - 0.020*np.sin(2*g*dg2rad) \
            + 2.466*np.sin(2*lm*dg2rad) - 0.053*np.sin(4*lm*dg2rad)
    # Greenwich hour angle
    gha = 15*ut - 180 + eqtime
    # declination of sun
    dec = np.arcsin(np.sin(ep*dg2rad) * np.sin(lm*dg2rad)) * rad2dg
    return gha, dec


def daynight_terminator(date, lons):
    """Calculates the latitude, Greenwich Hour Angle, and solar 
    declination from a given latitude and longitude.

    This routine is used by musicRTI for terminator calculations.

    Parameters
    ----------
    date : datetime.datetime
        UT date and time of terminator calculation.
    lons : np.array
        Longitudes of which to calculate the terminator.

    Returns
    -------
    lats : np.array
        Latitudes of solar terminator.
    tau : np.array
        Greenwhich Hour Angle.
    dec : np.array
        Solar declination.

    Notes
    -----
    Adapted from mpl_toolkits.basemap.solar by Nathaniel A. Frissell, Fall 2013

    """
    # import mpl_toolkits.basemap.solar as solar

    dg2rad = np.pi/180.
    # compute greenwich hour angle and solar declination
    # from datetime object (assumed UTC).
    tau, dec = epem(date)
    # compute day/night terminator from hour angle, declination.
    longitude = lons + tau
    lats = np.arctan(-np.cos(longitude*dg2rad)/np.tan(dec*dg2rad))/dg2rad
    return lats,tau,dec

def dateToYyyymmdd(my_date):
    """Takes a python datetime object and returns a string in yyyymmdd format

    Parameters
    ----------
    my_date : datetime
        a python datetime object

    Returns
    -------
    date_str : str
        a string in yyyymmdd format

    Example
    -------
        import datetime as dt
        date_str = utils.timeUtils.dateToYyyymmdd(dt.datetime(2012,7,10))

    Written by AJ 20120718
    Modified by ASR 20151120

    """
    from datetime import datetime

    assert(isinstance(my_date, datetime)), logging.error(
        'input must be type datetime')

    return my_date.strftime('%Y%m%d')


def yyyymmddToDate(date_str):
    """takes a string in yyyymmdd format and returns a python date object

    Parameters
    ----------
    date_str : str
        a string in yyyymmdd format

    Returns
    -------
    my_date : datetime
        a python datetime object

    Example
    -------
        my_date = utils.timeUtils.yyyymmddToDate('20120710')

    Written by AJ 20120718
    Modified by ASR 20151120

    """
    from datetime import datetime

    assert(isinstance(date_str, str)), logging.error(
        'input must be of type str')
    assert(len(date_str) == 8), logging.error(
        'input must be of yyyymmdd format')

    return datetime.strptime(date_str, '%Y%m%d')


def timeYrsecToDate(yrsec, year):
    """Converts time expressed in seconds from start of year to a python
    datetime object

    Parameters
    ----------
    yrsec : int
        seconds since start of year
    year : int
        year in YYYY

    Returns
    -------
    my_date : datetime
        a python datetime object.

    Example
        my_date = utils.timeUtils.timeYrsecToDate(1205304,2012)

    Written by Sebastien, Jul. 2012
    Modified by ASR 20151120

    """
    from datetime import datetime
    from datetime import timedelta

    assert(isinstance(yrsec, int)), logging.erorr(
        'yrsec must be of type int')
    assert(isinstance(year, int)), logging.error(
        'year must be of type int')

    return datetime(year, 1, 1) + timedelta(seconds=yrsec)


# def julToDatetime(ndarray):
#     """Convert a julian date to a datetime object.

#     Parameters
#     ----------
#     ndarray : float or list
#         single float64 or a numpy array of Julian Dates.

#     Returns
#     -------
#     dt : list
#         list of datetime objects

#     Example
#     -------
#         myDateList = utils.timeUtils.julToDatetime(2456118.5)

#     Created by Nathaniel Frissell 20120810 

#     """
#     import datetime
#     import dateutil.parser
#     import numpy
#     import spacepy.time as spt

#     t = spt.Ticktock(ndarray, 'JD')

#     dt = list()
#     for iso in t.ISO:
#         dt.append(dateutil.parser.parse(iso))

#     return dt


def datetimeToEpoch(my_date):
    """reads in a datetime and outputs the equivalent epoch time

    Parameters
    ----------
    my_date : datetime
        a datetime object

    Returns
    -------
    my_epoch : float or list of floats
        an epoch time equal to the datetime object

    Example
    -------
        import datetime as dt
        my_epoch = utils.timeUtils.datetimeToEpoch(dt.datetime(2012,7,10))

    Written by AJ 20120914
    Modified by Nathaniel Frissell 20130729 - Added list support.
    Modified by AJ 20130925 - fixed local time bug with conversion
    Modified by Nathaniel Frissell 20130930 - Fixed list support.
    Modified by ASR 20151120

    """
    from datetime import datetime
    import calendar
    import numpy as np

    assert(isinstance(my_date, (list, np.ndarray, datetime))), logging.error(
        'input must be of type datetime or list of datetimes')
    if isinstance(my_date, (list,np.ndarray)):
        for dt in my_date:
            assert(isinstance(dt, datetime)), logging.error(
                'each member of my_date list must be of type datetime')

    if isinstance(my_date, datetime):
        unx = calendar.timegm(my_date.timetuple()) + my_date.microsecond / 1e6
    else:
        unx = [calendar.timegm(dt.timetuple()) +
               dt.microsecond / 1e6 for dt in my_date]

    return unx


def dateToDecYear(date):
    """Convert (`datetime <http://tinyurl.com/bl352yx>`_) object to decimal
    year

    Parameters
    ----------
    date : datetime
        date and time

    Returns
    -------
    dyear : float
        decimal year

    Example
    -------
        import datetime as dt
        decYr = utils.timeUtils.dateToDecYear(dt.datetime(2012,7,10))

    written by Sebastien, 2013-02

    """
    from datetime import datetime as dt
    import time

    # returns seconds since epoch
    sE = lambda date: time.mktime(date.timetuple())

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year + 1, month=1, day=1)

    yearElapsed = sE(date) - sE(startOfThisYear)
    yearDuration = sE(startOfNextYear) - sE(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction


def parseDate(date):
    """ Parse YYYYMMDD dates in YYYY, MM, DD and vice versa

    Parameters
    ----------
    date : str or list
        experiment date in YYYYMMDD or [YYYY,MM,DD]

    Returns
    -------
    tdate : list or int
        experiment date in [YYYY,MM,DD] or YYYYMMDD

    Example
    -------
        tlist = utils.timeUtils.parseDate('20120710')
        ttime = utils.timeUtils.parseDate([2012,07,10])

    Created by Sebastien

    """
    # transform date into an array for testing
    if not isinstance(date, list):
        date = [date]

    # make sure we are getting integers
    for id in range(len(date)):
        if isinstance(date[id], str):
            date[id] = int(date[id])

    # parse date one way or another
    if len(date) == 3:
        tdate = date[0] * 10000 + date[1] * 100 + date[2]
    elif len(date) == 1:
        tdate = [date[0] / 10000, date[0] / 100 - date[0] /
                 10000 * 100, date[0] - date[0] / 100 * 100]
    else:
        logging.error('Invalid date format: ', date)
        return

    return tdate


def parseTime(time):
    """ Parse HHMM or HHMMSS dates in HH, MM, SS and vice versa

    Parameters
    ----------
    time : str or list
        time in HHMM or HHMMSS OR [HH,MM] or [HH,MM,SS]

    Returns
    -------
    ttime : list or int
        time in [HH,MM] or [HH,MM,SS] OR HHMM or HHMMSS

    Example
    -------
        tlist = utils.timeUtils.parseDate('065022')
        tstr = utils.timeUtils.parseDate([6,50,22])

    Created by Sebastien

    """
    # transform time into an array for testing
    if not isinstance(time, list):
        time = [time]

    # make sure we are getting integers
    for it in range(len(time)):
        if isinstance(time[it], str):
            time[it] = int(time[it])

    # parse time one way or another
    if len(time) == 3:
        ttime = time[0] * 10000 + time[1] * 100 + time[2]
    elif len(time) == 2:
        ttime = time[0] * 100 + time[1]
    elif len(time) == 1 and len(str(time[0])) > 4 and len(str(time[0])) <= 6:
        ttime = [time[0] / 10000, time[0] / 100 - time[0] /
                 10000 * 100, time[0] - time[0] / 100 * 100]
    elif len(time) == 1 and len(str(time[0])) >= 1 and len(str(time[0])) <= 4:
        ttime = [time[0] / 100, time[0] - time[0] / 100 * 100]
    else:
        logging.error('Invalid time format: ', time)
        return

    return ttime
