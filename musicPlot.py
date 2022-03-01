import numpy as np
import scipy as sp
import datetime

from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon
from matplotlib import dates as md
import matplotlib

# from mpl_toolkits.basemap import Basemap

# from davitpy import utils
# from davitpy.pydarn.radar.radUtils import getParamDict

from music import getDataSet
from radUtils import getParamDict
from plotUtils import genCmap
import logging

#Global Figure Size
figsize=(20,10)


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
    return gha, 


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


class musicRTI(object):
    """Class to create an RTI plot using a pydarn.proc.music.musicArray object as the data source.

    Parameters
    ----------
    dataObj : pydarn.proc.music.musicArray
        musicArray object
    dataSet : Optional[str]
        which dataSet in the musicArray object to plot
    beam : Optional[int]
        Beam number to plot.
    xlim : Optoinal[None or 2-element iterable of datetime.datetime]
        Limits for x-axis.
    ylim : Optional[None or 2-element iterable of floats]
        Limits for y-axis.
    axis : Optional[None or matplotlib.figure.axis]
        Matplotlib axis on which to plot.  If None, a new figure and axis will be created.
    scale : Optional[None or 2-Element iterable]
        Colorbar scale.  If None, the default scale for the current SuperDARN parameter will be used.
    plotZeros : Optional[bool]
        If True, plot data cells that are identically zero.
    max_sounding_time : Optional[None or datetime.timedelta]
        Do not allow data to be plotted for longer than this duration.
    xBoundaryLimits: Optional[None or 2-element iterable of datetime.datetime]
        Mark a region of times on the RTI plot.  A green dashed vertical line will be plotted
        at each of the boundary times.  The region of time outside of the boundary will be shaded gray.
        If set to None, this will automatically be set to the timeLimits set in the metadata, if they exist.
    yBoundaryLimits : Optional[None or 2-element iterable of floats]
        Mark a region of range on the RTI plot.  A green dashed horizontal line will be plotted
        at each of the boundary ranges.  The region of time outside of the boundary will be shaded gray.
        If set to None, this will automatically be set to the gateLimits set in the metadata, if they exist.
    yticks : Optional[list]
        Where to put the ticks on the y-axis.
    ytick_lat_format : Optional[str]
         %-style string format code for latitude y-tick labels
    autoScale : Optional[bool]
        If True, automatically scale the color bar for good data visualization. Keyword scale must be None when using autoScale.
        ax.set_xlim(xlim)
    plotTerminator : Optional[bool]
        If True, overlay day/night terminator on the RTI plot.  Every cell is evaluated for day/night and shaded accordingly.  Therefore,
        terminator resolution will match the resolution of the RTI plot data.
    axvlines : Optional[None or list of datetime.datetime]
        Dashed vertical lines will be drawn at each specified datetime.datetime.
    axvline_color : Optional[str]
        Matplotlib color code specifying color of the axvlines.
    secondary_coords : Optional[str]
        Secondary coordate system for RTI plot y-axis ('lat' or 'range')
    plot_info : Optional[bool]
        If True, plot frequency/noise plots
    plot_title : Optional[bool]
        If True, plot the title information
    plot_range_limits_label : Optoinal[bool]
        If True, plot the label corresponding to the range limits on the right-hand y-axis.
    cmap_handling : Optional[str]
        'superdarn' to use SuperDARN-style colorbars, 'matplotlib' for direct use of matplotlib's colorbars.
        'matplotlib' is recommended when using custom scales and the 'superdarn' mode is not providing a desirable result.
    plot_cbar : Optional[bool]
        If True, plot the color bar.
    cbar_ticks : Optional[list]
        Where to put the ticks on the color bar.
    cbar_shrink : Optional[float]
        fraction by which to shrink the colorbar
    cbar_fraction : Optional[float]
        fraction of original axes to use for colorbar
    cbar_gstext_offset : Optional[float]
        y-offset from colorbar of "Ground Scatter Only" text
    cbar_gstext_fontsize : Optional[float]
        fontsize of "Ground Scatter Only" text
    model_text_size : Optional[int]
        fontsize of model and coordinate indicator text
    **kwArgs :
        Keyword Arguments

    Attributes
    ----------
    cbar_info : list


    Written by Nathaniel A. Frissell, Fall 2013

    """
    def __init__(self,dataObject,
        dataSet                 = 'active',
        beam                    = 7,
        coords                  = 'gate',
        xlim                    = None,
        ylim                    = None,
        axis                    = None,
        scale                   = None,
        plotZeros               = False, 
        max_sounding_time       = datetime.timedelta(minutes=4),
        xBoundaryLimits         = None,
        yBoundaryLimits         = None,
        yticks                  = None,
        ytick_lat_format        = '.0f',
        autoScale               = False,
        plotTerminator          = True,
        axvlines                = None, 
        axvline_color           = '0.25',
        secondary_coords        = 'lat',
        plot_info               = True,
        plot_title              = True,
        plot_range_limits_label = True,
        cmap_handling           = 'superdarn',
        cmap                    = None,
        bounds                  = None,
        norm                    = None,
        plot_cbar               = True,
        cbar_ticks              = None,
        cbar_shrink             = 1.0,
        cbar_fraction           = 0.15,
        cbar_gstext_offset      = -0.075,
        cbar_gstext_fontsize    = None,
        model_text_size         = 'small',
        y_labelpad              = None,
        **kwArgs):

        from scipy import stats
        from rti import plot_freq,plot_nave,plot_skynoise,plot_searchnoise

        if axis is None:
            from matplotlib import pyplot as plt
            fig   = plt.figure(figsize=figsize)

        # Make some variables easier to get to...
        currentData = getDataSet(dataObject,dataSet)
        metadata    = currentData.metadata
        latFull     = currentData.fov.latFull
        lonFull     = currentData.fov.lonFull
        latCenter   = currentData.fov.latCenter
        lonCenter   = currentData.fov.lonCenter
        time        = currentData.time
        beamInx     = np.where(currentData.fov.beams == beam)[0]
        radar_lats  = latCenter[beamInx,:]
        nrTimes, nrBeams, nrGates = np.shape(currentData.data)

        # Calculate terminator. ########################################################
        if plotTerminator:
            daylight = np.ones([nrTimes,nrGates],np.bool)
            for tm_inx in range(nrTimes):
                tm                  = time[tm_inx]
                term_lons           = lonCenter[beamInx,:]
                term_lats,tau,dec   = daynight_terminator(tm,term_lons)

                if dec > 0: # NH Summer
                    day_inx = np.where(radar_lats < term_lats)[1]
                else:
                    day_inx = np.where(radar_lats > term_lats)[1]

                if day_inx.size != 0:
                    daylight[tm_inx,day_inx] = False

        # Translate parameter information from short to long form.
        paramDict = getParamDict(metadata['param'])
        if 'label' in paramDict:
            param     = paramDict['param']
            cbarLabel = paramDict['label']
        else:
            param = 'width' # Set param = 'width' at this point just to not screw up the colorbar function.
            cbarLabel = metadata['param']

        # Set colorbar scale if not explicitly defined.
        if(scale is None):
            if autoScale:
                sd          = stats.nanstd(np.abs(currentData.data),axis=None)
                mean        = stats.nanmean(np.abs(currentData.data),axis=None)
                scMax       = np.ceil(mean + 1.*sd)
                if np.min(currentData.data) < 0:
                    scale   = scMax*np.array([-1.,1.])
                else:
                    scale   = scMax*np.array([0.,1.])
            else:
                if 'range' in paramDict:
                    scale = paramDict['range']
                else:
                    scale = [-200,200]

        # See if an axis is provided... if not, set one up!
        if axis is None:
            axis    = fig.add_subplot(111)
        else:
            fig   = axis.get_figure()

        if np.size(beamInx) == 0:
            beamInx = 0
            beam    = currentData.fov.beams[0]

        # Plot the SuperDARN data!
        verts = []
        scan  = []
        data  = np.squeeze(currentData.data[:,beamInx,:])

#        The coords keyword needs to be tested better.  For now, just allow 'gate' only.
#        Even in 'gate' mode, the geographic latitudes are plotted along with gate.
#        if coords is None and metadata.has_key('coords'):
#            coords      = metadata['coords']
#
        if coords not in ['gate','range']:
            logging.warning('Coords "%s" not supported for RTI plots.  Using "gate".' % coords)
            coords = 'gate'

        if coords == 'gate':
            rnge  = currentData.fov.gates
        elif coords == 'range':
            rnge  = currentData.fov.slantRFull[beam,:]

        xvec  = [matplotlib.dates.date2num(x) for x in currentData.time]
        for tm in range(nrTimes-1):
            for rg in range(nrGates-1):
                if np.isnan(data[tm,rg]): continue
                if data[tm,rg] == 0 and not plotZeros: continue
                if max_sounding_time is not None:
                    if (currentData.time[tm+1] - currentData.time[tm+0]) > max_sounding_time: continue
                scan.append(data[tm,rg])

                x1,y1 = xvec[tm+0],rnge[rg+0]
                x2,y2 = xvec[tm+1],rnge[rg+0]
                x3,y3 = xvec[tm+1],rnge[rg+1]
                x4,y4 = xvec[tm+0],rnge[rg+1]
                verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))

        if (cmap_handling == 'matplotlib') or autoScale:
            if cmap is None:
                cmap = matplotlib.cm.jet
            if bounds is None:
                bounds  = np.linspace(scale[0],scale[1],256)
            if norm is None:
                norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
        elif cmap_handling == 'superdarn':
            colors  = 'lasse'
            cmap,norm,bounds = genCmap(param,scale,colors=colors)

        pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm,zorder=99)
        pcoll.set_array(np.array(scan))
        axis.add_collection(pcoll,autolim=False)

        # Plot the terminator! #########################################################
        if plotTerminator:
#            print 'Terminator functionality is disabled until further testing is completed.'
            term_verts = []
            term_scan  = []

            rnge  = currentData.fov.gates
            xvec  = [matplotlib.dates.date2num(x) for x in currentData.time]
            for tm in range(nrTimes-1):
                for rg in range(nrGates-1):
                    if daylight[tm,rg]: continue
                    term_scan.append(1)

                    x1,y1 = xvec[tm+0],rnge[rg+0]
                    x2,y2 = xvec[tm+1],rnge[rg+0]
                    x3,y3 = xvec[tm+1],rnge[rg+1]
                    x4,y4 = xvec[tm+0],rnge[rg+1]
                    term_verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))

            term_pcoll = PolyCollection(np.array(term_verts),facecolors='0.45',linewidth=0,zorder=99,alpha=0.25)
            axis.add_collection(term_pcoll,autolim=False)
        ################################################################################

        if axvlines is not None:
            for line in axvlines:
                axis.axvline(line,color=axvline_color,ls='--')

        if xlim is None:
            xlim = (np.min(time),np.max(time))
        axis.set_xlim(xlim)

        axis.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))
        axis.set_xlabel('Time [UT]')

        if ylim is None:
            ylim = (np.min(rnge),np.max(rnge))
        axis.set_ylim(ylim)

        if yticks is not None:
            axis.set_yticks(yticks)

        # Y-axis labeling ##############################################################
        if coords == 'gate':
            if secondary_coords:
                if secondary_coords == 'range':
                    if metadata['model'] == 'IS':
                        axis.set_ylabel('Range Gate\nSlant Range [km]',labelpad=y_labelpad)
                    elif metadata['model'] == 'GS':
                        axis.set_ylabel('Range Gate\nGS Mapped Range [km]',labelpad=y_labelpad)
                else:
                    geo_mag = 'Geographic' if currentData.fov.coords == 'geo' else 'Magnetic'
                    if metadata['model'] == 'IS':
                        axis.set_ylabel('Range Gate\n%s Latitude' % geo_mag,labelpad=y_labelpad)
                    elif metadata['model'] == 'GS':
                        axis.set_ylabel('Range Gate\nGS Mapped %s Latitude' % geo_mag,labelpad=y_labelpad)

                yticks  = axis.get_yticks()
                ytick_str    = []
                for tck in yticks:
                    txt = []
                    txt.append('%d' % tck)

                    rg_inx = np.where(tck == currentData.fov.gates)[0]
                    if np.size(rg_inx) != 0:
                        if secondary_coords == 'range':
                            rang = currentData.fov.slantRCenter[beamInx,rg_inx]
                            if np.isfinite(rang): 
                                txt.append('%d' % rang)
                            else:
                                txt.append('')
                        else:
                            lat = currentData.fov.latCenter[beamInx,rg_inx]
                            if np.isfinite(lat): 
                                txt.append(('%'+ytick_lat_format+'$^o$') % lat)
                            else:
                                txt.append('')
                    txt = '\n'.join(txt)
                    ytick_str.append(txt)
                axis.set_yticklabels(ytick_str,rotation=90,ma='center')
            else:
                axis.set_ylabel('Range Gate',labelpad=y_labelpad)
        elif coords == 'range':
            if secondary_coords == 'lat':
                # Use linear interpolation to get the latitude associated with a particular range.
                # Make sure we only include finite values in the interpolation function.
                finite_inx  = np.where(np.isfinite(currentData.fov.latCenter[beam,:]))[0]
                tmp_ranges  = currentData.fov.slantRCenter[beam,:][finite_inx]
                tmp_lats    = currentData.fov.latCenter[beam,:][finite_inx]
                tmp_fn      = sp.interpolate.interp1d(tmp_ranges,tmp_lats)

                yticks  = axis.get_yticks()
                ytick_str    = []
                for tck in yticks:
                    txt = []

                    # Append Latitude
                    try:
                        lat = tmp_fn(tck)
                        txt.append(('%'+ytick_lat_format+'$^o$') % lat)
                    except:
                        txt.append('')

                    # Append Range
                    txt.append('%d' % tck)
                    txt = '\n'.join(txt)

                    ytick_str.append(txt) # Put both lat and range on same string
                axis.set_yticklabels(ytick_str,rotation=90,ma='center') # Set yticklabels
                # Label y-axis
                geo_mag = 'Geographic' if currentData.fov.coords == 'geo' else 'Magnetic'
                if metadata['model'] == 'IS':
                    axis.set_ylabel('%s Latitude\nSlant Range [km]' % geo_mag,labelpad=y_labelpad)
                elif metadata['model'] == 'GS':
                    axis.set_ylabel('GS Mapped %s Latitude\nGS Mapped Range [km]' % geo_mag,labelpad=y_labelpad)
            else:
                if metadata['model'] == 'IS':
                    axis.set_ylabel('Slant Range [km]',labelpad=y_labelpad)
                elif metadata['model'] == 'GS':
                    axis.set_ylabel('GS Mapped Range [km]',labelpad=y_labelpad)

        axis.set_ylim(ylim)
        # Shade xBoundary Limits
        if xBoundaryLimits is None:
            if 'timeLimits' in currentData.metadata:
                xBoundaryLimits = currentData.metadata['timeLimits']

        if xBoundaryLimits is not None:
            gray = '0.75'
#            axis.axvspan(xlim[0],xBoundaryLimits[0],color=gray,zorder=150,alpha=0.5)
#            axis.axvspan(xBoundaryLimits[1],xlim[1],color=gray,zorder=150,alpha=0.5)
            axis.axvspan(xlim[0],xBoundaryLimits[0],color=gray,zorder=1)
            axis.axvspan(xBoundaryLimits[1],xlim[1],color=gray,zorder=1)
            axis.axvline(x=xBoundaryLimits[0],color='g',ls='--',lw=2,zorder=150)
            axis.axvline(x=xBoundaryLimits[1],color='g',ls='--',lw=2,zorder=150)

        # Shade yBoundary Limits
        if yBoundaryLimits is None:
            if 'gateLimits' in currentData.metadata and coords == 'gate':
                yBoundaryLimits = currentData.metadata['gateLimits']

            if 'rangeLimits' in currentData.metadata and coords == 'range':
                yBoundaryLimits = currentData.metadata['rangeLimits']

        if yBoundaryLimits is not None:
            gray = '0.75'
#            axis.axhspan(ylim[0],yBoundaryLimits[0],color=gray,zorder=150,alpha=0.5)
#            axis.axhspan(yBoundaryLimits[1],ylim[1],color=gray,zorder=150,alpha=0.5)
            axis.axhspan(ylim[0],yBoundaryLimits[0],color=gray,zorder=1)
            axis.axhspan(yBoundaryLimits[1],ylim[1],color=gray,zorder=1)
            axis.axhline(y=yBoundaryLimits[0],color='g',ls='--',lw=2,zorder=150)
            axis.axhline(y=yBoundaryLimits[1],color='g',ls='--',lw=2,zorder=150)
        
            for bnd_item in yBoundaryLimits:
                if coords == 'gate':
                    txt = []
                    txt.append('%d' % bnd_item)

                    rg_inx = np.where(bnd_item == currentData.fov.gates)[0]
                    if np.size(rg_inx) != 0:
                        lat = currentData.fov.latCenter[beamInx,rg_inx]
                        if np.isfinite(lat): 
                            txt.append('%.1f$^o$' % lat)
                        else:
                            txt.append('')
                    txt = '\n'.join(txt)
                else:
                    txt = '%.1f' % bnd_item
                if plot_range_limits_label:
                    axis.annotate(txt, (1.01, bnd_item) ,xycoords=('axes fraction','data'),rotation=90,ma='center')

        if plot_cbar:
            cbar = fig.colorbar(pcoll,orientation='vertical',shrink=cbar_shrink,fraction=cbar_fraction)
            cbar.set_label(cbarLabel)
            if cbar_ticks is None:
                labels = cbar.ax.get_yticklabels()
                labels[-1].set_visible(False)
            else:
                cbar.set_ticks(cbar_ticks)

            if 'gscat' in currentData.metadata:
                if currentData.metadata['gscat'] == 1:
                    cbar.ax.text(0.5,cbar_gstext_offset,'Ground\nscat\nonly',ha='center',fontsize=cbar_gstext_fontsize)

        txt = 'Model: ' + metadata['model']
        axis.text(1.01, 0, txt,
                horizontalalignment='left',
                verticalalignment='bottom',
                rotation='vertical',
                size=model_text_size,
                transform=axis.transAxes)

        # Get axis position information.
        pos = list(axis.get_position().bounds)

        # Plot frequency and noise information. ######################################## 
        if hasattr(dataObject,'prm') and plot_info:
            # Adjust current plot position to fit in the freq and noise plots.
            super_plot_hgt  = 0.06
            pos[3] = pos[3] - (2*super_plot_hgt)
            axis.set_position(pos)

            # Get current colorbar position and adjust it.
            cbar_pos = list(cbar.ax.get_position().bounds)
            cbar_pos[1] = pos[1]
            cbar_pos[3] = pos[3]
            cbar.ax.set_position(cbar_pos)

            curr_xlim   = axis.get_xlim()
            curr_xticks = axis.get_xticks()

            pos[1]      = pos[1] + pos[3]
            pos[3]      = super_plot_hgt
            freq_pos    = pos[:]

            pos[1]      = pos[1] + super_plot_hgt
            noise_pos   = pos[:]

            skynoise_ax = fig.add_axes(noise_pos, label='sky')
            searchnoise_ax = fig.add_axes(noise_pos, label='search', frameon=False)
            freq_ax = fig.add_axes(freq_pos, label='freq')
            nave_ax = fig.add_axes(freq_pos, label='nave', frameon=False)
#            cpid_ax = fig.add_axes(cpid_pos)
            plot_freq(freq_ax,dataObject.prm.time,dataObject.prm.tfreq,xlim=curr_xlim,xticks=curr_xticks)
            plot_nave(nave_ax,dataObject.prm.time,dataObject.prm.nave,xlim=curr_xlim,xticks=curr_xticks)

            plot_skynoise(skynoise_ax,dataObject.prm.time,dataObject.prm.noisesky,xlim=curr_xlim,xticks=curr_xticks)
            plot_searchnoise(searchnoise_ax,dataObject.prm.time,dataObject.prm.noisesearch,xlim=curr_xlim,xticks=curr_xticks)

        # Put a title on the RTI Plot. #################################################
        if plot_title:
            title_y = (pos[1] + pos[3]) + 0.015
            xmin    = pos[0]
            xmax    = pos[0] + pos[2]

            txt     = metadata['name']+'  ('+metadata['fType']+')'
            fig.text(xmin,title_y,txt,ha='left',weight=550)

            txt     = []
            txt.append(xlim[0].strftime('%Y %b %d %H%M UT - ')+xlim[1].strftime('%Y %b %d %H%M UT'))
            txt.append(currentData.history[max(currentData.history.keys())]) # Label the plot with the current level of data processing.
            txt     = '\n'.join(txt)
            fig.text((xmin+xmax)/2.,title_y,txt,weight=550,size='large',ha='center')

            txt     = 'Beam '+str(beam)
            fig.text(xmax,title_y,txt,weight=550,ha='right')

        cbar_info           = {}
        cbar_info['cmap']   = cmap
        cbar_info['bounds'] = bounds 
        cbar_info['norm']   = norm 
        cbar_info['label']  = cbarLabel
        cbar_info['ticks']  = cbar_ticks
        cbar_info['mappable']  = pcoll
        self.cbar_info      = cbar_info
