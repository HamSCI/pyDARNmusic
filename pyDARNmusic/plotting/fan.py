from matplotlib.collections import PolyCollection
import matplotlib

import numpy as np

from pyDARNmusic import getDataSet
from pyDARNmusic.utils.radUtils import getParamDict
from pyDARNmusic.utils.plotUtils import genCmap

#Global Figure Size
figsize=(20,18)

class musicFan(object):
    """Class to plot a fan plot using a pydarn.proc.music.musicArray object as the data source.

    Parameters
    ----------
    dataObj : pydarn.proc.music.musicArray
        musicArray object
    dataSet : Optional[str]
        Which dataSet in the musicArray object to plot
    time : Optional[None or datetime.datetime]
        Time scan plot.  If None, the first time in dataSet will be used.
    axis : Optional[None or matplotlib.figure.axis]
        Matplotlib axis on which to plot.  If None, a new figure and axis will be created.
    scale : Optional[None or 2-Element iterable]
        Colorbar scale.  If None, the default scale for the current SuperDARN parameter will be used.
    autoScale : Optional[bool]
        If True, automatically scale the color bar for good data visualization. Keyword scale must
        be None when using autoScale.
    plotZeros : Optional[bool]
        If True, plot cells that are exactly 0.
    markCell : Optional[None or 2-Element iterable]
        Mark the (beam, rangeGate) with black.
    markBeam : Optional[None or int]
        Mark a chosen beam.
    markBeam_dict : Optional[dict]
        dictionary of keywords defining markBeam line properties.
    plotTerminator : Optional[bool]
        If True, overlay day/night terminator on map.  Uses Basemap's nightshade.
    plot_title : Optional[bool]
        If True, plot the title information
    title : Optional[str]
        Overide default title text.
    parallels_ticks : Optional[list]
        Where to draw the parallel (latitude) lines
    meridians_ticks : Optional[list]
        Where to draw the meridian (longitude) lines
    zoom : Optional[float]
        Multiply the map height and width by this factor (bigger number shows more area).
    lat_shift : Optional[float]
        Add this number to the computed lat_0 sent to basemap.
    lon_shift : Optional[float]
        Add this number to the computed lon_0 sent to basemap.
    cmap_handling : Optional[str]
        'superdarn' to use SuperDARN-style colorbars, 'matplotlib' for direct use of matplotlib's colorbars.
        'matplotlib' is recommended when using custom scales and the 'superdarn' mode is not providing a desirable result.
    cmap : Optional[one or matplotlib colormap object]
        If Nonei and cmap_handling=='matplotlib', use jet.
    plot_cbar : Optional[bool]
        If True, plot the color bar.
    cbar_ticks : Optional[list]
        Where to put the ticks on the color bar.
    cbar_shrink : Optional[float]
        Fraction by which to shrink the colorbar
    cbar_fraction : Optional[float]
        Fraction of original axes to use for colorbar
    cbar_gstext_offset : Optional[float]
        y-offset from colorbar of "Ground Scatter Only" text
    cbar_gstext_fontsize : Optional[float]
        Fontsize of "Ground Scatter Only" text
    model_text_size : Optional[int]
        fontsize of model and coordinate indicator text
    draw_coastlines : Optional[bool]
        If True, draw the coastlines.
    basemap_dict : Optional[dict]
        Dictionary of keywords sent to the basemap invocation
    **kwArgs
        Keyword Arguments


    Attributes
    ----------
    map_obj

    pcoll

    Written by Nathaniel A. Frissell, Fall 2013
    Updated by: Francis Tholley, 2022
    """
    def __init__(self,dataObject,
        dataSet                 = 'active',
        time                    = None,
        fig                     = None,
        axis                    = None,
        subplot_tuple           = (1,1,1),
        scale                   = None,
        autoScale               = False,
        plotZeros               = False,
        cmap_handling           = 'superdarn',
        cmap                    = None,
        plot_cbar               = True,
        cbar_ticks              = None,
        cbar_shrink             = 1.0,
        cbar_fraction           = 0.15,
        cbar_gstext_offset      = -0.075,
        cbar_gstext_fontsize    = None,
        model_text_size         = 'small',
        plot_title              = True,
        title                   = None,
        **kwArgs):

        from matplotlib import pyplot as plt
        if fig is None:
            fig   = plt.figure(figsize=figsize)

        from scipy import nanstd, nanmean
        try:
            from cartopy.mpl import geoaxes
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature
            from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
            cartopyInstalled = True
        except Exception:
            cartopyInstalled = False
            
        from pydarn import (time2datetime, SuperDARNRadars,
                    Projs, Coords, Hemisphere, RangeEstimation)
        # Make some variables easier to get to...
        currentData = getDataSet(dataObject,dataSet)
        metadata    = currentData.metadata
        latFull     = currentData.fov["latFull"]
        lonFull     = currentData.fov["lonFull"]
        sdate       = currentData.time[0]
        coords      = metadata['coords']
        stid = metadata['stid']
        projs = Projs.GEO
        # Get center of FOV.
        #Determine center beam.
        ctrBeamInx  = len(currentData.fov["beams"])/2
        ctrGateInx  = len(currentData.fov["gates"])/2
        ctrLat      = currentData.fov["latCenter"][int(ctrBeamInx),int(ctrGateInx)]
        ctrLon      = currentData.fov["lonCenter"][int(ctrBeamInx),int(ctrGateInx)]
        # import ipdb;ipdb.set_trace()
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
                sd          = nanstd(np.abs(currentData.data),axis=None)
                mean        = nanmean(np.abs(currentData.data),axis=None)
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

        # Handles projection of data based on radar location
        deg_from_midnight = (sdate.hour + sdate.minute / 60) / 24 * 360
        hemisphere = SuperDARNRadars.radars[stid].hemisphere
        grid_lines = True
        if hemisphere == Hemisphere.North:
            pole_lat = 90
            noon = -deg_from_midnight
        else:
            pole_lat = -90
            noon = 360 - deg_from_midnight

        if stid:
            radar_lat = SuperDARNRadars.radars[stid].hardware_info.geographic.lat
            radar_lon = SuperDARNRadars.radars[stid].hardware_info.geographic.lon
        # handle none types or wrongly built axes
        # noon = noon + 12
        proj = ccrs.Orthographic(radar_lon,radar_lat)
        # proj = ccrs.SouthPolarStereo(noon,pole_lat)
        # import ipdb;ipdb.set_trace()
        if axis is None:
            axis = plt.subplot(*subplot_tuple, projection=proj, aspect='auto')

        else:
            fig   = axis.get_figure()

        extent = min(13e5,
                         (abs(proj.transform_point(noon, 30,
                                                   ccrs.PlateCarree())
                              [1])))

        axis.set_extent(extents=(-extent, extent, -extent, extent),
                          crs=proj)  
        if grid_lines:
                axis.gridlines(draw_labels=True,linewidth=1, color='black')
        axis.coastlines(resolution="110m")
        axis.add_feature(cfeature.LAND, color='lightgrey')
        axis.add_feature(cfeature.OCEAN, color = 'white')
        
        # Figure out which scan we are going to plot...
        if time is None:
            timeInx = 1
        else:
            timeInx = (np.where(currentData.time >= time))[0]
            # import ipdb;ipdb.set_trace()
            if np.size(timeInx) == 0:
                timeInx = -1
            else:
                timeInx = int(np.min(timeInx))
                

        # do some stuff in map projection coords to get necessary width and height of map
        lonFull,latFull = (np.array(lonFull)+360.)%360.,np.array(latFull)

        goodLatLon  = np.logical_and( np.logical_not(np.isnan(lonFull)), np.logical_not(np.isnan(latFull)) )
        goodInx     = np.where(goodLatLon)
        goodLatFull = latFull[goodInx]
        goodLonFull = lonFull[goodInx]

        # Plot the SuperDARN data!
        ngates = np.shape(currentData.data)[2]
        nbeams = np.shape(currentData.data)[1]
        data  = currentData.data[timeInx,:,:]
        # import ipdb;ipdb.set_trace()
        verts = []
        scan  = []
        # data  = currentData.data[timeInx,:,:]
        goodBmRg=[]
        geo = ccrs.Geodetic()
        for bm in range(nbeams):
            for rg in range(ngates):
                if goodLatLon[bm,rg] == False: continue
                if np.isnan(data[bm,rg]): continue
                if data[bm,rg] == 0 and not plotZeros: continue
                goodBmRg.append((bm,rg))
                scan.append(data[bm,rg])
                x1,y1 = proj.transform_point(lonFull[bm+0,rg+0],latFull[bm+0,rg+0],geo)
                x2,y2 = proj.transform_point(lonFull[bm+1,rg+0],latFull[bm+1,rg+0],geo)
                x3,y3 = proj.transform_point(lonFull[bm+1,rg+1],latFull[bm+1,rg+1],geo)
                x4,y4 = proj.transform_point(lonFull[bm+0,rg+1],latFull[bm+0,rg+1],geo)
                verts.append(((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x1,y1)))
        
        # np.savetxt('lat.csv', good_lat, delimiter=',')
        # np.savetxt('data.csv', good_data, delimiter=',')

        if (cmap_handling == 'matplotlib') or autoScale:
            if cmap is None:
                cmap    = matplotlib.cm.jet
            bounds  = np.linspace(scale[0],scale[1],256)
            norm    = matplotlib.colors.BoundaryNorm(bounds,cmap.N)
        elif cmap_handling == 'superdarn':
            colors  = 'lasse'
            cmap,norm,bounds = genCmap(param,scale,colors=colors)

#        pcoll = PolyCollection(np.array(verts),edgecolors='face',linewidths=0,closed=False,cmap=cmap,norm=norm,zorder=99)
        pcoll = PolyCollection(np.array(verts),edgecolors='face',closed=False,cmap=cmap,norm=norm,zorder=99)
        pcoll.set_array(np.array(scan))
        axis.add_collection(pcoll,autolim=False)
        radar_lat = []
        radar_lon = []

        # plt.scatter(radar_lon, radar_lat, c="black", s=20, transform=geo)
        
        dataName = currentData.history[max(currentData.history.keys())] # Label the plot with the current level of data processing.
        if plot_title:
            if title is None:
                axis.set_title(metadata['name']+' - '+dataName+currentData.time[timeInx].strftime('\n%Y %b %d %H%M UT')) 
            else:
                axis.set_title(title)

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
                    cbar.ax.text(0.5,cbar_gstext_offset,'Ground\nscat\nonly',ha='center',fontsize=cbar_gstext_fontsize,transform=cbar.ax.transAxes)

        txt = 'Coordinates: ' + metadata['coords'] +', Model: ' + metadata['model']
        axis.text(1.01, 0, txt,
                  horizontalalignment='left',
                  verticalalignment='bottom',
                  rotation='vertical',
                  size=model_text_size,
                  weight='bold',
                  transform=axis.transAxes)

        self.pcoll      = pcoll
