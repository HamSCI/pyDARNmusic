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

"""Range-time-intensity plotting

A module for generating rti plots.

Module author: AJ, 20130123

Functions
--------------------------------------------------
plot_rti            range-time-intensity plot
plot_freq           TX frequency data
plot_searchnoise    noise panel
plot_skynoise       sky noise panel
plot_cpid           control program ID panel
plot_nave           number of averges panel
rti_title           title an rti plot
draw_axes           draw empty axes
read_data           read data in
rti_panel           plot the main rti data
daynight_terminator calculate day/night terminator
--------------------------------------------------

"""
import logging

def plot_skynoise(ax, times, sky, xlim=None, xticks=None):
    """Plots a noise panel at position pos.

    Parameters
    ----------
    ax :
        a MPL axis object to plot to
    times : list
        a list of the times of the beam soundings
    sky: list
        a list of the noise.sky of the beam soundings
    search : list
        a list of the noise.search param
    xlim : Optional[list]
        2-element limits of the x-axis.  None for default.
    xticks : Optional[list]
        List of xtick poisitions.  None for default.

    Returns
    -------
    Nothing

    Example
    -------
        plot_skynoise(ax,times,sky)

    Written by AJ 20121002
    Modified by NAF 20131101
    Modified by ASR 20150916

    """

    from matplotlib.ticker import MultipleLocator
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    import numpy as np

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction='out')
    ax.set_ylim(bottom=0, top=6)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction='out', which='minor')

    # Plot the sky noise data.
    ax.plot_date(date2num(times), np.log10(sky), fmt='k-',
                 tz=None, xdate=True, ydate=False)

    # Format the xaxis.
    if xlim is not None: ax.set_xlim(xlim)
    if xticks is not None: ax.set_xticks(xticks)

    # Add labels to identify the noise axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(pos[0] - .01, pos[1] + .004, '10^0', ha='right', va='bottom',
             size=8)
    fig.text(pos[0] - .01, pos[1] + pos[3], '10^6', ha='right', va='top',
             size=8)
    fig.text(pos[0] - .07, pos[1] + pos[3] / 2., 'N.Sky', ha='center',
             va='center', size=8.5, rotation='vertical')
    l = Line2D([pos[0] - .06, pos[0] - .06], [pos[1] + .01,
               pos[1] + pos[3] - .01], transform=fig.transFigure,
               clip_on=False, ls='-', color='k', lw=1.5)
    ax.add_line(l)
    ax.set_xticklabels([' '])
    # Only use 2 major yticks.
    ax.set_yticks([0, 6])
    ax.set_yticklabels([' ', ' '])


def plot_searchnoise(ax, times, search, xlim=None, xticks=None,
                     ytickside='right'):
    """Plots a noise panel at position pos.

    Parameters
    ----------
    ax :
        a MPL axis object to plot to
    times : list
        a list of the times of the beam soundings
    sky : list
        a list of the noise.sky of the beam soundings
    search : list
        a list of the noise.search param
    xlim : Optional[list]
        2-element limits of the x-axis.  None for default.
    xticks : Optional[list]
        List of xtick poisitions.  None for default.
    ytickside : Optional[string]
        Default is right.

    Returns
    -------
    Nothing

    Example
    -------
        plot_searchnoise(ax,times,search)

    Written by AJ 20121002
    Modified by NAF 20131101
    Modified by ASR 20150916

    """

    from matplotlib.ticker import MultipleLocator
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D
    import numpy as np

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction='out')
    ax.set_ylim(bottom=0, top=6)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction='out', which='minor')

    # Plot the search noise data.
    ax.plot_date(date2num(times), np.log10(search),
                 fmt='k:', tz=None, xdate=True, ydate=False, lw=1.5)

    # Format the xaxis.
    if xlim is not None: ax.set_xlim(xlim)
    if xticks is not None: ax.set_xticks(xticks)

    # Add labels to identify the noise axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]

    fig.text(pos[0] + pos[2] + .01, pos[1] + .004, '10^0', ha='left',
             va='bottom', size=8)
    fig.text(pos[0] + pos[2] + .01, pos[1] + pos[3], '10^6', ha='left',
             va='top', size=8)
    fig.text(pos[0] + pos[2] + .06, pos[1] + pos[3] / 2., 'N.Sch', ha='center',
             va='center', size=8.5, rotation='vertical')

    l = Line2D([pos[0] + pos[2] + .07, pos[0] + pos[2] + .07],
               [pos[1] + .01, pos[1] + pos[3] - .01],
               transform=fig.transFigure, clip_on=False, ls=':',
               color='k', lw=1.5)
    ax.add_line(l)
    ax.set_xticklabels([' '])
    # use only 2 major yticks
    ax.set_yticks([0, 6])
    ax.set_yticklabels([' ', ' '])
    if ytickside == 'right':
        ax.yaxis.tick_right()


def plot_freq(ax, times, freq, xlim=None, xticks=None):
    """Plots the tx frequency data to an axis object.

    Parameters
    ----------
    ax :
        a MPL axis object to plot to
    times : list
        a list of the times of the beam soundings
    freq : list
        a list of the tfreq of the beam soundings
    xlim : Optional[list]
        2-element limits of the x-axis.  None for default.
    xticks : Optional[list]
        List of xtick poisitions.  None for default.

    Returns
    -------
    Nothing.


    Example
    -------
        plot_freq(ax, times, tfreq)

    Written by AJ 20121002
    Modified by NAF 20131101
    Modified by ASR 20150916

    """

    from matplotlib.ticker import MultipleLocator
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D

    # Format the yaxis.
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction='out')
    ax.set_ylim(bottom=8, top=20)
    ax.yaxis.set_minor_locator(MultipleLocator())
    ax.yaxis.set_tick_params(direction='out', which='minor')

    # Plot the TX frequency.
    ax.plot_date(date2num(times), freq, fmt='k-',
                 tz=None, xdate=True, ydate=False, markersize=2)

    # Format the xaxis.
    if xlim is not None: ax.set_xlim(xlim)
    if xticks is not None: ax.set_xticks(xticks)

    # Add labels to identify the frequency axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(pos[0] - .01, pos[1] + .005, '10', ha='right', va='bottom',
             size=8)
    fig.text(pos[0] - .01, pos[1] + pos[3] - .015, '16', ha='right', va='top',
             size=8)
    fig.text(pos[0] - .07, pos[1] + pos[3] / 2., 'Freq', ha='center',
             va='center', size=9, rotation='vertical')
    fig.text(pos[0] - .05, pos[1] + pos[3] / 2., '[MHz]', ha='center',
             va='center', size=7, rotation='vertical')
    l = Line2D([pos[0] - .04, pos[0] - .04], [pos[1] + .01,
               pos[1] + pos[3] - .01], transform=fig.transFigure,
               clip_on=False, ls='-', color='k', lw=1.5)
    ax.add_line(l)
    ax.set_xticklabels([' '])
    # use only 2 major yticks
    ax.set_yticks([10, 16])
    ax.set_yticklabels([' ', ' '])


def plot_nave(ax, times, nave, xlim=None, xticks=None, ytickside='right'):
    """Plots the number of averages (nave) data to an axis object.

    Parameters
    ----------
    ax :
        a MPL axis object to plot to
    times : list
        a list of the times of the beam soundings
    nave : list
        a list of the nave of the beam soundings
    xlim : Optional[list]
        2-element limits of the x-axis.  None for default.
    xticks : Optional[list]
        List of xtick poisitions.  None for default.
    ytickside : Optional[str]
        Default is right.

    Returns
    -------
    Nothing.

    Example
    -------
        plot_nave(ax, times, nave)

    Written by AJ 20121002
    Modified by NAF 20131101
    Modified by ASR 20150916

    """

    from matplotlib.ticker import MultipleLocator
    from matplotlib.dates import date2num
    from matplotlib.lines import Line2D

    # Format the yaxis
    ax.yaxis.tick_left()
    ax.yaxis.set_tick_params(direction='out')
    ax.set_ylim(bottom=0, top=80)
    ax.yaxis.set_minor_locator(MultipleLocator(base=5))
    ax.yaxis.set_tick_params(direction='out', which='minor')

    # Plot the number of averages.
    ax.plot_date(date2num(times), nave, fmt='k:',
                 tz=None, xdate=True, ydate=False, markersize=2)

    # Format the xaxis.
    if xlim is not None: ax.set_xlim(xlim)
    if xticks is not None: ax.set_xticks(xticks)

    # Add labels to identify the nave axis.
    fig = ax.get_figure()
    bb = ax.get_position()
    x0 = bb.x0
    y0 = bb.y0
    height = bb.height
    width = bb.width
    pos = [x0, y0, width, height]
    fig.text(pos[0] + pos[2] + .01, pos[1] - .004, '0', ha='left', va='bottom',
             size=8)
    fig.text(pos[0] + pos[2] + .01, pos[1] + pos[3], '80', ha='left', va='top',
             size=8)
    fig.text(pos[0] + pos[2] + .06, pos[1] + pos[3] / 2., 'Nave', ha='center',
             va='center', size=8.5, rotation='vertical')

    l = Line2D([pos[0] + pos[2] + .07, pos[0] + pos[2] + .07],
               [pos[1] + .01, pos[1] + pos[3] - .01],
               transform=fig.transFigure, clip_on=False, ls=':',
               color='k', lw=1.5)
    ax.add_line(l)
    ax.set_xticklabels([' '])
    # use only 2 major yticks
    ax.set_yticks([0, 80])
    ax.set_yticklabels([' ', ' '])
    if ytickside == 'right':
        ax.yaxis.tick_right()

