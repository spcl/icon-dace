#
# -*- coding: utf-8 -*-
#
# util_plot.py - Utilities for handling and plotting 2D data and maps
# 
# Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#

"""
###############################################################################
## Utility for handling and plotting 2D data and maps
#
###############################################################################
## Author: Stefan Hagemann (Hereon, Institute of Coastal Systems, GER)
##
"""

import os
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

import datetime as dt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes

from psy_simple.plotters import CMap, Bounds, MissColor
from psy_maps.plotters import Transform, MapPlot2D
from psyplot.plotter import Plotter

GeoAxes._pcolormesh_patched = Axes.pcolormesh


# ***************************************************
def plot_map(nlon, nlat, fdata, var, reg, lmask,
             flat_n, flat_s, flon_w, flon_e,
             opt, xmiss,
             dnout='map',
             plotpath='./',
             ext_plot='.pdf'):
    # ***************************************************
    #
    # Make 2D plot analogous to wreg.
    #
    # Basics
    now = dt.datetime.now()
    pdate = "Stefan Hagemann, Hereon, " + now.strftime("%Y-%m-%d")
    print(pdate)

    colmap = plt.get_cmap(opt.cmap, opt.ncol)
    #  colmap.set_over('gray')
    cproj = ccrs.PlateCarree()
    #
    # Deal with missing values
    #  mapdata = fdata.astype(float)
    mapdata = mask_missval(fdata, xmiss, opt.xlow)

    fig, ax = plt.figure(figsize=(10, 6)), {}
    if opt.ibarloc == 0:
        fig.subplots_adjust(left=0.05, right=0.95)
    elif opt.ibarloc == 1:
        fig.subplots_adjust(left=0.06, right=0.83)
    row, col = 1, 1
    im = 1
    ax = fig.add_subplot(row, col, im, projection=ccrs.PlateCarree())
    #
    # Define grid coordinates
    xlon, xlat = np.meshgrid(np.linspace(flon_w, flon_e, nlon), np.linspace(flat_n, flat_s, nlat))
    print(' arg: ', len(xlon), len(xlat), opt.xmin, opt.xmax, nlon, nlat)

    pmap = ax.pcolormesh(xlon, xlat, mapdata, transform=cproj, vmin=opt.xmin, vmax=opt.xmax, cmap=colmap,
                         rasterized=True, edgecolors='none')

    #  mapplot = ax.pcolormesh(lons, lats, mapdata, cmap=colmaps[var.split('_')[0]]['cmap'],
    #            norm=colmaps[var.split('_')[0]]['norm'], edgecolors='face', transform=ccrs.PlateCarree(), rasterized=True, zorder=0)

    print(nlon, nlat, cproj)
#
#   *** Set plot title, axis range & legends and background map settings
    set_plot_char(ax, cproj, opt, reg, var)

    if opt.imask == 1:
        draw_mask_borders_map(ax, lmask, flon_w, flon_e, flat_s, flat_n)

    #  cbar = mpl.colorbar.ColorbarBase(caxe, cmap=colmaps[var]['cmap'], norm=colmaps[var]['norm'],
    #      orientation='horizontal', ticks=colmaps[var]['bounds'], extend='both')
    #  cbar.set_label('Grid cell fraction [/] and cell average pond reservoir depth [m]')
    #
    # ******* Colorbar ********************
    valmin = mapdata.min()
    valmax = mapdata.max()
    print('Min: ', valmin, ' Max: ', valmax)
    cdum_ext = 'neither'
    if valmax > opt.xmax:
        cdum_ext = 'max'
    if valmin < opt.xmin and opt.xlow < opt.xmin:
        cdum_ext = 'min'
        if valmax > opt.xmax:
            cdum_ext = 'both'

    pmap.set_clim([opt.xmin, opt.xmax])  # set color limits
    if opt.ibarloc == 0:
        if opt.iclass == 0:
            barticks = get_barticks(opt.xmin, opt.xmax, opt.ncol)
            fig.colorbar(pmap, orientation="horizontal", fraction=0.07, pad=0.15, extend=cdum_ext, ticks=barticks,
                         label=var.unit)
        elif opt.iclass == 1:
            cb = fig.colorbar(pmap, orientation="horizontal", fraction=0.07, pad=0.15, extend=cdum_ext)
            labels = np.arange(0, opt.ncol, 1) + 1
            loc = labels - .5
            cb.set_ticks(loc)
            cb.set_ticklabels(labels)
        else:
            print('iclass > 1 for horizontal plot not implemented --> Abbruch')
            exit()
        if opt.pdate == 1:
            plt.figtext(0.95, 0.5, pdate, horizontalalignment='right', verticalalignment='center', rotation='vertical',
                        fontsize=10)

    elif opt.ibarloc == 1:
        if opt.iclass == 0:
            barticks = get_barticks(opt.xmin, opt.xmax, opt.ncol)
            cb = fig.colorbar(pmap, orientation="vertical", fraction=0.07, pad=0.05, extend=cdum_ext, ticks=barticks,
                              label=var.unit)
        else:

            # Now adding the colorbar: [left, bottom, width, height], fractions of plotting area
            cbaxes = fig.add_axes([0.82, 0.09, 0.03, 0.759])
            cb = plt.colorbar(pmap, cax=cbaxes)

            #      cb.ax.tick_params(labelsize=16)   good for short names
            cb.ax.tick_params(labelsize=11)
            #      cb=fig.colorbar(pmap,orientation="vertical", fraction=0.07, pad=0.02, extend=cdum_ext)
            labels = np.arange(0, opt.ncol, 1) + 1
            loc = labels - .5
            cb.set_ticks(loc)
            if opt.iclasstyp == 0:
                cb.set_ticklabels(labels)
            else:
                cb.set_ticklabels(opt.barlabel)
        if opt.pdate == 1:
            plt.figtext(0.95, 0.01, pdate, horizontalalignment='right', verticalalignment='center',
                        rotation='horizontal', fontsize=10)

    # adjustment of the text downwards:
    #  tdown = 0.05
    #  plt.figtext(0.9, -tdown, pdate, horizontalalignment='right', fontsize=10)

    #  plt.show()
    plotfile = plotpath + os.sep + dnout + ext_plot

    fig.savefig(plotfile, dpi=200, format='pdf')
    #  fig.savefig(plotfile, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=200)
    plt.close()
    print('Plot --> ' + plotfile)


# ***************************************************
def plot_rivdir(nlon, nlat, fdata, reg, lmask,
                flat_n, flat_s, flon_w, flon_e,
                opt, dnout='rdf', plotpath='./',
                ext_plot='.pdf'):
    # ***************************************************
    #
    # Make 2D plot analogous to wreg.
    #
    # Basics
    now = dt.datetime.now()
    pdate = "Stefan Hagemann, Hereon, " + now.strftime("%Y-%m-%d")
    print(pdate)

    reg.xticks = get_geoticks(reg.lonmin, reg.lonmax)
    reg.yticks = get_geoticks(reg.latmin, reg.latmax)

    cproj = ccrs.PlateCarree()

    deltalon = reg.lonmax - reg.lonmin
    deltalat = reg.latmax - reg.latmin
    if deltalon > deltalat:
        fig, ax = plt.figure(figsize=(12, 6)), {}
    elif deltalon < 0.75 * deltalat:
        fig, ax = plt.figure(figsize=(8, 12)), {}
    else:
        fig, ax = plt.figure(figsize=(10, 6)), {}

    fig.subplots_adjust(left=0.05, right=0.95)
    row, col = 1, 1
    im = 1
    ax = fig.add_subplot(row, col, im, projection=ccrs.PlateCarree())
    #
    # Define grid coordinates
    dx0 = (flon_e - flon_w) / nlon
    dy0 = (flat_n - flat_s) / nlat
    xlon, xlat = np.meshgrid(np.linspace(flon_w + dx0 / 2., flon_e - dx0 / 2., nlon),
                             np.linspace(flat_n - dy0 / 2., flat_s + dy0 / 2., nlat))

    print('Grid Origin: ', xlon.min(), xlat.max())
    print('Grid SE: ', xlon.max(), xlat.min())

    jbmin = nlat - 1
    jbmax = 0
    jlmin = nlon - 1
    jlmax = 0
    #
    # ************ Plot direction arrows
    if opt.idir == 1:
        #
        # Define arrows
        x = []
        y = []
        dx = []
        dy = []
        print(len(fdata), ' x', len(fdata[0]), len(xlon), ' x ', len(xlon[0]), dx0, dy0)

        for jb in range(nlat):
            for jl in range(nlon):
                if reg.lonmin <= xlon[jb][jl] <= reg.lonmax and reg.latmin <= xlat[jb][jl] <= reg.latmax:
                    jbmin = min(jbmin, jb)
                    jbmax = max(jbmax, jb)
                    jlmin = min(jlmin, jl)
                    jlmax = max(jlmax, jl)

                    if fdata[jb][jl] == 1 or fdata[jb][jl] == 4 or fdata[jb][jl] == 7:
                        x.append(xlon[jb][jl])
                        dx.append(-dx0)
                    elif fdata[jb][jl] == 2 or fdata[jb][jl] == 8:
                        x.append(xlon[jb][jl])
                        dx.append(0.)
                    elif fdata[jb][jl] == 3 or fdata[jb][jl] == 6 or fdata[jb][jl] == 9:
                        x.append(xlon[jb][jl])
                        dx.append(dx0)
                    if fdata[jb][jl] == 1 or fdata[jb][jl] == 2 or fdata[jb][jl] == 3:
                        y.append(xlat[jb][jl])
                        dy.append(-dy0)
                    elif fdata[jb][jl] == 4 or fdata[jb][jl] == 6:
                        y.append(xlat[jb][jl])
                        dy.append(0.)
                    elif fdata[jb][jl] == 7 or fdata[jb][jl] == 8 or fdata[jb][jl] == 9:
                        y.append(xlat[jb][jl])
                        dy.append(dy0)

        print('Land Region Origin: ', min(x), max(y), ', Index origin: ', jlmin, jbmin)
        print('    Land Region SE: ', max(x), min(y), ' Index SE: ', jlmax, jbmax)
        #
        plt.quiver(x, y, dx, dy, units='x', headlength=2, headaxislength=2, pivot='middle', color='blue',
               scale_units='xy', scale=1.1)
        print(nlon, nlat, cproj)
    elif opt.idir == 2 or opt.idir == 3: 
        for jb in range(nlat):
            for jl in range(nlon):
                if reg.lonmin <= xlon[jb][jl] <= reg.lonmax and reg.latmin <= xlat[jb][jl] <= reg.latmax:
                    jbmin = min(jbmin, jb)
                    jbmax = max(jbmax, jb)
                    jlmin = min(jlmin, jl)
                    jlmax = max(jlmax, jl)
                    ax.text(xlon[jb][jl], xlat[jb][jl], '{:d}'.format(int(fdata[jb][jl])), fontsize='small',
                          ha='center', va='center' )
    #
    # *** set Title
    ax.set_title(opt.ptitle, fontsize=18)
    #
    # *** High resolution coastline
    if opt.icoast == 1:
        ax.coastlines(linewidth=0.5, resolution='10m')
    #
    # add feature: rivers and lakes
    if opt.iriver == 1:
        rivers = cf.NaturalEarthFeature(
            category='physical',
            name='rivers_lake_centerlines',
            scale='10m',
            facecolor='none')
        ax.add_feature(rivers, edgecolor='blue', linewidth=0.5)
        lakes = cf.NaturalEarthFeature(
            category='physical',
            name='lakes',
            scale='10m',
            facecolor='none')
        ax.add_feature(lakes, edgecolor='blue', alpha=0.5, linewidth=0.5)
    #
    # add feature: Country Borders
    if opt.icountry == 1:
        country = cf.NaturalEarthFeature(
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='10m',
            facecolor='none')
        ax.add_feature(country, edgecolor='red', linewidth=0.2)

    ax.set_extent([reg.lonmin, reg.lonmax, reg.latmin, reg.latmax], crs=cproj)
    ax.set_xlabel('Lon', fontsize=16)
    ax.set_ylabel('Lat', fontsize=16)

    #
    # *** Tick business
    ax.set_xticks(reg.xticks, minor=False)
    ax.set_yticks(reg.yticks, minor=False)
    ksize = 12
    if reg.lonmin <= -100. or reg.lonmax >= 100:
        ksize = 10
    ax.tick_params(axis='x', which='major', length=5, labelsize=ksize)
    ax.tick_params(axis='y', which='major', labelsize=12)
    #  ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    #  ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    #  ax.minorticks_on       # no effect

    #
    # *** Plot mask boundaries
    if opt.imask == 1:
        regmask = lmask[jbmin:jbmax + 1, jlmin:jlmax + 1]  # Letzte Spalte/Index does not count
        nmas = regmask.sum()
        print('True mask values within regional zoom: ', nmas)
        if nmas == 0:
            print(' Zero Mask --> ABBRUCH!')
            exit()
        draw_mask_borders(ax, regmask, reg.lonmin, reg.lonmax, reg.latmax, reg.latmin)
    #
    # *** Plot locations, e.g. stations
    #  xstat=[] ; ystat=[]
    #  xstat.append(141.1061) ; ystat.append(-27.5967)
    #  xstat.append(138.8475) ; ystat.append(-25.5151)
    if opt.ispot == 1:
        for i in range(len(opt.xspot)):
            print(i + 1, '. Spot:', opt.xspot[i], opt.yspot[i])
        draw_circles(ax, opt.xspot, opt.yspot)
    #
    # *** Plot Shapefile
    #
    # add map from shape file
    if opt.ishp == 1:
        dn_shape = opt.dir_shp + os.sep + opt.dn_shp  # e.g. "/h/hagemans/try/seth/shp/Paajako"
        draw_shapefile(ax, dn_shape, opt.id_epsg)
    #    draw_shapefile(ax, dn_shape, id_epsg=3047)

    #
    # *** Plot grid lines
    if opt.igitt == 1:
        ax.grid(color='gray', alpha=0.5, linestyle='solid')
    #  cbar = mpl.colorbar.ColorbarBase(caxe, cmap=colmaps[var]['cmap'], norm=colmaps[var]['norm'],
    #      orientation='horizontal', ticks=colmaps[var]['bounds'], extend='both')
    #  cbar.set_label('Grid cell fraction [/] and cell average pond reservoir depth [m]')
    #

    # adjustment of the text downwards:
    if opt.pdate == 1:
        plt.figtext(0.95, 0.5, pdate, horizontalalignment='right', verticalalignment='center', rotation='vertical',
                    fontsize=10)

    #  plt.show()
    plotfile = plotpath + os.sep + dnout + ext_plot

    fig.savefig(plotfile, dpi=200, format='pdf')
    #  fig.savefig(plotfile, bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=200)
    plt.close()
    print('Plot --> ' + plotfile)


#
# ***************************************************
def set_plot_char(ax, cproj, opt, reg, var):
    #
    # ******* Plot title, Axis legends and background map features ********

    #    globmean = '('+str(round(mapdata.mean() * 100,2)) + '%)'
    ax.set_title(opt.ptitle + ': ' + var.longname, fontsize=18)
    #
    # *** High resolution coastline
    if opt.icoast == 1:
        ax.coastlines(linewidth=0.5, resolution='10m')
    #
    # add feature: rivers and lakes
    if opt.iriver == 1:
        rivers = cf.NaturalEarthFeature(
            category='physical',
            name='rivers_lake_centerlines',
            scale='10m',
            facecolor='none')
        ax.add_feature(rivers, edgecolor='blue', linewidth=0.5)
        lakes = cf.NaturalEarthFeature(
            category='physical',
            name='lakes',
            scale='10m',
            facecolor='none')
        ax.add_feature(lakes, edgecolor='blue', alpha=0.5, linewidth=0.5)
    #
    # add feature: Country Borders
    if opt.icountry == 1:
        country = cf.NaturalEarthFeature(
            category='cultural',
            name='admin_0_boundary_lines_land',
            scale='10m',
            facecolor='none')
        ax.add_feature(country, edgecolor='red', linewidth=0.2)
    #
    # *** set plot extent, projection and axis labels
    ax.set_extent([reg.lonmin, reg.lonmax, reg.latmin, reg.latmax], crs=cproj)

    ax.set_xlabel('Lon', fontsize=16)
    ax.set_ylabel('Lat', fontsize=16)
    #
    # *** Tick business
    reg.xticks = get_geoticks(reg.lonmin, reg.lonmax)
    reg.yticks = get_geoticks(reg.latmin, reg.latmax)

    ax.set_xticks(reg.xticks, minor=False)
    ax.set_yticks(reg.yticks, minor=False)
    ax.tick_params(axis='x', which='major', length=5, labelsize=12)
    ax.tick_params(axis='y', which='major', labelsize=12)
    #  ax.minorticks_on       # no effect

    #
    # *** Plot grid lines
    if opt.igitt == 1:
        ax.grid(color='gray', alpha=0.5, linestyle='solid')

#
# ***************************************************
def mask_missval(fdata, xmiss, xlow):
    # ***************************************************
    #
    # Converts Missing value from value xmiss to None in array fdata of size nlon * nlat.
    # It also sets all values below xlow (xlow is still plotted) to None if xlow < 1.E10
    #
    zeps = 1.E-10
    zmax = 1.E10

    masked_array = np.ma.masked_where(fdata <= xmiss + zeps, fdata)
    if xlow < zmax:
        masked_array = np.ma.masked_where(fdata <= xlow, masked_array)

    return masked_array


# ***************************************************
def get_colortype(icoltyp, linter=False):
    # ***************************************************
    #
    if linter:
        print('Choose Your Colormap No.')
        print("  0  ==> Jet           (b-t-g-y-o-r)")
        print("  1  ==> Jet reversed")
        print("  2  ==> Terrain       (b-g-y-br-gr-w)")
        print("  3  ==> Nipy Spectral (bl-v-b-g-y-0-r-gr)")
        print("  4  ==> NCAR (gist)   (bl-kh-b-t-g-y-o-r-p-v-te-w)")
        print("  5  ==> Red to White to Blue    - Good for differences")
        print("  6  ==> Blue to White to Red    - Good for Temp. differences")
        print("  7  ==> Gray scale")
        print("  8  ==> Yellow to Green         - good for b/w")
        print("  9  ==> Red to White to Black   - good? for b/w differences")
        print(" 10  ==> Earth (gist)   (bl-b-g-br-te-w)")
        print(" 11  ==> Ocean          (g-bl-b-t-w)")
        print(" 12  ==> Spectral       (r-y w-g-b)")
        print(" 13  ==> Rev. Spectral  (b-g-w-y-r)")
        print(" 14  ==> Prism  Peridoic (r-b-g-r ...)")
        print(" 15  ==> Pastel2")
        print(" 16  ==> Paired - 12 classes")
        print(" 17  ==> Tab20 - 20 classes")
        #       print, " 20  ==> Set WATCH specific colors: blue - red"
        #       print, " 21  ==> Setting specific separating colors"
        #       print, " 22  ==> Set reversed WATCH specific colors red - blue"
        #       print, " 25  ==> Manual setting of colors"

        #    idum = raw_input("Enter colormap number: ")   # Python 2.x
        idum = input("Enter colormap number: ")  # Python 3
        icoltyp = int(idum)

    if icoltyp == 0:
        cmap = 'jet'
    elif icoltyp == 1:
        cmap = 'jet_r'
    elif icoltyp == 2:
        cmap = 'terrain'
    elif icoltyp == 3:
        cmap = 'nipy_spectral'
    elif icoltyp == 4:
        cmap = 'gist_ncar'
    elif icoltyp == 5:
        cmap = 'RdBu'
    elif icoltyp == 6:
        cmap = 'RdBu_r'
    elif icoltyp == 7:
        cmap = 'gray'
    elif icoltyp == 8:
        cmap = 'YlGn'
    elif icoltyp == 9:
        cmap = 'RdGy'
    elif icoltyp == 10:
        cmap = 'gist_earth'
    elif icoltyp == 11:
        cmap = 'ocean'
    elif icoltyp == 12:
        cmap = 'Spectral'
    elif icoltyp == 13:
        cmap = 'Spectral_r'
    elif icoltyp == 14:
        cmap = 'prism'
    elif icoltyp == 15:
        cmap = 'Pastel2'
    elif icoltyp == 16:
        cmap = 'Paired'
    elif icoltyp == 17:
        cmap = 'tab20'
    else:
        print('ERROR in get_colortype! Type ', icoltyp, ' does not exist!')
        cmap = None

    #   cmaps: RdBu_r, rainbow (v-b-g-r), jet, terrain

    return cmap


#
# ***************************************************
def get_barticks(xmin, xmax, ncol):
    # ***************************************************
    #
    # Depending on minimum and maximum values,
    # a vector/list of ticks is generated that can be used as colorbar description

    xrange = xmax - xmin
    dcol = xrange / ncol

    nticks = ncol
    while nticks > 14:
        nticks = int(nticks / 2)
    dtick = dcol * int(ncol / nticks)

    nticks += 1
    barticks = [0] * nticks
    for i in range(nticks):
        barticks[i] = xmin + i * dtick

    return barticks


#
# ***************************************************
def get_geoticks(degmin, degmax):
    # ***************************************************
    #
    # Depending on minimum and maximum longitude (or latitude),
    # a vector/list of ticks is generated that can be used as axis description
    # to plot a map over the respective region
    #
    zeps = 0.00001

    drange = degmax - degmin
    if drange > 300:
        dtick = 60
    elif 210 < drange <= 300:
        dtick = 40
    elif 140 < drange <= 210:
        dtick = 30
    elif 80 < drange <= 140:
        dtick = 20
    elif 40 < drange <= 80:
        dtick = 10
    elif 10 < drange <= 40:
        dtick = 5
    elif 2 < drange <= 10:
        dtick = 1
    else:
        dtick = 0.1

    drest0 = degmin % dtick
    if drest0 < zeps or abs(drest0 - dtick) < zeps:
        drest0 = dtick
    dspan = drange - (dtick - drest0)

    nticks = int(dspan / dtick) + 1
    tick0 = degmin - drest0 + dtick

    degticks = [0] * nticks
    for i in range(nticks):
        degticks[i] = tick0 + i * dtick

    return degticks


class Region:
    number = ""

    def __init__(self, number):
        self.no = number
        self.CNAME = "Region-Name"
        self.CSHORT = "Short_Name"
        self.latmin = -90.
        self.latmax = 90.
        self.lonmin = -180.
        self.lonmax = 180.
        self.xticks = [-180, -120, -60, 0, 60, 120, 180]
        self.yticks = [-90, -60, -30, 0, 30, 60, 90]

    def get_region(self):

        # routine returns status: 0 = ok, -1 = region no. not found

        status = 0
        if self.no == 0:
            self.CNAME = "Globe"
            self.CSHORT = "globe"
            self.latmin = -90.
            self.latmax = 90.
            self.lonmin = -180.
            self.lonmax = 180.
        elif self.no == 1:
            self.CNAME = "South_America"
            self.CSHORT = "sa"
            self.latmin = -60.
            self.latmax = 20.
            self.lonmin = -90.
            self.lonmax = -30.
        elif self.no == 2:
            self.CNAME = "North_America"
            self.CSHORT = "na"
            self.latmin = 20.
            self.latmax = 90.
            self.lonmin = -170.
            self.lonmax = -50.
        elif self.no == 3:
            self.CNAME = "Europe"
            self.CSHORT = "euro"
            self.latmin = 35.
            self.latmax = 75.
            self.lonmin = -15.
            self.lonmax = 45.
        elif self.no == 4:
            self.CNAME = "Asia"
            self.CSHORT = "asia"
            self.latmin = -50.
            self.latmax = 90.
            self.lonmin = 20.
            self.lonmax = 180.
        elif self.no == 5:
            self.CNAME = "Africa"
            self.CSHORT = "afri"
            self.latmin = -40.
            self.latmax = 40.
            self.lonmin = -20.
            self.lonmax = 60.
        elif self.no == 6:
            self.CNAME = "Globe_without_Antartica"
            self.CSHORT = "globa_wa"
            self.latmin = -60.
            self.yticks = [-60, -30, 0, 30, 60, 90]
        elif self.no == 7:
            self.CNAME = "Arctic"
            self.CSHORT = "arctic"
            self.latmin = 40.
            self.latmax = 90.
            self.lonmin = -180.
            self.lonmax = 180.
        elif self.no == 8:
            self.CNAME = "Baltic"
            self.CSHORT = "baltic"
            self.latmin = 45.
            self.latmax = 70.
            self.lonmin = 5.
            self.lonmax = 40.
        elif self.no == 9:
            self.CNAME = "Australia"
            self.CSHORT = "aus"
            self.latmin = -45.
            self.latmax = -10.
            self.lonmin = 112.
            self.lonmax = 154.
        elif self.no == 10:
            self.CNAME = "Europe-5 Min."
            self.CSHORT = "euro5min"
            self.latmin = 27.
            self.latmax = 72.
            self.lonmin = -11.
            self.lonmax = 69.
        elif self.no == 11:
            self.CNAME = "Burdekin"
            self.CSHORT = "burde"
            self.latmin = -23.
            self.latmax = -18.
            self.lonmin = 144.
            self.lonmax = 149.
        elif self.no == 12:
            self.CNAME = "Elbe_estuar"
            self.CSHORT = "elbe"
            self.latmin = 53.
            self.latmax = 57.
            self.lonmin = 7.
            self.lonmax = 11.
        elif self.no == 13:
            self.CNAME = "Rhine_delta"
            self.CSHORT = "rhein"
            self.latmin = 51.
            self.latmax = 54.
            self.lonmin = 3.
            self.lonmax = 7.
        elif self.no == 14:
            self.CNAME = "NE Bottnian Bay"
            self.CSHORT = "kiimi"
            self.latmin = 64.
            self.latmax = 67.
            self.lonmin = 24.
            self.lonmax = 28.
        elif self.no == 19:
            self.CNAME = "Ems"
            self.CSHORT = "ems"
            self.latmin = 51.5
            self.latmax = 54.
            self.lonmin = 7.
            self.lonmax = 9.
        elif self.no == 20:
            self.CNAME = "Europe-mHm"
            self.CSHORT = "euromhm"
            self.latmin = 35.
            self.latmax = 72.
            self.lonmin = -11.
            self.lonmax = 41.
        elif self.no == 21:
            self.CNAME = "Iceland"
            self.CSHORT = "iceland"
            self.latmin = 63.
            self.latmax = 67.
            self.lonmin = -25.
            self.lonmax = -13.
        elif self.no == 22:
            self.CNAME = "Turkey"
            self.CSHORT = "turkey"
            self.latmin = 30.
            self.latmax = 50.
            self.lonmin = 20.
            self.lonmax = 45.
        elif self.no == 23:
            self.CNAME = "Baltic_and_North_Sea"
            self.CSHORT = "baltic_nsea"
            self.latmin = 45.
            self.latmax = 70.
            self.lonmin = -5.
            self.lonmax = 35.
        elif self.no == 40:
            self.CNAME = "Mexico"
            self.CSHORT = "mexico"
            self.latmin = 13.
            self.latmax = 33.
            self.lonmin = -115.
            self.lonmax = -90.
        elif self.no == 50:
            self.CNAME = "Leichhardt"
            self.CSHORT = "leich"
            self.latmin = -23.
            self.latmax = -15.
            self.lonmin = 135.
            self.lonmax = 146.
        elif self.no == 51:
            self.CNAME = "Rio_Negro"
            self.CSHORT = "neg_a"
            self.latmin = -45.
            self.latmax = -37.
            self.lonmin = -72.
            self.lonmax = -62.
        elif self.no == 99:
            print("Self-defined region - no predefined bounds", self.no)
        else:
            print("Unknown region no. ", self.no)
            status = -1
        return status

    def set_region(self):

        # Select region number and retrieve region information
        # Routine returns status: 0 = ok, -1 = region no. not found

        status = 0
        print('\n   Latitude Axis: ', self.latmin, ' to ', self.latmax)
        print('  Longitude Axis: ', self.lonmin, ' to ', self.lonmax)
        print("Set new region no.")
        print("      99 = Define own region")
        print("       0 = Globe")
        print("       1 = South America")
        print("       2 = North America")
        print("       3 = Europe")
        print("       4 = Asia")
        print("       5 = Africa")
        print("       6 = Globe_without_Antartica")
        print("       7 = Arctic")
        print("       8 = Baltic")
        print("      10 = Europe-5 Min.")
        print("      21 = Iceland")
        print("      22 = Turkey")
        idum = input("Enter menu point number: ")  # Python 3
        self.no = int(idum)

        if self.no == 99:
            print("\nSet new region bounds for plotting")
            idum = input("Enter minimum latitude: ")
            self.latmin = float(idum)
            idum = input("Enter maximum latitude: ")
            self.latmax = float(idum)
            idum = input("Enter minimum longitude: ")
            self.lonmin = float(idum)
            idum = input("Enter maximum longitude: ")
            self.lonmax = float(idum)
        else:
            status = self.get_region()

        return status


# ***************************************************
class PlotOptions:

    def __init__(self):           # Default Initialization
        self.icoast = 1  
        self.igitt = 0 
        self.imask = 0 
        self.iloga = 0
        self.pdate = 1
        self.ispot = 0 ; self.ishp = 0 ; self.iriver = 0 ; self.iclass = 0
        self.icountry = 0
        self.ufak = 1.
        self.xmin = 0. ; self.xmax = 100. ; self.xlow = 0. ; self.ncol = 10
        self.icoltyp = 1
        self.cmap = get_colortype(self.icoltyp, False)
        self.ibarloc = 0
        self.ptitle = 'Test'
        self.ifile = 0 ; self.iout = 0 
        self.dninp = 'test.nc' ; self.dnout = 'test.pdf'

# ***************************************************

# ***************************************************
def plot_menu(opt, reg, lmap=False):
    # ***************************************************
    #
    # Menu for pl_map.py
    #
    lmenu = True
    while lmenu:
        print("\n      Plot title: ", opt.ptitle)
        print('   Latitude Axis: ', reg.latmin, ' to ', reg.latmax)
        print('  Longitude Axis: ', reg.lonmin, ' to ', reg.lonmax)
        if lmap:
            print('Value bounds set: ', opt.xmin, ' to ', opt.xmax)
        #    print ('   for data range ', valmin, ' to ', valmax)
        print(" ")
        print("  -1 = Exit")
        print("   1 ==> Make plot with current setting")
        if lmap:
            print("   2 ==> Set value bounds for plotting")
            print("   3 ==> Set amount of colors ncol = ", opt.ncol)
            print("   4 ==> Set threshold below which values are treated as missing values = ", opt.xlow)
            print("   5 ==> Change colormap: ", opt.cmap)
        print("   6 ==> Set region bounds for plotting")
        print("   8 ==> Change plot title ", opt.ptitle)
        print("   9 ==> Draw coast lines (0 = off): ", opt.icoast)
        #    print, "  10 ==> Daten auslesen"
        #    print ("  11 ==> Change character size of CO system ", BCHAR
        print("  12 ==> Change input filename: ", opt.dninp)
        print("  13 ==> Change output filename: ", opt.dnout)
        if opt.igitt == 0:
            print("  14 ==> Switch ON: Plotting coordinate grid lines")
        else:
            print("  14 ==> Switch OFF: Plotting coordinate grid lines")
        if opt.imask == 0:
            print("  15 ==> Switch ON: Plotting/Application of mask. e.g. land sea mask")
        else:
            print("  15 ==> Change or switch OFF: Plotting/Application of mask")
        if opt.ispot == 0:
            print("  16 ==> Switch ON: Plotting selected location spots")
        else:
            print("  16 ==> Change or switch OFF: Plotting selected location spots")
        if lmap:
            if opt.iloga == 0:
                print("  17 ==> Switch ON: Logarithmic plot")
            else:
                print("  17 ==> Change or switch OFF: Logarithmic plot, iloga = ", opt.iloga)
            print("  18 ==> Change unit factor: ", opt.ufak)
        if opt.iriver == 0:
            print("  19 ==> Switch ON: Plotting river network")
        else:
            print("  19 ==> Switch OFF: Plotting river network")
        if opt.pdate == 0:
            print("  20 ==> Switch ON: Add name and date to plot")
        else:
            print("  20 ==> Switch OFF: Add name and date to plot")
        if opt.ishp == 0:
            print("  21 ==> Switch ON: Plotting shapefile")
        else:
            print("  21 ==> Switch OFF: Plotting shapefile")
            print("  22 ==> Change shapefile name: ", opt.dn_shp)
            print("  23 ==> Change shapefile directory: ", opt.dir_shp)
            print("  24 ==> Change shapefile projection: EPSG = ", opt.id_epsg)

        #      print, "   4 ==> Eingabe einer Farbnummer fuer Weiss, z.Zt. = ", IWEISS
        #      print, "  26 ==> Achsen-Shift eingeben und ausfuehren, L= ", $
        #                       XLSHIFT, " B=", XBSHIFT
        #      print, "   1 ==> IDRUCK aendern, z.Zt. = ", IDRUCK
        #     IF (IZOOM EQ 0) THEN print, "  27 ==> Zooming ANschalten"
        #      IF (IZOOM EQ 1) THEN print, "  27 ==> Zooming AUSschalten"
        #      print, "  28 ==> Treatment of Missing values, z.Zt.: ", CMISS, XMISS

        idum = input("Enter menu point number: ")  # Python 3
        if idum == '':
            imenu = 0
        else:
            imenu = int(idum)

        if imenu == -1:
            exit()
        elif imenu == 1:
            lmenu = False
        elif imenu == 2 and lmap:
            print('   Current value bounds: ', opt.xmin, ' to ', opt.xmax)
            print("   Set new value bounds for plotting")
            idum = input("Enter lower bound: ")  # Python 3
            opt.xmin = float(idum)
            idum = input("Enter upper bound: ")  # Python 3
            opt.xmax = float(idum)
        elif imenu == 3 and lmap:
            print('   Current Number of colors = ', opt.ncol)
            idum = input('Enter new number of colors: ')
            opt.ncol = int(idum)
        elif imenu == 4 and lmap:
            print('   Current Threshold xlow = ', opt.xlow)
            print("   Set threshold for which values below are treated as mssing values")
            idum = input("Enter threshold: ")  # Python 3
            opt.xlow = float(idum)
        elif imenu == 5 and lmap:  # set color map
            opt.cmap = get_colortype(opt.icoltyp, lmenu)
        elif imenu == 6:  # set_region
            status = -1
            while status == -1:
                status = reg.set_region()
        elif imenu == 8:
            opt.ptitle = input("Enter new plot title: ")
        elif imenu == 9:
            if opt.icoast == 0:
                opt.icoast = 1
            else:
                opt.icoast = 0
        elif imenu == 12:
            opt.dninp = input("Enter new input file name: ")
        elif imenu == 13:
            opt.dnout = input("Enter new output file name: ")
        elif imenu == 14:
            if opt.igitt == 0:
                opt.igitt = 1
            else:
                opt.igitt = 0
        elif imenu == 15:
            if opt.imask == 0:
                opt.imask = 1
            else:
                opt.imask = 0
        elif imenu == 16:
            opt.ispot = get_spotlist(opt)
        elif imenu == 17 and lmap:
            if opt.iloga == 0:
                opt.iloga = 1
            else:
                opt.iloga = 0
        elif imenu == 19:
            if opt.iriver == 0:
                opt.iriver = 1
            else:
                opt.iriver = 0
        elif imenu == 20:
            if opt.pdate == 0:
                opt.pdate = 1
            else:
                opt.pdate = 0
        elif imenu == 21:
            if opt.ishp == 0:
                opt.ishp = 1
            else:
                opt.ishp = 0
        elif imenu == 22:
            print('\033[1m' + " Content              Shapefile" + '\033[0m')
            print(" Australia            aus/australia_basins  aus/australia_streams  aus/rbasin_polygon")
            print(" Canada               na/watershed_p_v2  na/canada_basins  na/canada_water")
            print(" CCM21 data           river_ or cat_balticum _caucasus _dnj _donau _neurope _weurope _see")
            print(" Diva data            Balkan_water baltic_states_ Caucasus_TR_ B_D_NL_ ITA_ rus_europe_ Ukraine_")
            print(" Greek rivers         potamoi_eper_gr")
            print(" Finnish catchments   Paajako")
            print(" Finnish rivers       Uoma10 (BIG, extended until mouth) or VHSJoki2016 (without lakes)")
            print(" Irish rivers         irl/WATER_RiverBasin & irl/WATER_RivNetRoutes ")
            print(" Italian rivers       waterways ")
            print(" Mexican rivers       sa/Rios_principales_2020")
            print(" Norway cat/Rivers    Nedborfelt_Vassdragsomr  Elv_Hovedelv")
            print(" Spanish rivers       Rius_CARACT")
            print(" Swedish cat/rivers   smhi_huvudavrinningsomraden_SVAR_2012_2  smhi_vattenforekomster_vattendrag_SVAR2012_2")
            print(" UK catchments/rivers ccm21/cat_uk  uk/uk_watercourse")
            print(" US rivers            na/watershed_p_v2  na/us_rivers  na/us_major_rivers")
            opt.dn_shp = input("    Enter new shapefile name without .shp: ")
            if opt.dn_shp == "na/watershed_p_v2":
                opt.id_epsg = 2163
            elif opt.dn_shp == "na/us_rivers" or opt.dn_shp == "na/canada_water" or opt.dn_shp == "na/canada_basins":
                opt.id_epsg = 4169
            elif opt.dn_shp == "sa/Rios_principales_2020":
                opt.id_epsg = 6362
 
        elif imenu == 23:
            opt.dir_shp = input("Enter new shapefile directory: ")
        elif imenu == 24:
            print("Set new region no.")
            print("       0 = Standard lat/lon = GCS_WGS_1984 = 4326 ")
            print("    2100 = Greek_Grid -> Greek rivers")
            print("    2163 = Sphere_ARC_INFO_Lambert_Azimuthal_Equal_Area --> na/watershed_p_v2")
            print("    3006 = SWEREF99 TM -> Swedish rivers")
            print("    3035 = ETRS_1989_LAEA_L52_M10 -> ERC catchments")
            print("    3043 = ETRS_1989_UTM_Zone_31N -> Spanish rivers")
            print("    3047 = EUREF_FIN_TM35FIN -> Finnish rivers")
            print("    4202 = Australian Geodetic Datum 1966 -> Australia rbasin")
            print("    4269 = GCS_North_American_1983 -> North America (Canada/USA)")
            print("    4283 = GCS_GDA_1994 -> Australia")
            print("    4326 = GCS_WGS_1984, e.g. CCM21, Norway & Italy rivers, Natural Earth data, Mex. cat.")
            print("    6362 = CCL_ITRF_1992 -> Mexican rivers")
            print("   27700 = UK rivers")
            print("   29900 = TM65_Irish_Grid --> Irish rivers")
            idum = input("Enter shapefile projection EPSG no.: ")
            opt.id_epsg = int(idum)
        else:
            print('\n Option does not exists! \n')
        #
    #  end of while loop
    return


# ***************************************************
def draw_mask_borders_map(fig, lmask, x0, x1, y0, y1):
    # ***************************************************
    #
    # Draw the outline of a (0/1)-mask, e.g. a land sea mask.
    # Boolean mask is provided as lmask
    #
    print('size of lmask: ', len(lmask), ' x ', len(lmask[0]))

    ver_seg = np.where(lmask[:, 1:] != lmask[:, :-1])  # vertical segments
    hor_seg = np.where(lmask[1:, :] != lmask[:-1, :])  # horizontal segments

    # if we have a horizontal segment at 7,2, it means that it must be drawn between pixels
    #   (2,7) and (2,8), i.e. from (2,8)..(3,8)
    # in order to draw a discountinuous line, we add Nones in between segments
    j = []
    for p in zip(*hor_seg):
        j.append((p[1], p[0] + 1))
        j.append((p[1] + 1, p[0] + 1))
        j.append((np.nan, np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        j.append((p[1] + 1, p[0]))
        j.append((p[1] + 1, p[0] + 1))
        j.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(j)

    # now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points
    segments[:, 0] = x0 + (x1 - x0) * segments[:, 0] / lmask.shape[1]
    segments[:, 1] = -y0 - (y1 - y0) * segments[:, 1] / lmask.shape[0]

    print('Lon mask segments Min ', np.nanmin(segments[:, 0]), ' Max ', np.nanmax(segments[:, 0]))
    print('Lat mask segments Min ', np.nanmin(segments[:, 1]), ' Max ', np.nanmax(segments[:, 1]))

    fig.plot(segments[:, 0], segments[:, 1], color='gray', linewidth=1)


# ***************************************************
def draw_mask_borders(fig, lmask, x0, x1, y0, y1):
    # ***************************************************
    #
    # Draw the outline of a (0/1)-mask, e.g. a land sea mask.
    # Boolean mask is provided as lmask
    #
    print('size of lmask: ', len(lmask), ' x ', len(lmask[0]))

    ver_seg = np.where(lmask[:, 1:] != lmask[:, :-1])  # vertical segments
    hor_seg = np.where(lmask[1:, :] != lmask[:-1, :])  # horizontal segments

    # if we have a horizontal segment at 7,2, it means that it must be drawn between pixels
    #   (2,7) and (2,8), i.e. from (2,8)..(3,8)
    # in order to draw a discountinuous line, we add Nones in between segments
    j = []
    for p in zip(*hor_seg):
        j.append((p[1], p[0] + 1))
        j.append((p[1] + 1, p[0] + 1))
        j.append((np.nan, np.nan))

    # and the same for vertical segments
    for p in zip(*ver_seg):
        j.append((p[1] + 1, p[0]))
        j.append((p[1] + 1, p[0] + 1))
        j.append((np.nan, np.nan))

    # now we transform the list into a numpy array of Nx2 shape
    segments = np.array(j)

    # now we need to know something about the image which is shown
    #   at this point let's assume it has extents (x0, y0)..(x1,y1) on the axis
    #   drawn with origin='lower'
    # with this information we can rescale our points
    segments[:, 0] = x0 + (x1 - x0) * segments[:, 0] / lmask.shape[1]
    #  segments[:,1] = -y0 - (y1-y0) * segments[:,1] / lmask.shape[0]
    segments[:, 1] = y0 + (y1 - y0) * segments[:, 1] / lmask.shape[0]

    print('Lon mask segments Min ', np.nanmin(segments[:, 0]), ' Max ', np.nanmax(segments[:, 0]))
    print('Lat mask segments Min ', np.nanmin(segments[:, 1]), ' Max ', np.nanmax(segments[:, 1]))

    fig.plot(segments[:, 0], segments[:, 1], color='gray', linewidth=1)


# ***************************************************
def get_spotlist(opt):
    # ***************************************************

    print('\n   Define Spot List or switch off plotting spots')

    if hasattr(opt, 'xspot'):
        if opt.ispot > 0:
            for i in range(len(opt.xspot)):
                print(i + 1, '. Spot:', opt.xspot[i], opt.yspot[i])
        print("\nType desired acion:")
        print("   0 = Switch OFF plotting spots")
        print("   1 = Define new spot list")
        print("   2 = Append to existing spot list")
        print("   3 = Use and plot existing spot list")
        iaction = input("Enter menu point number: ")  # Python 3
        idum = int(iaction)
    else:
        idum = 1

    ispot = 1
    if idum == 0:
        ispot = 0
    elif idum == 1 or idum == 2:
        if idum == 1:
            opt.xspot = []
            opt.yspot = []
        loopspot = True
        print("Set new spot coordinates for plotting (-999 0 for exit)")
        while loopspot:
            cnspot = str(len(opt.xspot) + 1)
            xdum, ydum = input(cnspot + ". Spot Enter <longitude latitude>: ").split()
            if float(xdum) == -999.:
                loopspot = False
            if loopspot:
                opt.xspot.append(float(xdum))
                opt.yspot.append(float(ydum))
    elif idum == 3:
        print('Existing Spotlist will be used!')
    else:
        print('Menu point does not exist --> Plotting Switched OFF ')
        ispot = 0

    return ispot


# ***************************************************
def draw_circles(fig, x, y):
    # ***************************************************
    #
    # Draw circles for every pair of x and y coordinates
    #
    print('size of x: ', len(x), ' and y ', len(y))

    fig.plot(x, y, 'mo', fillstyle='none')


#
# ***************************************************
def mask_classes(fdata, classlist, iclass):
    # ***************************************************
    #
    # Masks data array and keeps only those numbers that are included in classlist
    #
    zeps = 1.E-10
    icl = 0
    masked_array = None

    for cl in classlist:
        mask1 = (abs(fdata - float(cl)) < zeps)
        if icl == 0:
            if iclass == 2:
                masked_array = np.where(mask1, fdata, 0)
            elif iclass == 3:
                masked_array = np.where(mask1, icl + 0.5, 0)
        else:
            if iclass == 2:
                masked_array = np.where(mask1, fdata, masked_array)
            elif iclass == 3:
                masked_array = np.where(mask1, icl + 0.5, masked_array)
        icl += 1

    return masked_array


# ***************************************************
def draw_shapefile(fig, dn_shape,
                   id_epsg=3047):
    # ***************************************************
    #
    # Reads ad draws the content of the shapfile dn_shape.shp
    # Note that dn_shape comprises the path+filename of the shapefile without extension.
    #
    import shapefile
    from pyproj import Proj, transform

    #
    # *** 1) read with PyShp
    sf = shapefile.Reader(dn_shape + ".shp", encodingErrors="replace")  # encoding="ascii", utf-8, ansi,...
    input_projection = None
    output_projection = None

    if id_epsg > 0 and id_epsg != 4326:
##        input_projection = Proj(init="epsg:" + str(id_epsg))  # Finnish catchment shapefile projection: 3047
##        output_projection = Proj(init="epsg:4326")  # Normal geogrpahical lon/lat projection
        input_projection = Proj("epsg:" + str(id_epsg))  # Finnish catchment shapefile projection: 3047
        output_projection = Proj("epsg:4326")  # Normal geogrpahical lon/lat projection


        x = None
        y = None
        for shape in sf.shapeRecords():
            i = 0
            i_start = shape.shape.parts[i]
            if i == len(shape.shape.parts) - 1:
                i_end = len(shape.shape.points)
            else:
                i_end = shape.shape.parts[i + 1]
            x = [i[0] for i in shape.shape.points[i_start:i_end]]
            y = [i[1] for i in shape.shape.points[i_start:i_end]]

        lon, lat = transform(input_projection, output_projection, x, y, always_xy=True)
        print(x[0], y[0])
        print(' --> ', float(lon[0]), float(lat[0]))

    # Note that the propertie "parts" of a shape returns the starting indexes of differents geometries inside a feature.
    for shape in sf.shapeRecords():
        for i in range(len(shape.shape.parts)):
            i_start = shape.shape.parts[i]
            if i == len(shape.shape.parts) - 1:
                i_end = len(shape.shape.points)
            else:
                i_end = shape.shape.parts[i + 1]
            x = [i[0] for i in shape.shape.points[i_start:i_end]]
            y = [i[1] for i in shape.shape.points[i_start:i_end]]
            if id_epsg > 0 and id_epsg != 4326:
                lon, lat = transform(input_projection, output_projection, x, y, always_xy=True)
            else:
                lon = x
                lat = y

            fig.plot(lon, lat, linewidth=0.5, color='red')

    print("Shapefile " + dn_shape + " drawn")

# ***************************************************
def plot_mapping_arrows(xlon, xlat, xlon_dst, xlat_dst, reg, cproj):
# ***************************************************
    """
    # Plot arrows from Source coordinates (1D arrays xlon, xlat) to Destination coordinates (xlon_dst, xlat_dst)
    #   for all coordinates that are located within the target region reg.
    """
    #
    # Define arrows
    x = []
    y = []
    dx = []
    dy = []

    narr = len(xlon)
    iarrow = -1
    for i in range(narr):
        if 41 < xlat[i] < 43 and 10 < xlon[i] < 12: 
            print(i, 'plot_mapping: ', xlon[i], xlat[i], ' -> ', xlon_dst[i], xlat_dst[i]) 

        if reg.lonmin <= xlon[i] <= reg.lonmax and reg.latmin <= xlat[i] <= reg.latmax and    \
                reg.lonmin <= xlon_dst[i] <= reg.lonmax and reg.latmin <= xlat_dst[i] <= reg.latmax:
            iarrow += 1

            x.append(xlon[i])
            dx.append(xlon_dst[i]-xlon[i])
            y.append(xlat[i])
            dy.append(xlat_dst[i]-xlat[i])

        #
    if iarrow == -1:
        print('*** plot_mapping_arrows: No mapping arrow within selected region!')
        print('*** Selected Region: Xlon', reg.lonmin, ' - ', reg.lonmax, ', xlat: ',  reg.latmin, ' - ', reg.latmax)
        return

    print('SRC Region: xlon', min(x), ' - ', max(x), ', xlat: ',  min(y), ' - ', max(y))
    plt.quiver(x, y, dx, dy, units='xy', headlength=2, headaxislength=2, color='blue',
               angles='xy', scale_units='xy', scale=1, transform=cproj)

# ***************************************************
class Iconmap(Plotter):
# ***************************************************
#
    # Specify the defaults
    rc = {
          "cmap": "Reds",
          "norm": None,
          "transform": "cf",
          "plot": "poly",
          "miss_color": "brown",
         }

    # ***************************************************
    def convert_coordinate(self, coord, *variables):
        if coord.attrs.get("units", "").startswith("radian") or any(
            var.attrs.get('units', '').startswith('radian')
            for var in variables
        ):
            coord = coord.copy(data=coord * 180. / np.pi)
            coord.attrs["units"] = "degrees"
        return coord

    # Specify the formatoptions
    cmap = CMap("cmap", bounds="norm")
    norm = Bounds("norm")
    transform = Transform("transform")
    miss_color = MissColor("miss_color")

    plot = MapPlot2D("plot", bounds="norm")

# ***************************************************
def plot_icon(data, ax, **formatoptions):
    return Iconmap(data, ax=ax, **formatoptions)

