import numpy as np
import healpy as hp
import math
from healpy import _healpy_pixel_lib as pixlib
from matplotlib.collections import PolyCollection
from matplotlib import cm

from position_utils import latlon2xy,to_spherical
from time_utils import gps_to_HMS,find_peak

def make_beam(lats,lons,alts,spec_raw,lat0=0.0,lon0=0.0,
              nsides=8,volts=False,normalize=False,freq_chan=0):
    # Convert lat/lon to x/y
    x,y = latlon2xy(lats,lons,lat0,lon0)
    # Obtain spherical coordinates for x, y, and alt
    rs,thetas,phis = to_spherical(x,y,alts)

    #freq_chan = np.argmax(spec_raw[0,:])
    z = spec_raw[:,freq_chan].copy()
    if volts:
        # Distance normalization
        r0 = 100 # reference position for distance normalization (unit: meters)
        z = 10*np.log10((2*z**2)*(rs/r0)**2)
    if normalize:
        z -= z.max() # Scaled on [-infty,0]

    # Set binsize (used in function grid_data)
    # Affects the apparent size of the pixels on the plot created below.
    binsize=5
    # Obtain gridded data
    grid,bins,rmsBins,binloc,xg,yg,gcounts,grms = grid_data(x,y,z,binsize=binsize)

    # Healpix things
    nPixels = hp.nside2npix(nsides)
    hpx_beam = np.zeros(nPixels)
    hpx_counts = np.zeros(nPixels)
    hpx_rms = np.zeros(nPixels)
    # Find pixel # for a given theta and phi
    pixInd = hp.ang2pix(nsides,thetas,phis,nest=False)
    # Set pixel values at pixInd to power values
    hpx_beam[pixInd] = z
    hpx_counts[pixInd] = gcounts
    hpx_rms[pixInd] = grms
    # Grey out pixels with no measurements
    hpx_beam[hpx_beam == 0] = np.nan
    hpx_counts[hpx_counts == 0] = np.nan
    hpx_rms[hpx_rms == 0] = np.nan

    return np.ma.masked_invalid(hpx_beam),\
           np.ma.masked_invalid(hpx_counts),\
           np.ma.masked_invalid(hpx_rms)


def grid_data(x, y, z, binsize=0.01, retbin=True, retloc=True, retrms=True):
    # Get extrema values.
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    # Make coordinate arrays.
    xi = np.arange(xmin, xmax+binsize, binsize)
    yi = np.arange(ymin, ymax+binsize, binsize)
    xi, yi = np.meshgrid(xi,yi)

    # Make the grid.
    grid = np.zeros(xi.shape, dtype=x.dtype)
    nrow, ncol = grid.shape
    if retbin: bins = np.copy(grid)
    if retrms: rmsBins = np.copy(grid)

    # Make arrays to store counts/rms for Healpix data
    gcounts = np.zeros_like(z)
    grms = np.zeros_like(z)

    # Create list in same shape as grid to store indices
    if retloc:
        wherebin = np.copy(grid)
        wherebin = wherebin.tolist()
    # Fill in the grid.
    for row in range(nrow):
        for col in range(ncol):
            xc = xi[row, col]    # x coordinate.
            yc = yi[row, col]    # y coordinate.

            # Find the position that xc and yc correspond to.
            posx = np.abs(x - xc)
            posy = np.abs(y - yc)
            ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)[0]
            ind  = np.where(ibin == True)[0]

            # Fill the bin.
            bin = z[ibin]
            gcounts[ibin] = np.sum(ibin) # Update counts at each x position
            if retloc: wherebin[row][col] = ind
            if retbin: bins[row, col] = bin.size
            if bin.size != 0:
                binval = np.median(bin)
                grid[row, col] = binval
                if retrms:
                    rmsBins[row,col] = np.std(bin)
                    grms[ibin] = np.std(bin)
            else:
                grid[row, col] = np.nan   # Fill empty bins with nans.
                if retrms:
                    rmsBins[row,col] = np.nan
                    grms[ibin] = np.nan

    # Return the grid
    if retbin:
        if retloc:
            if retrms:
                return np.ma.masked_invalid(grid), bins, np.ma.masked_invalid(rmsBins),\
                            wherebin, xi, yi, gcounts, grms
            else:
                return np.ma.masked_invalid(grid), bins, wherebin, xi, yi, gcounts, grms
        else:
            return np.ma.masked_invalid(grid), bins, xi, yi, gcounts, grms
    else:
        if retloc:
            if retrms:
                return np.ma.masked_invalid(grid), np.ma.masked_invalid(rmsBins),\
                            wherebin, xi, yi, gcounts, grms
            else:
                return np.ma.masked_invalid(grid), wherebin, xi, yi, gcounts, grms
        else:
            if retrms:
                return np.ma.masked_invalid(grid), np.ma.masked_invalid(rmsBins),\
                            xi, yi, gcounts, grms
            else:
                return np.ma.masked_invalid(grid), xi, yi, gcounts, grms


def animate_spectrum(i,spec_plot,spec_line,spec_raw):
    spec_line.set_ydata(spec_raw[i,:])
    #spec_plot.set_ylim([spec_raw[i,:].min(),spec_raw[i,:].max()])
    return


def make_polycoll(hpx_beam,plot_lim=[-90,-50],nsides=8):
    pix = np.where(np.isnan(hpx_beam)==False)[0]
    boundaries = hp.boundaries(nsides,pix)
    verts = np.swapaxes(boundaries[:,0:2,:],1,2)
    coll = PolyCollection(verts, array=hpx_beam[np.isnan(hpx_beam)==False],\
                                    cmap=cm.gnuplot,edgecolors='none')
    return coll


def get_interp_val(m,theta,phi,nest=False):
    m2=m.ravel()
    nside=hp.pixelfunc.npix2nside(m2.size)
    if nest:
        r=pixlib._get_interpol_nest(nside,theta,phi)
    else:
        r=pixlib._get_interpol_ring(nside,theta,phi)
    p = np.array(r[0:4])
    w = np.array(r[4:8])
    w = np.ma.array(w)
    w.mask = m2[p].mask
    del r
    return np.ma.sum(m2[p]*w/np.ma.sum(w,0),0)


def add_diagram(axs,xys,xytexts,colors,labels=None):
    for i in range(len(xys)):
        ax = axs[i]
        xy = xys[i]
        xytext = xytexts[i]
        color = colors[i]
        ax.annotate("",
                    xy=xy,
                    xycoords='axes fraction',
                    xytext=xytext,
                    textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='-',
                                    color=color,
                                    shrinkA=5, shrinkB=5,
                                    patchA=None,
                                    patchB=None,
                                    connectionstyle="arc3,rad=0.",
                                    ),
                    size=20)
        if not labels is None:
            label = labels[i]
            xyl = xy[:]
            if i == 0:
                xyl[0] += 0.05
            elif i == 1:
                xyl[1] += 0.05
            elif i == 2:
                xyl = (xy[0]+0.2,
                      xyl[1]+0.05)
            ax.annotate(label,
                        xy=xyl,
                        xycoords='axes fraction',
                        xytext=xyl,
                        textcoords='axes fraction',
                        arrowprops=dict(arrowstyle="-",
                                        color=color,
                                        shrinkA=5, shrinkB=5,
                                        patchA=None,
                                        patchB=None,
                                        connectionstyle="arc3,rad=0.",
                                        alpha=0.0
                                        ),
                        size=16)
