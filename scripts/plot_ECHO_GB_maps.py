#! /usr/bin/env python
from matplotlib.pyplot import *
from matplotlib import gridspec
import numpy as np
import sys,os,optparse
import healpy as hp
from ECHO.plot_utils import project_healpix
from ECHO.read_utils import read_map
o = optparse.OptionParser()
o.set_description('make quad plot of orbcomm beams')
o.add_option('--max',type=float,help='max plotted value')
o.add_option('--min',type=float,help='min plotted value')
o.add_option('--savefig',default=None,type=str)
o.add_option('--mode',type='str',
    help="'usecal'= use the calibrated maps, 'err'= plot the rms files ")
opts,args = o.parse_args(sys.argv[1:])



files = ['../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits']
if opts.mode=='usecal':
    for i in xrange(len(files)):
        files[i] = files[i].replace('.fits','-bicolog_legs_360.fits')
elif opts.mode=='err':
    for i in xrange(len(files)):
        files[i] = files[i].replace('beam','rms')
#contour labeling borrowed from
#http://matplotlib.org/examples/pylab_examples/contour_label_demo.html
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
fmt = '%r$^\circ$'
fig=figure()
rows = int(np.ceil(np.sqrt(len(files))))
cols = int(np.ceil(len(files)/rows))
#gs1 = gridspec.GridSpec(rows, cols)
#gs1.update(bottom=0.2)
axs = []
for i,filename in enumerate(files):
    print i
    #ax = subplot(gs1[i])
    ax = subplot(rows,cols,i+1)
    axs.append(ax)
    M = read_map(filename)
    THETA,PHI,M_proj = project_healpix(M)
    imax =imshow(M_proj)
    CS = contour(THETA*180/np.pi,[20,40,60],colors='k')
    # Recast levels to new class
    CS.levels = [nf(val) for val in CS.levels]
    clabel(CS, inline=1, fontsize=10,fmt=fmt)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
#tight_layout()
subplots_adjust(bottom=0.15,wspace=0)
cax1=axes([.14,0.09,0.75,0.03])
colorbar(imax,cax=cax1,orientation='horizontal',label='dB')
axs[0].set_title('Calibrator A')
axs[1].set_title('Calibrator B')
axs[0].set_ylabel('North-South pol')
axs[2].set_ylabel('East-West pol')
if not opts.savefig is None:
    savefig(opts.savefig)
