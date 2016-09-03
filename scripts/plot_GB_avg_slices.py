import healpy as hp
import numpy as n
import sys,optparse
from pylab import *
from scipy import optimize
from ECHO.plot_utils import get_interp_val
from ECHO.read_utils import read_map
from matplotlib import rcParams
style.use('./echo.mplstyle')
o = optparse.OptionParser()
o.set_description('plot slices of GB models on averaged data ')
opts,args = o.parse_args(sys.argv[1:])


rmsfiles = ['../data/GB_rx_model_rms_NS.fits',
            '../data/GB_rx_model_rms_EW.fits',
            '../data/GB_rx_model_rms.fits']
beamfiles = ['../data/GB_rx_model_beam_NS.fits',
    '../data/GB_rx_model_beam_EW.fits',
    '../data/GB_rx_model_beam.fits']




fig,axarr = subplots(1,3,sharey=True,sharex=True,figsize=(12,4))
labels=['NS','EW','all']
#cross cut coordinates
print "plotting E and H plane slices"
alt = n.linspace(-n.pi,n.pi,num=100)
az = n.zeros_like(alt)

#load the model
modelfile = '../data/bradley_amr_nside32.fits'
B = read_map(modelfile)
B_slice_E = get_interp_val(B,alt,az)
B_slice_H = get_interp_val(B,alt,az+np.pi/2)


ls = ['-',':','--']
EHcolors = [['r','b'],['b','r'],['b','r']]
for i,(filename,rmsfile) in enumerate(zip(beamfiles,rmsfiles)):
    M = read_map(filename)
    M_slice_E = get_interp_val(M,alt,az)
    M_slice_H = get_interp_val(M,alt,az+np.pi/2)

    M_err = read_map(rmsfile)
    M_slice_E_err = get_interp_val(M_err,alt,az)
    M_slice_H_err = get_interp_val(M_err,alt,az+np.pi/2)
    E_color,H_color = EHcolors[i]
    axarr[i].errorbar(alt*180/np.pi,M_slice_E,yerr=M_slice_E_err,color=E_color,fmt='-',label=labels[i])
    axarr[i].errorbar(alt*180/np.pi,M_slice_H,yerr=M_slice_H_err,color=H_color,fmt='-')
    #plot the model
    axarr[i].plot(alt*180/np.pi,B_slice_E,':',color='b',label='model')
    axarr[i].plot(alt*180/np.pi,B_slice_H,':',color='r')

    axarr[i].set_xlabel('$\\theta$ [deg]')
#legend(loc='best')
axarr[0].set_ylabel('dB V/m')
axarr[0].set_xlim(-80,80)
axarr[0].set_ylim(-25,3)
for ax in axarr.ravel():
    ax.grid()
tight_layout()
savefig('../figures/GB_model_avg_slice.png')
