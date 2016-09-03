import healpy as hp
import numpy as n,os,sys
import sys,optparse
from pylab import *
from scipy import optimize
from ECHO.plot_utils import get_interp_val
from ECHO.read_utils import read_map
style.use('./echo.mplstyle')

o = optparse.OptionParser()
o.set_description('plot slices of healpix maps')
o.add_option('--savefig',type=str,help='output figure name')
#o.add_option('--exp',action='store_true',help='do experimental stuff')
opts,args = o.parse_args(sys.argv[1:])

filenames = ['../data/acc_MWAtile_2016_EWtx_EWrx_8_beam-bicolog_legs_360.fits']

#rx_model = '../data/ff137_5_xpol.fits'  #old bradley model?
#rx_model = '../data/bradley_amr.fits'  #rerun of the bradley model using amr
rx_model = '../data/MWA_Tile_nside32_X.fits'



#cross cut coordinates
print "plotting E and H plane slices"
alt = n.linspace(-n.pi/2,n.pi/2)
az = n.zeros_like(alt)

#get the rx model
rx_map = read_map(rx_model)
rx_map -= np.mean(rx_map[:3])
rx_slice_E = get_interp_val(rx_map,alt,az)
rx_slice_H = get_interp_val(rx_map,alt,az+np.pi/2)

fig,axarr = subplots(1,2,figsize=(10,6),sharey=True,sharex=True)
#plot the total model
axarr[0].plot(alt*180/np.pi,rx_slice_E,'-k',lw=2)
axarr[1].plot(alt*180/np.pi,rx_slice_H,'-k',lw=2)

E_color = 'k'
H_color = 'k'
ls = ['-',':','--','-.']
for i,filename in enumerate(filenames):
    rot=0
    if filename.find('NS')!=-1:rot=np.pi/2
    M = read_map(filename)
    M -= np.ma.mean(M[:3])
    M_slice_E = get_interp_val(M,alt,az+rot)
    M_slice_H = get_interp_val(M,alt,az+np.pi/2+rot)

    err_filename = filename.replace('beam','rms')
    ME = read_map(err_filename)*2
    ME_slice_E = get_interp_val(ME,alt,az+rot)
    ME_slice_H = get_interp_val(ME,alt,az+np.pi/2+rot)
    axarr[0].errorbar(alt*180/np.pi,M_slice_E,yerr=ME_slice_E,
            color=E_color,fmt=ls[i])
    axarr[1].errorbar(alt*180/np.pi,M_slice_H,yerr=ME_slice_H,
            color=H_color,fmt=ls[i])
#legend(loc='best')
axarr[0].set_ylabel('[dB V/m]')
axarr[0].set_xlabel('$\\theta$ [deg]')
axarr[0].set_title('E plane')
axarr[1].set_title('H plane')
axarr[1].set_xlabel('$\\theta$ [deg]')
axarr[0].grid()
axarr[1].grid()

if not opts.savefig is None:
    print "plotting to", opts.savefig
    savefig(opts.savefig)
else:
    show()
