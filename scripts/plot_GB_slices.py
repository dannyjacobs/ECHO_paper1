import healpy as hp
import numpy as n,os,sys
import sys,optparse
from pylab import *
from scipy import optimize
from ECHO.plot_utils import get_interp_val,add_cut_glyph
from ECHO.read_utils import read_map
style.use('./echo.mplstyle')
o = optparse.OptionParser()
o.set_description('plot slices of GB healpix maps')
o.add_option('--savefig',type=str,help='output figure name')
o.add_option('--exp',action='store_true',help='do experimental stuff')
opts,args = o.parse_args(sys.argv[1:])

# filenames = ['../data/8_power_Nant_EWdipole_EWtransmitter.fits',
#             '../data/8_power_Sant_EWdipole_EWtransmitter.fits',
#             '../data/8_power_Nant_NSdipole_NStransmitter.fits',
#             '../data/8_power_Sant_NSdipole_NStransmitter.fits']
filenames = ['../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits']
#labels = ['Flight 8,','Flight 12,','Flight 14,','Flight 15,']
labels = ['station A, EW pol','station B, NS pol','station B, EW pol','station A, NS pol']
rots = {'../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits':0,
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits':0.367,
'../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits':6.63,
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits':4.89}
#rx_model = '../data/ff137_5_xpol.fits'  #old bradley model?
rx_model = '../data/bradley_amr.fits'  #rerun of the bradley model using amr
tx_model = '../data/bicolog_legs_360.fits'
rxtx_model = '../data/ECHO_ORBCOMM_model.fits'

#Try a short dipole model.

#cross cut coordinates
print "plotting E and H plane slices"
alt = n.linspace(-n.pi/2,n.pi/2)
az = n.zeros_like(alt)

#get the rx model
rx_map = read_map(rx_model)
rx_map -= np.mean(rx_map[:3])
rx_slice_E = get_interp_val(rx_map,alt,az)
rx_slice_H = get_interp_val(rx_map,alt,az+np.pi/2)

#get the tx model
tx_map = read_map(tx_model)
tx_map -= np.mean(tx_map[:3])
tx_slice_E = get_interp_val(tx_map,alt,az)
tx_slice_H = get_interp_val(tx_map,alt,az+np.pi/2)

rxtx_map = read_map(rxtx_model)
rxtx_slice_E = get_interp_val(rxtx_map,alt,az)
rxtx_slice_H = get_interp_val(rxtx_map,alt,az+np.pi/2)


#altfig = figure(figsize=(8,6))
fig,axarr = subplots(2,2,figsize=(10,6),sharey=True,sharex=True)

# todo:
# plot E and H in seperate panels
# use colors for data
# subtract the tx model
# plot the rx model

"""
filenames = ['../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits']
"""

E_color = 'k'
H_color = 'k'
ls = ['s','d','d','s'] #file order is north, south, south, north
for i,filename in enumerate(filenames):
    rot=0  #rotate the NS maps to our EW model frame
    if filename.find('NS')!=-1:rot=np.pi/2
    #apply rotation calibrations
    if opts.exp: rot += rots[filename]*np.pi/180
    #load the tx subtracted model
    files = []
    if True: #use the fitted transmitter map
        if filenames[i].find('Nant_NS')!=-1:
            beam_filename = filenames[i].replace('.fits','-137MHz_tx_farfield_rot10deg_dist9cm_rot.fits')
        else:
            beam_filename = filenames[i].replace('.fits','-137MHz_tx_farfield_dist6.8cm_rot.fits')
    else:
        beam_filename = filenames[i].replace('.fits','-bicolog_legs_360_rot.fits')
    #beam_filename = filename.replace('.fits','-bicolog_legs_360.fits')
    M = read_map(beam_filename)
    M -= np.ma.mean(M[:3])
    M_slice_E = get_interp_val(M,alt,az+rot)
    M_slice_H = get_interp_val(M,alt,az+np.pi/2+rot)

    err_filename = filename.replace('beam','rms')

    ME = read_map(err_filename)*2
    ME_slice_E = get_interp_val(ME,alt,az+rot)
    ME_slice_H = get_interp_val(ME,alt,az+np.pi/2+rot)
    if rot==0:
        col=1 #NS
    else:
        col=0 #EW
    axarr[0,col].errorbar(alt*180/np.pi,M_slice_E,yerr=ME_slice_E,
            label=labels[i],color=E_color,fmt=ls[i])
    axarr[1,col].errorbar(alt*180/np.pi,M_slice_H,yerr=ME_slice_H,
            color=H_color,fmt=ls[i])

#plot the rx model
axarr[0,1].plot(alt*180/np.pi,rx_slice_E,'-k',lw=2)
axarr[0,0].plot(alt*180/np.pi,rx_slice_E,'-k',lw=2)
axarr[1,1].plot(alt*180/np.pi,rx_slice_H,'-k',lw=2)
axarr[1,0].plot(alt*180/np.pi,rx_slice_H,'-k',lw=2)
#legend(loc='best')
axarr[0,0].set_title('NS')
axarr[0,1].set_title('EW')
axarr[0,0].set_ylabel('E plane\n [dB V/m]')
axarr[1,0].set_ylabel('H plane\n [dB V/m]')
axarr[1,0].set_xlabel('$\\theta$ (deg)')
axarr[1,1].set_xlabel('$\\theta$ (deg)')

subplots_adjust(wspace=0,hspace=0)
#add the E/H polarization slice glyph
add_cut_glyph(parent_axes=axarr[0,0],parent_fig=fig,cut='NS',pol='NS')  #E plane EW pol
add_cut_glyph(parent_axes=axarr[1,0],parent_fig=fig,cut='EW',pol='NS')  #H plane EW pol
add_cut_glyph(parent_axes=axarr[0,1],parent_fig=fig,cut='EW',pol='EW')  #E plane NS pol
add_cut_glyph(parent_axes=axarr[1,1],parent_fig=fig,cut='NS',pol='EW')  #H plane NS pol

#subplots_adjust(hspace=0,wspace=0)
axarr[0,0].set_ylim(-35,2)
for ax in axarr.ravel():
    ax.grid(which='both') #turn on everyones grid
    ax.xaxis.set_ticks(np.arange(-80,100,20))

if not opts.savefig is None:
    print "plotting to", opts.savefig
    fig.savefig(opts.savefig)
else:
    show()
