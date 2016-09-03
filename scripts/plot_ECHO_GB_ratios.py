#! /usr/bin/env python
from matplotlib.pyplot import *
from matplotlib import gridspec
import numpy as np
import sys,os,optparse
import healpy as hp
from ECHO.plot_utils import project_healpix,rotate_hpm,nf,fmt,get_interp_val,make_polycoll
from ECHO.read_utils import read_map
o = optparse.OptionParser()
o.set_description('plot the ratio of all pairs of N/S antennas')
opts,args = o.parse_args(sys.argv[1:])
style.use('./echo.mplstyle')


infiles = ['../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits']
titles = ['A NS','B NS','A EW','B EW']


#load the tx subtracted models
files = []
if False: #use the fitted transmitter map
    for i in xrange(len(infiles)):
        if infiles[i].find('Nant_NS')!=-1:
            files.append(infiles[i].replace('.fits','_txbeamfit.fits'))
        else:
            files.append(infiles[i].replace('.fits','-bicolog_legs_360.fits'))
else:
    files = [f.replace('.fits','-bicolog_legs_360_rot.fits') for f in infiles]
print "loading:"
print '\n'.join(files)
#load all the files
maps = [read_map(filename) for filename in files]
#rotate the NS maps
maps[0] = rotate_hpm(maps[0],90,0)
maps[1] = rotate_hpm(maps[1],90,0)
#make sure the maps are zero normalized
for i in xrange(len(maps)):
    maps[i] -= np.ma.mean(maps[i][:3])

fig,axarr = subplots(4,4,sharex=True,sharey=True,subplot_kw={'visible':False,'aspect':'equal'})
imaxes = []
for i in xrange(len(files)):
    for j in xrange(i,len(files)):
        print files[i],'/',files[j]
        if i==j:
            R = maps[i]
            vmin=-15
            vmax=0
        else:
            R = maps[i]-maps[j]
            vmin=-5
            vmax=5
        THETA,PHI,R_proj = project_healpix(R)
        axarr[j,i].set_visible(True)
        myax =axarr[j,i].imshow(R_proj,vmin=vmin,vmax=vmax)
        imaxes.append(myax)
        if i==j:
            axarr[j,i].set_title(titles[i])
        if i==j and j==len(files)-1:#add colorbars to the diagonals
            colorbar(imaxes[-1],ax=axarr[i,j],orientation='vertical')
        if i==0 and j==1:
            imax = imaxes[-1]
        CS = axarr[j,i].contour(THETA*180/np.pi,[20,40,60],colors='k')
        # Recast levels to new class
        CS.levels = [nf(val) for val in CS.levels]
        clabel(CS, inline=1, fontsize=10,fmt=fmt)
        axarr[j,i].set_xticklabels([])
        axarr[j,i].set_yticklabels([])
        axarr[j,i].set_ylim([-10,99])
        axarr[j,i].set_xlim([-10,99])
        axarr[j,i].xaxis.set_ticks(np.arange(-100,100,20))

subplots_adjust(bottom=0.15,wspace=0,hspace=0)
cax1=axes([.12,0.09,0.59,0.03])
colorbar(imax,cax=cax1,orientation='horizontal',label='dB')
savefig('../figures/GB_ratio_maps.png')





#make a slice plot of some of the ratios
errfiles = [f.replace('beam','rms') for f in infiles]
errmaps = [read_map(f) for f in errfiles]
fig,axarr = subplots(2,2,sharex=True,sharey=True)
ratios = [[0,1],[2,3]]
theta = np.linspace(-np.pi/2,np.pi/2,num=20)
phi = np.zeros_like(theta)
Hcolor='r'
Ecolor='b'
ls = ['-','--']
labels = ['NS','EW']
colors=['b','b']
for k,(i,j) in enumerate(ratios):
    #difference
    R = maps[i] - maps[j]
    #slices
    R_slice_E = get_interp_val(R,theta,phi)
    R_slice_H = get_interp_val(R,theta,phi+np.pi/2)
    #the error
    R_E = np.sqrt(errmaps[i]**2 + errmaps[j]**2)
    R_slice_E_err = get_interp_val(R_E,theta,phi)
    R_slice_H_err = get_interp_val(R_E,theta,phi+np.pi/2)
    #plot
    axarr[0,k].errorbar(theta*180/np.pi,R_slice_E,yerr=R_slice_E_err,label=labels[k],fmt=ls[k],color=colors[k])
    axarr[1,k].errorbar(theta*180/np.pi,R_slice_H,yerr=R_slice_H_err,fmt=ls[k],color=colors[k])
#get the Neben et al ratio results
#R_OC_NS = read_map('../data/null4_north_over_south_ratio_ns.fits') #the old data from aug 2015
#R_OC_NS = rotate_hpm(R_OC_NS,angle=90)
R_OC_NS = read_map('../data/null4_north_over_south_ratio_ns_nside64.fits')
#R_OC_EW = read_map('../data/null4_north_over_south_ratio_ew.fits')
R_OC_EW = read_map('../data/null4_north_over_south_ratio_ew_nside64.fits')

R_slice_E_NS = get_interp_val(R_OC_NS,theta,phi)
R_slice_H_NS = get_interp_val(R_OC_NS,theta,phi+np.pi/2)

axarr[0,0].plot(theta*180/np.pi,R_slice_E_NS,ls=':',label='ORBCOMM NS',color='k')
axarr[1,0].plot(theta*180/np.pi,R_slice_H_NS,ls=':',color='k')

#try out some ways of plotting estimated orbcomm uncertainty
OC_err = 0.5
axarr[0,0].fill_between(theta*180/np.pi,R_slice_E_NS-OC_err,y2=R_slice_E_NS+OC_err,alpha=0.2,color='k')
axarr[1,0].fill_between(theta*180/np.pi,R_slice_H_NS-OC_err,y2=R_slice_H_NS+OC_err,alpha=0.2,color='k')

R_slice_E_EW = get_interp_val(R_OC_EW,theta,phi)
R_slice_H_EW = get_interp_val(R_OC_EW,theta,phi+np.pi/2)
axarr[0,1].plot(theta*180/np.pi,R_slice_E_EW,'.',label='ORBCOMM EW',color='k')
axarr[1,1].plot(theta*180/np.pi,R_slice_H_EW,'.',color='k')

axarr[0,1].fill_between(theta*180/np.pi,R_slice_E_EW-OC_err,y2=R_slice_E_EW+OC_err,alpha=0.2,color='k')
axarr[1,1].fill_between(theta*180/np.pi,R_slice_H_EW-OC_err,y2=R_slice_H_EW+OC_err,alpha=0.2,color='k')


#axarr[0].legend(loc='best',fontsize=9)
axarr[0,0].set_ylabel('E plane\n [dB V/m]')
axarr[1,0].set_ylabel('H plane\n [dB V/m]')
axarr[0,0].set_title('NS')
axarr[0,1].set_title('EW')
axarr[1,0].set_xlabel('$\\theta$')
axarr[1,1].set_xlabel('$\\theta$')
for ax in axarr.ravel():
    ax.grid(which='both') #turn on everyones grid
    ax.xaxis.set_ticks(np.arange(-80,100,20))
subplots_adjust(hspace=0,wspace=0)
savefig('../figures/GB_ratio_slices.png')
axarr[0,0].set_ylim(-5,5)
savefig('../figures/GB_ratio_slice_5dB.png')






#plot a comparison of the ratio maps
fig,axarr = subplots(2,2,sharex=True,sharey=True)
mn,mx = -1,1
#         NS - EW
#Orbcomm
#ECHO
nside=hp.npix2nside(len(R_OC_NS))
coll = make_polycoll(R_OC_NS,nsides=nside,cmap=cm.jet)
coll.set_clim(mn,mx)
axarr[0,0].add_collection(coll)
axarr[0,0].autoscale_view()

nside=hp.npix2nside(len(R_OC_EW))
coll = make_polycoll(R_OC_EW,nsides=nside,cmap=cm.jet)
coll.set_clim(mn,mx)
axarr[0,1].add_collection(coll)
axarr[0,1].autoscale_view()

#maps = [read_map(filename) for filename in files]

nside=hp.npix2nside(len(maps[0]))
coll = make_polycoll(maps[0]-maps[1],nsides=nside,cmap=cm.jet)
coll.set_clim(mn,mx)
axarr[1,0].add_collection(coll)
axarr[1,0].autoscale_view()


coll = make_polycoll(maps[2]-maps[3],nsides=nside,cmap=cm.jet)
coll.set_clim(mn,mx)
axarr[1,1].add_collection(coll)
axarr[1,1].autoscale_view()

subplots_adjust(bottom=0.15,wspace=0,hspace=0)
cax1=axes([.1,0.09,0.9,0.03])
colorbar(coll,cax=cax1,orientation='horizontal',label='dB')
show()
