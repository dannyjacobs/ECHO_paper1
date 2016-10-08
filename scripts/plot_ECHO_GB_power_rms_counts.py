from ECHO.plot_utils import make_polycoll,project_healpix,nf,fmt,cmap_discretize
import numpy as np,optparse
import healpy as hp
from matplotlib.pyplot import *
from matplotlib import ticker,cm

o = optparse.OptionParser()
o.set_description('make three panel plot of orbcomm beams')
o.add_option('--max',type=float,help='max plotted value')
o.add_option('--min',type=float,help='min plotted value')
o.add_option('--savefig',default=None,type=str)
opts,args = o.parse_args(sys.argv[1:])


#load the data
beamfile='../data/8_power_Nant_NSdipole_NStransmitter.fits'
rmsfile = beamfile.replace('power','rms')
countsfile = beamfile.replace('power','counts')
beam = np.ma.masked_invalid(hp.read_map(beamfile))
rms = np.ma.masked_invalid(hp.read_map(rmsfile))
counts = np.ma.masked_invalid(hp.read_map(countsfile))

THETA,PHI,IM = project_healpix(beam)
X,Y = np.meshgrid(
        np.linspace(-1,1,num=THETA.shape[0]),
        np.linspace(-1,1,num=THETA.shape[1]))

figure(figsize=(10,4))


#make the beam subplot
ax1 = subplot(131)
axis('equal')
beamcoll = make_polycoll(beam,cmap=cm.jet)
beamcoll.set_clim(-40,0)
ax1.add_collection(beamcoll)
CS = ax1.contour(X,Y,THETA*180/np.pi,[20,40,60],colors='k')
CS.levels = [nf(val) for val in CS.levels]
clabel(CS, inline=1, fontsize=10,fmt=fmt)
ax1.autoscale_view()
ax1.set_yticklabels([])
ax1.set_xticklabels([])
ax1.set_title('Gridded Power')
cb = colorbar(beamcoll, ax=ax1,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('dB')


#make the std plot
ax2 = subplot(132)
axis('equal')
rmscoll = make_polycoll(rms,cmap=cm.jet)
rmscoll.set_clim(0,2)
ax2.add_collection(rmscoll)
ax2.autoscale_view()
ax2.set_title('Standard Deviation')
ax2.set_yticklabels([])
ax2.set_xticklabels([])
CS = contour(X,Y,THETA*180/np.pi,[20,40,60],colors='w')
clabel(CS, inline=1, fontsize=10,fmt=fmt)
cb = colorbar(rmscoll, ax=ax2,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('dB')

#make the count plot
ax3 = subplot(133)
axis('equal')
#quantize the counts into log bins
# bins = np.digitize(counts,[0,1,10,100])
# counts[bins] = np.array([0,1,10,100])
print counts.min(),counts.max()
cmap = cmap_discretize(cm.bone_r,10)
#counts = np.log(counts)
countscoll = make_polycoll(counts,cmap=cmap)
countscoll.set_clim(1,10)
ax3.add_collection(countscoll)
ax3.autoscale_view()
ax3.set_title('Sample Count')
ax3.set_yticklabels([])
ax3.set_xticklabels([])
CS = contour(X,Y,THETA*180/np.pi,[20,40,60],colors='w')
clabel(CS, inline=1, fontsize=10,fmt=fmt)
cb = colorbar(countscoll, ax=ax3,orientation='horizontal')
# tick_locator = ticker.MaxNLocator(nbins=5)
# cb.locator = tick_locator
# cb.update_ticks()
# cb.set_ticks(np.log2(np.array([2,4,8,16])))
# cb.set_ticklabels([2,4,8,16])
# cb.set_ticks(np.arange(2**counts.min(),16,
#         np.ceil((counts.max()-counts.min())/6)))
cb.set_label('sample count')
if not opts.savefig is None:
    print "saving to",opts.savefig
    savefig(opts.savefig)
