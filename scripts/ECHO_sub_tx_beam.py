#!//anaconda/bin/python

import healpy as hp
import numpy as np
from matplotlib.pyplot import *
import sys,os,optparse
from ECHO.read_utils import read_map,write_map
from ECHO.plot_utils import project_healpix,rotate_hpm

o = optparse.OptionParser()
o.set_description('ECHO_tx_beam_cal.py ECHOMap.fits receiver_ant_beam_model.fits ')
o.add_option('--theta',type=float,help='theta rotation (applied after phi, deg)')
o.add_option('--phi',type=float,help='phi rotation (deg)')
opts,args = o.parse_args(sys.argv[1:])



ECHOmapFile = args[0]
TXmodelFile = args[1]



ECHOmap = read_map(ECHOmapFile)
ECHOmapnside = hp.get_nside(ECHOmap.data)



#load the model
TXmodel = read_map(TXmodelFile)
#rotate the data to match the model
if ECHOmapFile.find('NS')!=-1:
    pol = 'NS'
    print "rotating TXmodel to match NS beam"
    TXmodel = rotate_hpm(TXmodel,90,0,pol=pol)
else:
    pol = 'EW'
print("assuming beams stored in dB, and assuming maps stored in RING")

print("applying rotation to TX model theta={theta}, phi={phi}".format(theta=opts.theta,phi=opts.phi))
TXmodel = rotate_hpm(TXmodel,opts.phi,opts.theta,pol=pol)

#power received = power transmitted
#                +receiver beam gain (which is negative)
#                +transmitter beam gain (also negative)
#assume receiver beam gain = transmitter beam gain = 0 at zenith
if True:
    print("normalizng the data map to the mean of the top few healpix pixels ")
    print ECHOmap[:5]==hp.UNSEEN,ECHOmap[:5]
    print "subtracting",np.ma.average(ECHOmap[:5])
    echonorm = np.mean(ECHOmap[:5])
    ECHOmap -= echonorm
TXmodel -= np.mean(TXmodel[:5])


#ECHOmap = receiver beam gain (RX) + transmitter beam gain (TX)
#solving for the ECHObeam map we get
ECHObeam = ECHOmap - TXmodel

print "saving output difference",
outfile=ECHOmapFile[:-5]+'-'+os.path.basename(TXmodelFile[:-5])+'_rot.fits'
print outfile
write_map(outfile,ECHObeam)

THETA,PHI,ECHOmapflat = project_healpix(ECHOmap)
THETA,PHI,TXmapflat = project_healpix(TXmodel)
THETA,PHI,ECHObeamflat = project_healpix(ECHObeam)

figure(figsize=(15,6))
suptitle(outfile[:-5])
subplot(131)
imshow(ECHOmapflat,vmin=-20,vmax=0)
title("Measured Power")
subplot(132)
imax = imshow(TXmapflat,vmin=-20,vmax=0)
title("Reciever Antenna Model")

subplot(133)
diffax = imshow(ECHObeamflat,vmin=-20,vmax=0)
title('Difference')
tight_layout()
subplots_adjust(bottom=0.2,top=0.9)
cax1=axes([.07,0.09,0.55,0.05])
colorbar(imax,cax=cax1,orientation='horizontal',label='dB')
cax=axes([.7,0.09,0.29,0.05])
colorbar(diffax,cax=cax,orientation='horizontal',label='dB')
draw()
savefig('../figures/'+os.path.basename(outfile)[:-5]+'.png')
print "done."
