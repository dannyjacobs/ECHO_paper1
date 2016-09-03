from mwapy.pb import primary_beam
import optparse,sys
import healpy as hp
import numpy as np


def dB(x):
    #for converting power to dB
    return 10*np.log10(x)
def dB20(x):
    #for converting V/m to dB
    return 20*np.log10(x)

o = optparse.OptionParser()
o.set_description('read in the latest mwa beam model and output as a healpix file')
o.add_option('--nside',type=int,
    help='nside of the output map')
o.add_option('--freq',type=float,
    help='desired frequency in MHz')
o.add_option('--delays',type=str,default="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0",
    help='16 comma delimited integer delays, default is all 0s')
o.add_option('--outname',default='mwatile_beam',
    help='output file prefix')
opts,args = o.parse_args(sys.argv[1:])

delays = map(int,opts.delays.split(','))
pixes = np.arange(hp.nside2npix(opts.nside))
theta,phi = hp.pix2ang(opts.nside,pixes)
#the mwa primary beam function expects a meshgrid,
#   so we need to make our pixels 2d
theta = np.reshape(theta,(1,-1))
phi = np.reshape(phi,(1,-1))

#data file copied from MWA_Tools/pb/ v. mwapy-v1.1.0_103_g7763af8
# on August 7, 2016
h5filepath='../data/mwa_full_embedded_element_pattern.h5'
rX,rY=primary_beam.MWA_Tile_full_EE(h5filepath, theta, phi,
                                     freq=opts.freq*1e6, delays=delays,
                                     zenithnorm=True, power=True )
rX = dB(rX.squeeze()) #convert to power dB
print "subtracting ",rX[0]
rX -= rX[0] #normalize to peak
rY = dB(rY.squeeze())
rY -= rY[0]

hp.write_map(opts.outname+'_X.fits',rX)
hp.write_map(opts.outname+'_Y.fits',rY)
