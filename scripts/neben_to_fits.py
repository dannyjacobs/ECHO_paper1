"""Load the raw data sent by abraham and save to fits files."""
import numpy as np
import healpy as hp
from ECHO.plot_utils import rotate_hpm,get_interp_val
from ECHO.read_utils import write_map


def mirror_hpm(hpm):
    """Mirror the map about Y (phi=0) axis"""
    nside = hp.npix2nside(len(hpm))
    indexes = np.arange(len(hpm))
    theta,phi = hp.pix2ang(nside,indexes)      #get the input map angles
    #interpolate to a higher nside, rotate that thing
    interp_map = hp.ud_grade(hpm,nside*2)
    theta_interp,phi_interp = hp.pix2ang(nside*2,np.arange(len(hpm)*4))
    phi_interp *= -1  #apply the azimuthal rotation
    #get the rotated interpolated map
    rotated_interpolated_beam = get_interp_val(hpm,theta_interp,phi_interp)
    #interpolate back to the original map pixels
    rotated_hpm = get_interp_val(rotated_interpolated_beam,theta,phi)
    return rotated_hpm


for n in [1,4]:
    nullfile = '../data/null{n}_north_over_south_ratio_ew_ns_nside64.csv'.format(n=n)
    nulldata = np.loadtxt(nullfile,delimiter=',')
    counts = nulldata[:,2]
    NSbeam = np.ma.masked_where(counts==0,nulldata[:,0])
    EWbeam = np.ma.masked_where(counts==0,nulldata[:,1])
    NSbeam = rotate_hpm(mirror_hpm(NSbeam),-90)
    EWbeam = rotate_hpm(mirror_hpm(EWbeam),-90)
    outNS = "../data/null{n}_north_over_south_ratio_ns_nside64.fits".format(n=n)
    outEW = "../data/null{n}_north_over_south_ratio_ew_nside64.fits".format(n=n)
    write_map(outNS,NSbeam)
    write_map(outEW,EWbeam)
