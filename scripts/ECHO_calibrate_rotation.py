#! /usr/bin/env python
from matplotlib.pyplot import *
import numpy as np
import sys,os,optparse
import healpy as hp
from ECHO.read_utils import read_map
from ECHO.plot_utils import project_healpix,rotate_hpm
from scipy.optimize import minimize
o = optparse.OptionParser()
o.set_description('find a common beam model')
opts,args = o.parse_args(sys.argv[1:])


files = ['../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits']
pols = [0,0,90,90]
def mfunc(angle,mapA,mapB,mapA_err,mapB_err):
    #function to minimize
    THETA,PHI,IMA = project_healpix(mapA)
    _,_,IMA_err = project_healpix(mapA_err)
    _,_,IMB = project_healpix(mapB,angle)
    _,_,IMB_err = project_healpix(mapB_err,angle)
    return np.sqrt(np.mean(np.abs(IMA-IMB)**2/(IMA_err**2 + IMB_err**2)))

rots = np.zeros(shape=(len(files),len(files)))
for i in xrange(len(files)):
    for j in xrange(i,len(files)):
        if i==j:continue
        mapA = read_map(files[i])
        mapA_err = read_map(files[i].replace('beam','rms'))
        mapB = read_map(files[j])
        mapB_err = read_map(files[j].replace('beam','rms'))
        #remove polarization angles
        mapA = rotate_hpm(mapA,pols[i])
        mapB = rotate_hpm(mapB,pols[j])
        res = minimize(mfunc,np.array([0]),args=(mapA,mapB,mapA_err,mapB_err))
        rots[i,j] = res.x
print rots
