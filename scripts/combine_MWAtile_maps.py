from ECHO.read_utils import read_map,write_map
import sys,numpy as np


files = ['../data/acc_MWAtile_2016_run2_EWtx_EWrx_8_beam.fits',
        '../data/acc_MWAtile_2016_run1_EWtx_EWrx_8_beam.fits']
rmsfiles = [f.replace('beam','rms') for f in files]
countfiles = [f.replace('beam','counts') for f in files]
offsets = [2,0] #use these to level the different attenuators used on the runs

beams,errs,counts = [],[],[]
for i,filename in enumerate(files):
    beam= read_map(filename)
    err = read_map(rmsfiles[i])
    count = read_map(countfiles[i])
    beam += offsets[i]
    beams.append(beam)
    errs.append(err)
    counts.append(count)
beams = np.ma.array(beams)
errs = np.ma.array(errs)
counts = np.ma.array(counts)
#average things together,
#   acounting for healpix cells that have been measured in both maps
sumcounts = np.ma.sum(counts,axis=0)
meanbeam = np.ma.sum(beams*counts,axis=0)/sumcounts
meanerr = np.sqrt(np.ma.sum((errs**2 + beams**2)*counts,axis=0)/sumcounts - meanbeam**2)

#normalizations going to be shot by the scalings, lets renormalize
print "subtracting: ",np.ma.mean(meanbeam[:5])
meanbeam -= np.ma.mean(meanbeam[:5])

#write it all out
outbeam  = '../data/acc_MWAtile_2016_EWtx_EWrx_8_beam.fits'
write_map(outbeam,meanbeam)
write_map(outbeam.replace('beam','rms'),meanerr)
write_map(outbeam.replace('beam','counts'),sumcounts)
