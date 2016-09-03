from ECHO.read_utils import read_map,write_map
from ECHO.plot_utils import rotate_hpm
import numpy as np

def avg_beams(beams,errors):
    """
    input:
        beams = list of beams (numpy arrays)
        errors = list of matching error maps
    return averaged_beam,average_beam_error
    """
    beams = np.ma.array(beams)
    beams -= np.reshape(np.ma.mean(beams[:,:5],axis=1),(-1,1))
    errors = np.ma.array(errors)
    #inverse variance weighting
    avg_beam = np.ma.mean(beams/errors**2,axis=0)/np.ma.mean(1/errors**2,axis=0)

    #avg_beam_err = 1/np.sqrt(np.ma.mean(1./errors**2,axis=0))
    avg_beam_std = np.ma.mean((beams/errors**2)**2,axis=0)/np.ma.mean(1/errors**4,axis=0)
    avg_beam_std -= avg_beam**2
    avg_beam_std = np.sqrt(avg_beam_std)

    return avg_beam,avg_beam_std
def read_maps(filenames):
    """
    input:
        list of beam filenames
    return:
        list of maps
    """
    beams = []
    errors = []
    for filename in filenames:
        print "loading:",filename
        bm = read_map(filename)
        beams.append(bm)
    return beams
NS_filenames  = ['../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits',
'../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits']
EW_filenames = ['../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits',
'../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits']
#combine the matching polarizations between the two antennas
pol = ['NS','EW']
for i,infiles in enumerate([NS_filenames,EW_filenames]):
    rmsfiles = [f.replace('beam','rms') for f in infiles]
    filenames = [f.replace('.fits','-bicolog_legs_360.fits') for f in infiles]
    beams = read_maps(filenames)
    rmsbeams = read_maps(rmsfiles)
    avg_beam,avg_beam_std = avg_beams(beams,rmsbeams)

    write_map('../data/GB_rx_model_beam_{pol}.fits'.format(pol=pol[i]),avg_beam)
    write_map('../data/GB_rx_model_rms_{pol}.fits'.format(pol=pol[i]),avg_beam_std)


#make a map of the full combined dataset
infiles = NS_filenames+EW_filenames
rmsfiles = [f.replace('beam','rms') for f in infiles]
filenames = [f.replace('.fits','-bicolog_legs_360.fits') for f in infiles]
beams = read_maps(filenames)
rmsbeams = read_maps(rmsfiles)

#rotate NS beams to match EW beams
for i,filename in enumerate(infiles):
    if filename.find('NS')!=-1:
        beams[i] = rotate_hpm(beams[i],angle=90)
        rmsbeams[i] = rotate_hpm(rmsbeams[i],angle=90)
avg_beam,avg_beam_std = avg_beams(beams,rmsbeams)
write_map('../data/GB_rx_model_beam.fits',avg_beam)
write_map('../data/GB_rx_model_rms.fits',avg_beam_std)
