"""fit for a pointing offset between the two NS runs"""
from ECHO.read_utils import read_map,write_map
from ECHO.plot_utils import rotate_hpm,get_interp_val,project_healpix
from matplotlib.pyplot import *
from scipy.optimize import fmin
import healpy as hp

A_file = '../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits'
map_A   = read_map(A_file)
map_A_err = read_map(A_file.replace('beam','rms'))
map_A_counts = read_map(A_file.replace('beam','counts'))

B_file = '../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits'
map_B   = read_map(B_file)
map_B_err   = read_map(B_file.replace('beam','rms'))
map_B_counts = read_map(B_file.replace('beam','counts'))
TXmodel = read_map('../data/bicolog_legs_360.fits')
pol = 'NS'
TXmodel = rotate_hpm(TXmodel,90,0,pol=pol)

#heavily flag the maps (only need to flag the model, the flags will propagate when we subtract)
TXmodel = np.ma.masked_where(np.sqrt(map_A_err**2+map_B_err**2)>1,TXmodel)
TXmodel = np.ma.masked_where(map_A_counts<3,TXmodel)
TXmodel = np.ma.masked_where(map_B_counts<3,TXmodel)


def chisq(angles,map_A,map_A_err,map_B,map_B_err,TXmodel,pol):
    theta = angles[0]
    phi = angles[1]
    #phi=0
    RXA = map_A - rotate_hpm(TXmodel,phi,theta,pol)
    RXB = map_B  - TXmodel
    return np.ma.mean((RXA-RXB)**2/(2*(map_A_err**2+map_B_err**2)))
#Calibrate the transmitter angle of map A against map B
print "fitting for a transmitter rotation"
result = fmin(chisq,[0,0],args=(map_A,map_A_err,map_B,map_B_err,TXmodel,pol))
print "theta (deg) = ",np.round(result[0],2)
print "phi (deg) = ",np.round(result[1],2)
TXmodel.mask = False

#make a tx subtracted
#choosing -5 from the E slice
RXA = map_A - rotate_hpm(TXmodel,0,-10,pol)
#Set the top to 0
RXA -= np.ma.mean(RXA[:5])

# RXB = map_B - TXmodel
# RXB -= np.ma.mean(RXB[:5])
#
# TH,PH,ANTflat = project_healpix(RXA-RXB)
# imshow(ANTflat,vmin=-5,vmax=5)


A_file_corrected = A_file.replace('.fits','_txbeamfit.fits')
print "outputing beam subtracted map:", A_file_corrected
write_map(A_file_corrected,RXA)



sys.exit()
#test code kept for posterity
TXmodel = np.ma.masked_where(np.sqrt(map_A_err**2+map_B_err**2)>1,TXmodel)
TXmodel = np.ma.masked_where(map_A_counts<3,TXmodel)
TXmodel = np.ma.masked_where(map_B_counts<3,TXmodel)
if False:
    for phi in np.arange(80,110,1):
        RXA = map_A - rotate_hpm(TXmodel,phi)
        TH,PH,ANTflat = project_healpix(RXA-RXB)
        print phi,np.std(RXA-RXB)
        show(block=True)
if False:
    thetas = np.arange(-45,45,5)
    theta_slice = np.linspace(-np.pi/2,np.pi/2,num=20)
    phi_slice = np.zeros_like(theta_slice)
    figure()
    phi=1.15
    CHISQ = []
    RXB = map_B - TXmodel
    for theta in thetas:
        RXA = map_A - rotate_hpm(TXmodel,phi,theta_angle=theta)
        R_slice_E = get_interp_val(RXA-RXB,theta_slice,phi_slice + np.pi/2)
        plot(theta_slice,R_slice_E)
        #TH,PH,ANTflat = project_healpix(RXA-RXB)
        #imshow(ANTflat,vmin=-5,vmax=5)
        #CHISQ.append(chisq([theta,phi],map_A,map_A_err,map_B,map_B_err,TXmodel))
        print theta, chisq([theta,phi],map_A,map_A_err,map_B,map_B_err,TXmodel)
        show(block=True)
    #plot(thetas,CHISQ)
if False:
    TXmodel = rotate_hpm(TXmodel,90)
    thetas = np.arange(-45,45,0.25)
    phis = np.arange(80,110,5)

    CHISQ = []
    for phi in phis:
        for th in thetas:
            CHISQ.append(chisq([th,phi],map_A,map_A_err,map_B,map_B_err,TXmodel))
        plot(thetas,[chisq([th,phi],map_A,map_A_err,map_B,map_B_err,TXmodel) for th in thetas])
    THETA,PHI = np.meshgrid(thetas,phis)
    CHISQ = np.reshape(np.array(CHISQ),THETA.shape)
    print THETA.ravel()[CHISQ.argmin()],PHI.ravel()[CHISQ.argmin()]
    show()
    sys.exit()
