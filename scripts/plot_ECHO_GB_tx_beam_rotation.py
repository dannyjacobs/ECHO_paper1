"""fit for a pointing offset between the two NS runs"""
from ECHO.read_utils import read_map,write_map
from ECHO.plot_utils import rotate_hpm,get_interp_val,project_healpix
from matplotlib.pyplot import *
from scipy.optimize import fmin
import healpy as hp

#color mapping stuff
import matplotlib.colors as Colors
import matplotlib.cm as cmx
from matplotlib import colorbar as Colorbar

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
#TXmodel = np.ma.masked_where(np.sqrt(map_A_err**2+map_B_err**2)>1,TXmodel)
TXmodel = np.ma.masked_where(map_A_counts<3,TXmodel)
TXmodel = np.ma.masked_where(map_B_counts<3,TXmodel)

#normalize maps
map_B -= np.ma.mean(map_B[:3])
map_A -= np.ma.mean(map_A[:3])

theta_slice = np.linspace(-np.pi/2,np.pi/2,num=20)
phi_slice = np.zeros_like(theta_slice) + np.pi/2
theta_rots = np.array([-10])#np.arange(-20,10,5)
phi_rots = np.array([0])

#setup the line coloring
cNorm  = Colors.Normalize(vmin=np.min(theta_rots), vmax=np.max(theta_rots))
cm = get_cmap('jet')
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)

R_slice_E_err = get_interp_val(np.ma.sqrt(map_A_err**2+map_B_err**2),theta_slice,phi_slice)
print "theta,phi"
figure()
for phi_rot in phi_rots:
    for theta_rot in theta_rots:
        RXA = map_A - rotate_hpm(TXmodel,phi_rot,theta_rot,pol=pol)
        RXB = map_B - TXmodel
        TH,PH,IM = project_healpix(RXA-RXB)
        imshow(IM,vmin=-1,vmax=1)
        print theta_rot,phi_rot
        write_map('test.fits',RXA-RXB)
        show(block=True)
        #colorVal = scalarMap.to_rgba(theta_rot)
        #R_slice_E = get_interp_val(RXA-RXB,theta_slice,phi_slice)
        #errorbar(theta_slice*180/np.pi,R_slice_E,
        #             yerr=R_slice_E_err,color=colorVal)
        #plot(theta_slice*180/np.pi,R_slice_E,color=colorVal)
sys.exit()
subplots_adjust(hspace=0,wspace=0,bottom=0.33)
ax3 = axes([0.1, 0.14, 0.85, 0.05])
cb1 = Colorbar.ColorbarBase(ax3, cmap=cm,
                                norm=cNorm,
                                orientation='horizontal')
show()
