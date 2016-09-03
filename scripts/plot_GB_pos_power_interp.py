import numpy as np
import sys,optparse
from matplotlib.pyplot import *
from ECHO.read_utils import read_orbcomm_spectrum,read_apm_logs,interp_rx,channel_select
from glob import glob
rx_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/satpowerflight12.0*')
assert(len(rx_files)>0)
apm_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/apm_Aug13_v*log')
positiontimes,positions,angletimes,angles = read_apm_logs(apm_files)

colors=['c','b','k','m']
i=0
for ant in ['N','S']:
    for pol in ['NS','EW']:
        rxtimes,freqs,rxspectrum = read_orbcomm_spectrum(rx_files,ant,pol)
        #get the power in our current channel
        rx_power = channel_select(freqs,rxspectrum,137.5)
        #interpolate the rx data down to match the gps times
        rx_interp = interp_rx(positiontimes,rxtimes,rx_power)
        plot_date(rxtimes.plot_date,rx_power,'.',color=colors[i])
        plot_date(positiontimes.plot_date,rx_interp,'-',color=colors[i])
        i+=1
gcf().autofmt_xdate()
show()
