from matplotlib.pyplot import *
import sys
from glob import glob
import numpy as np
from astropy.time import Time
from ECHO.read_utils import read_apm_logs,flag_angles,apply_flagtimes,read_orbcomm_spectrum,channel_select,interp_rx,dB
import matplotlib.dates as mdates

#find all our files
rx_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/satpowerflight12.0*')
assert(len(rx_files)>0)
apm_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/apm_Aug13_v*log')



#get data apm data
positiontimes,positions,angletimes,angles,cmdtimes,cmds = read_apm_logs(apm_files)
print 'gps start-end',positiontimes[0].iso,positiontimes[-1].iso
print 'att start-end',angletimes[0].iso,angletimes[-1].iso

#get orbcomm data
rxtimes,freqs,rxspectrum = read_orbcomm_spectrum(rx_files,'S','NS')
#get the power in our current channel
rx_power = channel_select(freqs,rxspectrum,137.5)
#interpolate the rx data down to match the gps times
rx_interp = interp_rx(positiontimes,rxtimes,rx_power)

fig = figure(figsize=(15,10))
title('South dipoole, NS pol')
#flagging the yaws
yawmask,badyawtimes = flag_angles(angletimes,angles,1)
print "found {n} bad yaws".format(n=len(badyawtimes))
yawmask = np.reshape(yawmask,(1,-1))
angles = np.ma.masked_where(yawmask,angles)
plot_date(angletimes.plot_date,angles[0],'.k') #plot the flagged yaws
ylabel('heading [deg]')

#plot the data
twinx()
#applying the yaw flags to the data
posmask = apply_flagtimes(positiontimes,badyawtimes,1.0)
rx_interp = np.ma.masked_where(posmask,rx_interp)
plot_date(positiontimes.plot_date,dB(rx_interp),'-k') #plot the flagged altitudes
ylabel('power [dB]')

#set some nice time formatting
hours = mdates.HourLocator()
gca().xaxis.set_major_locator(hours)
hourFmt = mdates.DateFormatter('%H:%M')
gca().xaxis.set_major_formatter(hourFmt)
minutes = mdates.MinuteLocator(interval=5)
gca().xaxis.set_minor_locator(minutes)
minFmt = mdates.DateFormatter('%M')
gca().xaxis.set_minor_formatter(minFmt)
fig.autofmt_xdate()
grid()
xlabel('time [minutes]')
show()
