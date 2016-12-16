from matplotlib.pyplot import *
import sys
from glob import glob
import numpy as np
from astropy.time import Time
from ECHO.read_utils import read_apm_logs,flag_angles,apply_flagtimes,read_orbcomm_spectrum,channel_select,interp_rx,dB
import matplotlib.dates as mdates
style.use('./echo.mplstyle')
#find all our files
rx_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/satpowerflight12.0*')
assert(len(rx_files)>0)
apm_files = glob('/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/apm_Aug13_v*log')



#get data apm data
positiontimes,positions,angletimes,angles,cmdtimes,cmds = read_apm_logs(apm_files)
#get orbcomm data
rxtimes,freqs,rxspectrum = read_orbcomm_spectrum(rx_files,'S','NS')
#get the power in our current channel
rx_power = channel_select(freqs,rxspectrum,137.5)
#interpolate the rx data down to match the gps times
rx_interp = interp_rx(positiontimes,rxtimes,rx_power)


fig = figure(figsize=(15,10))
title('Antenna A, NS pol')


# plot_date(cmdtimes.plot_date,np.ones(len(cmdtimes)),'.k') #plot the flagged yaws

# ylabel('waypoint')


#flagging the yaws
yawmask,badyawtimes = flag_angles(angletimes,angles,2)
print "found {n} bad yaws".format(n=len(badyawtimes))
#plot the data
#twinx()
#applying the yaw flags to the data
posmask = apply_flagtimes(positiontimes,badyawtimes,1.0)
rx_interp = np.ma.masked_where(posmask,rx_interp)
print "total flags after yaw flagging:",np.sum(rx_interp.mask)
#flagging the waypoints
cmdmask = apply_flagtimes(positiontimes,cmdtimes,0.5)
plot_date(positiontimes.plot_date,rx_interp,'.r',ms=10)
rx_interp = np.ma.masked_where(cmdmask,rx_interp)
print "total flags after cmd flagging:",np.sum(rx_interp.mask)

plot_date(positiontimes.plot_date,rx_interp,'.k',ms=10) #plot the flagged altitudes
ylabel('power [dB]')

#set some nice time formatting
seconds = mdates.SecondLocator(interval=10)
gca().xaxis.set_minor_locator(seconds)
secondFmt = mdates.DateFormatter('%S')
gca().xaxis.set_minor_formatter(secondFmt)

minutes = mdates.MinuteLocator(interval=5)
gca().xaxis.set_major_locator(minutes)
minFmt = mdates.DateFormatter('%M:00')
gca().xaxis.set_major_formatter(minFmt)
fig.autofmt_xdate()

for x in cmdtimes.plot_date:
    axvline(x,color='k')

grid()
if True:
    twinx()
    yaws = angles.squeeze()-np.mean(angles)
    plot_date(angletimes.plot_date,yaws)
    ylim(-5,5)

xlabel('time [minutes]')
tight_layout()
show()
