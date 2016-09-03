from matplotlib.pyplot import *
import sys,glob
import numpy as np
from astropy.time import Time
from ECHO.read_utils import read_apm_logs,flag_angles,apply_flagtimes
import matplotlib.dates as mdates

infile = sys.argv[1]
outfile = sys.argv[2]

apm_files = glob.glob(infile)
assert(len(apm_files)>0)

#get data
postimes,positions,angletimes,angles,cmdtimes,cmds = read_apm_logs(apm_files)
print 'gps start-end',postimes[0].iso,postimes[-1].iso
print 'att start-end',angletimes[0].iso,angletimes[-1].iso



fig = figure(figsize=(15,10))
title(infile)
plot_date(angletimes.plot_date,angles[0],',k')
#flagging the yaws
yawmask,badyawtimes = flag_angles(angletimes,angles,2)
print "found {n} bad yaws".format(n=len(badyawtimes))
yawmask = np.reshape(yawmask,(1,-1))
angles = np.ma.masked_where(yawmask,angles)
plot_date(angletimes.plot_date,angles[0],'ok') #plot the flagged yaws
ylabel('heading [deg]')

#plot the altitude
twinx()
plot_date(postimes.plot_date,positions[2],',k')
#applying the yaw flags to the positions
posmask = apply_flagtimes(postimes,badyawtimes,1.0)
#extend the time mask to the data
posmask = np.reshape(posmask,(1,-1))
posmask = np.repeat(posmask,3,axis=0)
positions = np.ma.masked_where(posmask,positions)
plot_date(postimes.plot_date,positions[2],'.k') #plot the flagged altitudes
ylabel('alt [m]')

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
savefig(outfile)




'''
spec_times = []
spec_raw = []
spec_files = glob.glob(specfile)
for spec_file in spec_files:
    print 'Reading in %s...' %spec_file
    lines = open(spec_file).readlines()
    if not len(lines) == 0:
        for line in lines[1:]:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n').split(',')
            if len(line) == (1025): # Make sure line has finished printing
                spec_times.append(float(line[0]))
                spec_raw.append(map(float,line[1:]))
spec_times = Time(spec_times,format='unix')
spec_raw = np.array(spec_raw)
print spec_raw.shape
freq_chan = np.where(spec_raw[0,:]==137.500)

spec = spec_raw[1:,freq_chan]
spec_times = Time(spec_raw[:,0],format='unix')


ax1 = fig.add_subplot(211)
ax1.plot(spec_times.gps,spec_raw[:,freq_chan])
ax1.set_ylabel('dB')

fig = figure(figsize=(16,9))
ax2 = fig.add_suplot(212)
ax2.plot(ATT_times.gps,yaws)
ax2.set_ylabel('deg')
'''
