from matplotlib.pyplot import *
import sys,glob
import numpy as np
from astropy.time import Time

infile = sys.argv[1]
specfile = sys.argv[2]

SEC_PER_WEEK = 604800

lats,lons,alts = [],[],[]
ATT_times,yaws = [],[]
weektimes = []
isGPS = True
apm_files = glob.glob(infile)
for apm_file in apm_files:
   print 'Reading in %s...' %apm_file
   lines = open(apm_file).readlines()
   if not len(lines) == 0:
      for line in lines:
         if line.startswith('GPS'):
            if isGPS:
                startTime = float(line.split(',')[-1])
                isGPS = False
            lats.append(map(float,line.split(',')[7:8]))
            lons.append(map(float,line.split(',')[8:9]))
            alts.append(map(float,line.split(',')[9:10]))
            weektimes.append(map(float,line.split(',')[3:5])) #ms and week number
         if line.startswith('ATT'):
            ATT_times.append(map(float,[line.split(',')[1].strip(' ')]))
            yaws.append(map(float,[line.split(',')[7].strip(' ')]))
weektimes = np.array(weektimes)
apm_times = weektimes[:,1]*SEC_PER_WEEK+weektimes[:,0]/1000.
apm_times = Time(apm_times,format='gps')
ATT_times = np.array(ATT_times)
ATTGPSseconds = ATT_times/1000.+apm_times.gps[0]-startTime/1000.
ATT_times = Time(ATTGPSseconds, format = 'gps')
lats = np.array(lats).squeeze()
lons = np.array(lons).squeeze()
alts = np.array(alts).squeeze()
yaws = np.array(yaws).squeeze()

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

plot(ATT_times.gps,yaws)
show()
