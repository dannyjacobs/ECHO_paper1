import numpy as np
import glob
from astropy.time import Time

from time_utils import flight_time_filter,waypt_time_filter

SEC_PER_WEEK = 604800
APMLOG_SEC_PER_TICK = 1.0e-6
NFFT = 1024


def get_data(infile,filetype=None,freqs=[],freq=0.0,freq_chan=None,
                   ant=None,dip=None,width=100,times=None,waypts=None):#isList=False,

    if filetype == 'gps':
        gps_arr = []
        gps_files = glob.glob(infile)
        for gps_file in gps_files:
            lines = open(gps_file).readlines()
            if not len(lines) == 0:
                for line in lines:
                    if line.startswith('#'):
                        continue
                    line = line.rstrip('\n').split(',')
                    if len(line) == 4: # Make sure line has finished printing
                        gps_arr.append(map(float,line))
        gps_arr = np.array(gps_arr)
        gps_times,lats,lons,alts = np.split(gps_arr,
                                            [1,2,3],
                                            axis=1)
        gps_times = Time(gps_times,format='gps')
        lats = lats.squeeze()
        lons = lons.squeeze()
        alts = alts.squeeze()
        return gps_times,lats.squeeze(),lons.squeeze(),alts.squeeze()

    elif filetype == 'apm':
        lats,lons,alts = [],[],[]
        weektimes = []
        apm_files = glob.glob(infile)
        for apm_file in apm_files:
            print 'Reading in %s...' %apm_file
            lines = open(apm_file).readlines()
            if not len(lines) == 0:
                for line in lines:
                    if line.startswith('GPS'):
                        lats.append(map(float,line.split(',')[7:8]))
                        lons.append(map(float,line.split(',')[8:9]))
                        alts.append(map(float,line.split(',')[9:10]))
                        weektimes.append(map(float,line.split(',')[3:5])) #ms and week number
        weektimes = np.array(weektimes)
        apm_times = weektimes[:,1]*SEC_PER_WEEK+weektimes[:,0]/1000.
        apm_times = Time(apm_times,format='gps')
        lats = np.array(lats).squeeze()
        lons = np.array(lons).squeeze()
        alts = np.array(alts).squeeze()
        return apm_times,lats,lons,alts

    elif filetype == 'sh':
        spec_times = []
        spec_raw = []
        spec_files = glob.glob(infile)
        for spec_file in spec_files:
            print 'Reading in %s...' %spec_file
            lines = open(spec_file).readlines()
            if not len(lines) == 0:
                if len(freqs) == 0:
                    freqs = np.array(map(float,lines[1].rstrip('\n').split(',')[1:]))
                    # Get index of freq for gridding
                    freq_chan = np.where(np.abs(freqs-freq).min()==np.abs(freqs-freq))[0]
                    # Filter freqs around freq_chan
                    freqs = freqs[freq_chan-width:freq_chan+width]
                for line in lines[2:]:
                    if line.startswith('#'):
                        continue
                    line = line.rstrip('\n').split(',')
                    if len(line) == (NFFT+1): # Make sure line has finished printing
                        spec_times.append(float(line[0]))
                        spec_raw.append(map(float,line[freq_chan-width+1:freq_chan+width+1]))
        spec_times = Time(spec_times,format='unix')
        spec_raw = np.array(spec_raw)
        freqs = np.array(freqs).squeeze()
        return spec_times,spec_raw,freqs,freq_chan

    elif filetype == 'echo':
        all_Data = []
        freqs = []
        echo_files = glob.glob(infile)
        for echo_file in echo_files:
            lines = open(echo_file).readlines()
            lat0,lon0 = map(float,lines[2].rstrip('\n').split(':')[1].strip(' ').split(','))
            freqs = map(float,lines[3].rstrip('\n').split(':')[1].strip(' ').split(','))
            #freqs = np.array(freqs)
            for line in lines:
                if line.startswith('#'):
                    continue
                line = line.rstrip('\n').split(',')
                if len(line) == (len(freqs)+4):
                    if not line[1] == '-1':
                        all_Data.append(map(float,line))
        all_Data = np.array(all_Data)
        spec_times,lats,lons,alts,spec_raw = (all_Data[:,0],all_Data[:,1],\
                                         all_Data[:,2],all_Data[:,3],\
                                         all_Data[:,4:])
        spec_times = Time(spec_times,format='gps')
        lats = lats.squeeze(); lats = np.insert(lats,0,lat0)
        lons = lons.squeeze(); lons = np.insert(lons,0,lon0)
        alts = alts.squeeze()
        freqs = np.array(freqs)
        return spec_times,spec_raw,freqs,lats,lons,alts#,lat0,lon0

    elif filetype == 'orbcomm':
        all_Data = []
        lines = open(infile).readlines()
        for line in lines: # Data begins on fifth line of accumulated file
            if line.startswith('#'):
                continue
            elif not line.split(',')[1] == '-1':
                all_Data.append(map(float,line.rstrip('\n').split(',')))
        all_Data = np.array(all_Data)
        spec_times,lats,lons,alts,yaws = (all_Data[:,1],all_Data[:,2],\
                                          all_Data[:,3],all_Data[:,4],\
                                          all_Data[:,5])
        spec_raw = all_Data[:,6:]
        spec_times = Time(spec_times,format='gps')
        lats = lats.squeeze()
        lons = lons.squeeze()
        alts = alts.squeeze()
        yaws = yaws.squeeze()
        return spec_times,spec_raw,lats,lons,alts,yaws

    elif filetype == 'start-stop':
        time_ranges = []
        lines = open(infile).readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            line = line.split(' ')
            if not len(line) == 0:
                time_ranges.append(map(float,line[0:2]))
        time_ranges = np.array(time_ranges)[0]
        return time_ranges

    elif filetype == 'waypts':
        waypt_times = []
        lines = open(infile).readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            waypt_times.append(float(line))
        waypt_times = np.array(waypt_times).squeeze()
        return waypt_times


    else:
        print '\nNo valid filetype found for %s' %infile
        print 'Exiting...\n\n'
        sys.exit()


def get_start_stop_times(infile):
    # infile can be filename or glob

    start_stop_times = []
    apm_files = glob.glob(infile)
    for apm_file in apm_files:
        lines=open(filename).readlines()
        weektimes = []
        for line in lines:
           if line.startswith('GPS'):
              weektimes.append(map(float,line.split(',')[3:5]))
        weektimes = np.array(weektimes)
        seconds = weektimes[:,1]*SEC_PER_WEEK + weektimes[:,0]/1000.
        times = Time(seconds, format='gps')
        mintime = times.gps.min()
        maxtime = times.gps.max()
        start_stop_times.append([mintime,maxtime])
    return start_stop_times


def get_way(infile):
    lines=open(infile).readlines()
    GPS_weektimes,GPS_arm,CMD_time,CMD_num =[],[],[],[]
    for line in lines[630:]:
        if line.startswith('GPS'):
            GPS_weektimes.append(map(float,line.split(',')[3:5]))
            GPS_arm.append(float(line.split(',')[1]))
        if line.startswith('CMD'):
            CMD_time.append(float(line.split(',')[1].strip()))
            CMD_num.append(int(line.split(',')[3].strip()))
    GPS_weektimes = np.array(GPS_weektimes)
    GPS_seconds = GPS_weektimes[:,1]*SEC_PER_WEEK + GPS_weektimes[:,0]/1000.
    GPS_arm= Time((np.array(GPS_arm[0])*APMLOG_SEC_PER_TICK), format = 'gps')
    GPS_time = Time(GPS_seconds, format='gps')
    CMD_time = (np.array(CMD_time).astype(float))*APMLOG_SEC_PER_TICK
    return GPS_time, GPS_arm, np.array(CMD_num), CMD_time


def get_filter_times(infile,first_waypt=3,waypts=False):
    # infile can be filename or glob
    waypoint_times = []
    start_stop_times = []
    apm_files = glob.glob(infile)
    for apm_file in apm_files:
        GPS_times,GPS_arm,CMD_num,CMD_times = get_way(apm_file)
        CMD_times = Time((CMD_times+(GPS_times.gps[0]-GPS_arm.gps)),
                         format='gps')

        for k,CMD in enumerate(CMD_num):
            if CMD==first_waypt:
                start = int(np.round((CMD_times.gps[k]),0))
            if CMD==CMD_num.max():
                stop = int(np.ceil((CMD_times.gps[k])))
                duration =  int(np.ceil((CMD_times.gps[k]))) - start

        start_stop_times.append([start,stop,duration])
        if waypts:
            for i in range(1,CMD_times.shape[0]):
                waypoint_times.append(CMD_times[i].gps)
    print start_stop_times
    start_stop_times = np.array(start_stop_times)
    if waypts:
        waypoint_times = np.array(waypoint_times)
        return start_stop_times,waypoint_times
    else:
        return start_stop_times
