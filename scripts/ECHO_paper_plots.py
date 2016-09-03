import numpy as np
import healpy as hp
import optparse,sys,glob,warnings
import matplotlib.gridspec as gridspec

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d

from ECHO.read_utils import *
from ECHO.position_utils import *
from ECHO.plot_utils import *
from ECHO.time_utils import *

o = optparse.OptionParser()
o.set_description('Creates plots for ECHO paper, yo.')

o.add_option('--acc_file',type=str,
    help='File containing times, positions, and spectra')
o.add_option('--freq',type=float,
    help='Frequency of interest for plotting')
o.add_option('--apm_file',action='append',
    help='Filename or glob for APM data')
o.add_option('--spec_file',type=str,
    help='Filename or glob for Signal Hound data')
o.add_option('--gb_data',action='store_true',
    help='If passed, use Green Bank data format')
o.add_option('--gb_waypts',type=str,
    help='Waypoint file for Green Bank flight of S ant, EW dipole with EW trans')
o.add_option('--gb_flight_times',type=str,
    help='Flight start and stop times file')
o.add_option('--fits',type=str,
    help='Filename or glob for FITS beams plotting')

opts,args = o.parse_args(sys.argv[1:])

FONT_SIZE = 18

if opts.fits:
    def init_plot(fig,gs,coll,unit,cmap,fontsize):
        ax = fig.add_subplot(gs,aspect='equal')
        coll.set_cmap(cmap)
        ax.add_collection(coll)
        ax.autoscale_view()
        # Position colorbar next to plot with same height as plot
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        if unit == 'Counts':
            max = coll.get_array().max()
            ticks = np.linspace(1,max,max)
            if 15 < max <= 30:
                ticks = ticks[::2]
            elif 30 < max:
                ticks = ticks[::4]
            cbar = fig.colorbar(coll,
                                        cax=cax,
                                        use_gridspec=True,
                                        ticks=ticks)
        else:
            cbar = fig.colorbar(coll,
                                        cax=cax,
                                        use_gridspec=True)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.set_label(label=unit,size=fontsize)
        for radius_deg in [20,40,60,80]:
            r = np.sin(radius_deg*np.pi/180.)
            x = np.linspace(-r,r,100)
            ax.plot(x,np.sqrt(r**2-x**2),'w-',linewidth=3)
            ax.plot(x,-np.sqrt(r**2-x**2),'w-',linewidth=3)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        return ax

    infiles = glob.glob(opts.fits)
    for i in range(len(infiles)):
        infile = infiles[i]
        print '\nReading in %s...' %infile

        info = infile.split('_')[5:]
        info[-1] = info[-1].rstrip('.fits')

        beam = hp.read_map(infile)
        beam_coll = make_polycoll(beam)
        counts = hp.read_map(infile.replace('power','counts'))
        counts_coll = make_polycoll(counts)
        rms = hp.read_map(infile.replace('power','rms'))
        rms_coll = make_polycoll(rms)

        gs = gridspec.GridSpec(1,6)
        fig = figure(figsize=(16,5))
        title = ' '.join(info)
        fig.suptitle(title,size=FONT_SIZE)
        beam_plot = init_plot(fig,
                                        gs[:2],
                                        beam_coll,
                                        'dB',
                                        cm.gnuplot,
                                        FONT_SIZE)
        rms_plot = init_plot(fig,
                                      gs[2:4],
                                      rms_coll,
                                      'dB',
                                      cm.gnuplot,
                                      FONT_SIZE)
        counts_plot = init_plot(fig,
                                         gs[4:],
                                         counts_coll,
                                         'Counts',
                                         cm.bone,
                                         FONT_SIZE)

        for ax in fig.axes:
            ax.tick_params(axis='both', labelsize=FONT_SIZE)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            gs.tight_layout(fig, rect=[None, None, None, 0.97])

        outfile = '_'.join(info)+'.png'
        print 'Saving %s...' %outfile
        savefig(outfile)

        print '\n'

    show()







elif opts.gb_data:
    spec_times,spec_raw,lats,lons,alts,yaws = get_data(opts.acc_file,
                                                                                filetype='orbcomm')
    '''spec_times = spec_times.gps
    spec_times -= spec_times.min()
    print spec_times.min(),spec_times.max()
    sys.exit()


    print spec_times.shape
    print spec_times.min(),spec_times.max()
    print np.diff(spec_times).min(),np.diff(spec_times).max()
    print np.argmax(np.diff(spec_times) == np.diff(spec_times).max())
    print spec_times[14718:14721]
    sys.exit()'''

    spec_times = np.expand_dims(spec_times.gps,axis=1)
    lats = np.expand_dims(lats,axis=1)
    lons = np.expand_dims(lons,axis=1)
    alts = np.expand_dims(alts,axis=1)
    yaws = np.expand_dims(yaws,axis=1)

    for i in range(spec_raw.shape[0]):
        for j in range(spec_raw.shape[1]):
            spec_raw[i,j] = 10*np.log10(spec_raw[i,j])
    spec_raw -= spec_raw.max()

    all_Data = np.concatenate((spec_times,lats,lons,alts,yaws,spec_raw),
                                            axis=1)

    #
    N_lat0,N_lon0 = (38.4248532,-79.8503723)
    N_NS_spec = spec_raw[:,6:12] # N antenna, NS dipole
    N_EW_spec = spec_raw[:,18:] # N antenna, EW dipole
    S_lat0,S_lon0 = (38.4239235,-79.8503418)
    S_NS_spec = spec_raw[:,0:6] # S antenna, NS dipole
    S_EW_spec = spec_raw[:,12:18] # S antenna, EW dipole
    #N_NS_spec=N_NS_spec[:,3]
    #N_EW_spec=N_EW_spec[:,3]
    #S_NS_spec=S_NS_spec[:,3]
    #S_EW_spec=S_EW_spec[:,3]

    # Read in waypoint times for S ant EW dipole EW trans flight
    waypt_times = np.loadtxt(opts.gb_waypts,skiprows=1)
    print '\nRead in ' + str(waypt_times.shape[0]) + ' waypoints...'

    start_stop_times = np.loadtxt(opts.gb_flight_times,skiprows=1,delimiter=' ')
    print '\nRead in %s flights' %start_stop_times.shape[0]

    #print spec_times.shape,S_EW_spec.shape
    freq_chan = 3
    S_EW_spec_interp = interp1d(spec_times.squeeze(),S_EW_spec[:,freq_chan],kind='zero')
    S_EW_spec_interp = S_EW_spec_interp(waypt_times)
    #print S_EW_spec.shape,S_EW_spec_interp.shape

    start_stop_flag_inds = flight_time_filter(start_stop_times,spec_times).squeeze()
    print '\nShape of S_EW_spec before flight time filtering: ' +str(S_EW_spec.shape)
    S_EW_spec_filtered = S_EW_spec[start_stop_flag_inds]
    spec_times_filtered = spec_times[start_stop_flag_inds]
    print 'Shape of S_EW_spec after flight time filtering: ' +str(S_EW_spec_filtered.shape)


    # Get waypt filtered spectrum
    waypt_flag_inds = waypt_time_filter(waypt_times,spec_times_filtered)
    print waypt_flag_inds.shape
    print '\nShape of S_EW_spec before waypoint filtering: ' +str(S_EW_spec_filtered.shape)
    S_EW_spec_filtered = S_EW_spec_filtered[waypt_flag_inds]
    spec_times_filtered = spec_times_filtered[waypt_flag_inds]
    print 'Shape of S_EW_spec after waypoint filtering: ' +str(S_EW_spec_filtered.shape)

    #freq_chan = 2 # Green bank spectra are five channels wide, with 137.5 MHz in the middle

    # Setup figure 1 for time_series plots
    fig1 = figure(figsize=((16,9)),facecolor='w',edgecolor='w')
    gs1 = gridspec.GridSpec(3,1)

    # Make time series plot
    initmin = spec_times.min()
    waypt_times -= initmin
    spec_times_filtered -= initmin
    spec_times -= initmin
    start_stop_times -= initmin
    minspec = S_EW_spec.min()
    maxspec = S_EW_spec.max()
    ylim = [minspec,maxspec]
    mintime = spec_times.min()
    maxtime = spec_times.max()
    xlim = [mintime,maxtime]

    #zoom_cube_x = [112353200,112353225]
    zoom_cube_x = [1123531980,1123532060]
    zoom_cube_x -= initmin
    zoom_cube_y = [-10.5,-8]
    len_x = zoom_cube_x[1]-zoom_cube_x[0]
    len_y = np.abs(zoom_cube_y).max()-np.abs(zoom_cube_y).min()

    time_series_plot = fig1.add_subplot(gs1[0])
    time_series_plot.plot(spec_times,S_EW_spec[:,freq_chan],
                                    c='k',
                                    marker='o',
                                    linestyle='None',
                                    label='Raw Spectrum')
    '''time_series_plot.vlines(waypt_times,-100,100,
                                    color='0.5',
                                    lw=2,
                                    label='Waypoints')
    time_series_plot.vlines(waypt_times-3,-100,100,
                                    color='c',
                                    lw=2,
                                    label='Waypoints - 3')
    time_series_plot.vlines(waypt_times+3,-100,100,
                                    color='m',
                                    lw=2,
                                    label='Waypoints + 3')'''
    time_series_plot.vlines(start_stop_times[:,0],minspec,maxspec,
                                       colors='b',
                                       lw=2,
                                       label='Flight Start')
    time_series_plot.vlines(start_stop_times[:,1],minspec,maxspec,
                                       colors='g',
                                       lw=2,
                                       label='Flight Stop')
    '''
    time_series_plot.vlines(zoom_cube_x,minspec,maxspec,colors='k')
    time_series_plot.add_patch(Rectangle((zoom_cube_x[0],
                                                            minspec),
                                                            len_x,
                                                            maxspec-minspec,
                                                            facecolor='0.75',
                                                            edgecolor='0.75'))
    '''
    time_series_plot.legend(loc='best')
    time_series_plot.set_xlim(xlim)
    time_series_plot.set_ylim(ylim)
    time_series_plot.set_ylabel('Power [dB]')
    time_series_plot.set_xlabel('Time [s]')

    # Make filtered spectrum plot
    filtered_spec_plot = fig1.add_subplot(gs1[1])#,
                                                            #sharex=time_series_plot,
                                                            #sharey=time_series_plot)
    filtered_spec_plot.plot(spec_times,S_EW_spec[:,freq_chan],
                                     color='0.75',
                                     marker='o',
                                     linestyle='None',
                                     label='Raw Spectrum')
    filtered_spec_plot.plot(spec_times_filtered,S_EW_spec_filtered[:,freq_chan],
                                     color='k',
                                     marker='o',
                                     linestyle='None',
                                     label='Filtered Spectrum')
    filtered_spec_plot.vlines(waypt_times,-100,100,
                                    color='r',
                                    lw=2,
                                    label='Waypoints')
    filtered_spec_plot.set_ylabel('Power [dB]')
    filtered_spec_plot.set_xlabel('Time [s]')
    filtered_spec_plot.legend(loc='best')
    filtered_spec_plot.set_xlim(zoom_cube_x)
    filtered_spec_plot.set_ylim(zoom_cube_y)

    yaw_plot = fig1.add_subplot(gs1[2])
    yaw_plot.plot(spec_times,yaws,
                        color='b',
                        lw=2)
    yaw_plot.set_xlabel('Time [s]')
    yaw_plot.set_ylabel('Yaw [deg]')
    yaw_plot.set_xlim(zoom_cube_x)

    # Adjust subplots
    tight_layout()
    fig1.subplots_adjust(hspace=0.5)
    #setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)

    show()

else:
    def get_yaws(infile):
        yaw_times,yaws,yaw_errs = [],[],[]
        weektimes = []
        print 'Reading in %s...' %infile
        isGPS = True
        lines = open(infile).readlines()
        if not len(lines) == 0:
            for i,line in enumerate(lines):
                if line.startswith('GPS'):
                    weektimes.append(map(float,line.split(',')[3:5])) #ms and week number
                    if isGPS:
                        start_time = float(line.split(',')[1])/1000.
                        isGPS = False
                        print start_time
                if line.startswith('ATT') and not isGPS:
                   yaw_times.append(map(float,[line.split(',')[1].strip(' ')]))
                   yaws.append(map(float,[line.split(',')[7].strip(' ')]))
                   yaw_errs.append(float(line.split(',')[9]))
        weektimes = np.array(weektimes)
        apm_times = weektimes[:,1]*SEC_PER_WEEK+weektimes[:,0]/1000.
        apm_times = Time(apm_times,format='gps')
        yaw_times = np.array(yaw_times)*1e-6
        yaw_times = yaw_times+apm_times.gps[0]-start_time*1.0e-6
        yaw_times = Time(yaw_times, format = 'gps')
        yaws = np.array(yaws).squeeze()
        yaw_errs = np.array(yaw_errs).squeeze()
        return yaw_times,yaws,yaw_errs


    def get_way(infile):
        lines=open(infile).readlines()
        GPS_weektimes,GPS_arm,CMD_time,CMD_num =[],[],[],[]
        for line in lines:
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

        start_stop_times = np.array(start_stop_times)
        if waypts:
            waypoint_times = np.array(waypoint_times)
            return start_stop_times,waypoint_times
        else:
            return start_stop_times



    spec_times,spec_raw,freqs,freq_chan = get_data(opts.spec_file,filetype='sh')
    yt_1,y_1,ye_1 = get_yaws(opts.apm_file[0])
    yt_2,y_2,ye_2 = get_yaws(opts.apm_file[1])
    yt_3,y_3,ye_3 = get_yaws(opts.apm_file[2])
    '''
    print '\n\n'
    print yt_1.gps.min(),yt_2.gps.min(),yt_3.gps.min()
    print yt_1.gps.max(),yt_2.gps.max(),yt_3.gps.max()
    print '\n\n'
    '''
    s_1,w_1 = get_filter_times(opts.apm_file[0],waypts=True)
    s_2,w_2 = get_filter_times(opts.apm_file[1],waypts=True)
    s_3,w_3 = get_filter_times(opts.apm_file[2],waypts=True)
    yaw_times = np.concatenate([yt_1.gps,yt_2.gps,yt_3.gps])
    yaws = np.concatenate([y_1,y_2,y_3])
    yaw_errs = np.concatenate([ye_1,ye_2,ye_3])
    start_stop_times = np.concatenate([s_1,s_2,s_3])
    print ' '
    print start_stop_times.shape
    print ' '
    waypoint_times = np.concatenate([w_1,w_2,w_3])
    print ' '
    print waypoint_times.shape
    print waypoint_times.min(),waypoint_times.max()
    print ' '

    freq = 137.500
    freq_chan = np.where(np.abs(freqs-freq).min()==np.abs(freqs-freq))[0]
    #yaws = interp1d(yaw_times,yaws,kind='0')
    #yaws = yaws(spec_times)

    fig = figure(figsize=(16,9),facecolor='w',edgecolor='w')
    ax1 = fig.add_subplot(211)
    #ax1.plot(yt_1.gps-yt_1.gps.min(),y_1)
    #ax1.errorbar(yaw_times,yaws,yerr=yaw_errs,color='k',linestyle='-')
    ax1.plot(yaw_times,yaws,'k.-')
    ax1.vlines(waypoint_times,0,360,color='r')
    ax1.set_ylabel('yaw [deg]')


    print spec_times.gps.min(),spec_times.gps.max()
    print yaw_times.min(),yaw_times.max()

    print "\nFreq: %s\n" %freqs[freq_chan]

    ax2 = fig.add_subplot(212,sharex=ax1)
    minspec,maxspec = spec_raw[:,freq_chan].min(),spec_raw[:,freq_chan].max()
    ax2.plot(spec_times.gps,spec_raw[:,freq_chan],'b.-')
    ax2.vlines(waypoint_times,minspec,maxspec,color='r')
    ax2.set_ylabel('power [dB]')

    show()
