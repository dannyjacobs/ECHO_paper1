import numpy as np
import healpy as hp
import optparse,sys
import matplotlib.gridspec as gridspec

from matplotlib.pyplot import *
from scipy.interpolate import interp1d

from read_utils import *
from position_utils import *
from plot_utils import *
from time_utils import *

o = optparse.OptionParser()
o.set_description('Creates plots for ECHO paper, yo.')

o.add_option('--acc_file',type=str,
    help='File containing times, positions, and spectra')
o.add_option('--freq',type=float,
    help='Frequency of interest for plotting')
o.add_option('--gb_data',action='store_true',
    help='If passed, use Green Bank data format')
o.add_option('--gb_waypts',type=str,
    help='Waypoint file for Green Bank flight of S ant, EW dipole with EW trans')
o.add_option('--gb_flight_times',type=str,
    help='Flight start and stop times file')

opts,args = o.parse_args(sys.argv[1:])

if opts.gb_data:
    spec_times,spec_raw,lats,lons,alts,yaws = get_data(opts.acc_file,
                                                                                filetype='orbcomm')
    spec_times = np.expand_dims(spec_times.gps,axis=1)
    lats = np.expand_dims(lats,axis=1)
    lons = np.expand_dims(lons,axis=1)
    alts = np.expand_dims(alts,axis=1)
    yaws = np.expand_dims(yaws,axis=1)

    for i in range(spec_raw.shape[0]):
        for j in range(spec_raw.shape[1]):
            spec_raw[i,j] = 10*np.log10(spec_raw[i,j])
    all_Data = np.concatenate((spec_times,lats,lons,alts,yaws,spec_raw),
                                            axis=1)

    #
    N_lat0,N_lon0 = (38.4248532,-79.8503723)
    N_NS_spec = all_Data[:,11:16] # N antenna, NS dipole
    N_EW_spec = all_Data[:,23:28] # N antenna, EW dipole
    S_lat0,S_lon0 = (38.4239235,-79.8503418)
    S_NS_spec = all_Data[:,5:10] # S antenna, NS dipole
    S_EW_spec = all_Data[:,17:22] # S antenna, EW dipole

    # Read in waypoint times for S ant EW dipole EW trans flight
    waypt_times = get_data(opts.gb_waypts,filetype='waypts')
    print '\nRead in ' + str(waypt_times.shape[0]) + ' waypoints...'
    #print spec_times.shape,S_EW_spec.shape
    S_EW_freq_chan = np.argmax(S_EW_spec[0,:])
    S_EW_spec_interp = interp1d(spec_times.squeeze(),S_EW_spec[:,S_EW_freq_chan],kind='zero')
    S_EW_spec_interp = S_EW_spec_interp(waypt_times)
    #print S_EW_spec.shape,S_EW_spec_interp.shape

    # Zipper data together for filtering

    # Read in start and stop times for flights
    start_stop_times = get_data(opts.gb_flight_times,filetype='start-stop')
    print '\nRead in %s flights' %start_stop_times.shape[0]

    ''' Only reading in two flights out of six.  WTF DUDE??? '''

    start_stop_flag_inds = flight_time_filter(start_stop_times,spec_times)
    print '\nShape of S_EW_spec before flight time filtering: ' +str(S_EW_spec.shape)
    print S_EW_spec.shape
    S_EW_spec_filtered = S_EW_spec[start_stop_flag_inds,:]
    spec_times_filtered = spec_times[start_stop_flag_inds]
    print 'Shape of S_EW_spec after flight time filtering: ' +str(S_EW_spec_filtered.shape)

    # Get waypt filtered spectrum
    waypt_flag_inds = waypt_time_filter(waypt_times,spec_times_filtered)
    print '\nShape of S_EW_spec before waypoint filtering: ' +str(S_EW_spec_filtered.shape)
    S_EW_spec_filtered = S_EW_spec_filtered[waypt_flag_inds,:]
    spec_times_filtered = spec_times_filtered[waypt_flag_inds]
    print 'Shape of S_EW_spec after waypoint filtering: ' +str(S_EW_spec_filtered.shape)

    #freq_chan = 2 # Green bank spectra are five channels wide, with 137.5 MHz in the middle

    # Get beam of S ant, EW dipole for beam
    freq_chan = np.argmax(S_EW_spec_filtered[0,:])
    hpx_beam,hpx_counts,hpx_rms = make_beam(lats,lons,alts,S_EW_spec_filtered,
                                                                         lat0=S_lat0,
                                                                         lon0=S_lon0,
                                                                         freq_chan=freq_chan)
    # Setup figure 1 for time_series plots
    fig1 = figure(figsize=((16,9)),facecolor='w',edgecolor='w')
    gs1 = gridspec.GridSpec(2,1)

    # Make time series plot
    minspec = S_EW_spec.min()
    maxspec = S_EW_spec.max()
    time_series_plot = fig1.add_subplot(gs1[0])
    time_series_plot.plot(spec_times,S_EW_spec[:,freq_chan],
                                    c='k',
                                    marker='o',
                                    linestyle='None',
                                    label='Raw Spectrum')
    time_series_plot.plot(waypt_times,S_EW_spec_interp,
                                    color='r',
                                    marker='x',
                                    linestyle='None',
                                    label='Waypoints')
    time_series_plot.vlines(start_stop_times[:,0],-100,100,
                                       colors='b',
                                       lw=2,
                                       label='Flight Start')
    time_series_plot.vlines(start_stop_times[:,1],-100,100,
                                       colors='g',
                                       lw=2,
                                       label='Flight Stop')
    time_series_plot.legend(loc='best')
    time_series_plot.set_ylim([minspec,maxspec])

    # Make filtered spectrum plot
    filtered_spec_plot = fig1.add_subplot(gs1[1],
                                                            sharex=time_series_plot,
                                                            sharey=time_series_plot)
    filtered_spec_plot.plot(spec_times_filtered,S_EW_spec_filtered[:,freq_chan],
                                     color='k',
                                     marker='o',
                                     linestyle='None',
                                     label='Filtered Spectrum')
    filtered_spec_plot.legend(loc='best')
    filtered_spec_plot.set_ylim([minspec,maxspec])

    # Adjust subplots
    fig1.subplots_adjust(hspace=0)
    setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)

    show()
