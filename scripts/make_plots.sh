#! /bin/bash

error_alert()
{
  if [ $? -ne 0 ] ;
  then
  echo "------------------------------------"
  echo "        ${1}"
  echo "------------------------------------"
  return
  fi
}

#make the plots for the ECHO paper

#GB DATA
#an example position time series
#power,rms,counts for orbcomm,North,NS


#yaw and alt plots
python plot_yaw.py '/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/apm_Aug13_v*log' ../figures/Sant_NSdipole_NStrans_yawalt.png
python plot_yaw.py "/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/EW_transmitter_polarization/apm_Aug*log" ../figures/Sant_EWdipole_EWtrans_yawalt.png
python plot_yaw.py "/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/NS_transmitter_polarization/apm_Aug*log" ../figures/Nant_NSdipole_NStrans_yawalt.png
python plot_yaw.py "/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/EW_transmitter_polarization/apm_Aug*log" ../figures/Nant_EWdipole_EWtrans_yawalt.png

#make a plot showing what interpolating on the power signal looks like
# note, zoom in manually the trace and save to GB_power_interpolation_zoom.png

python plot_GB_pos_power_interp.py
error_alert "error plotting power interpolation"
#make a plot showing yaw and waypoint flagging
# note, zoom in manually the trace and save to GB_waypoint_flagging_zoom.png
#python plot_GB_data_cmds.py



echo "-----------------------------"
echo making the accumalation files
echo "-----------------------------"
#make the ECHO accumulation file
 python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/satpowerflight12.0*' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/NS_transmitter_polarization/apm_Aug13_v*log' --lat0=38.4239235 --lon0=-79.8503418 --freq=137.5 --pol=NS_NS --rxant=S --acc_file=../data/acc_GB_2015_Sant_NStx_NSrx.txt
error_alert "error zipping together files"

python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/EW_transmitter_polarization/satpowerflight*' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/South_dipole/EW_transmitter_polarization/apm_*log' --lat0=38.4239235 --lon0=-79.8503418 --freq=137.5 --pol=EW_EW --rxant=S --acc_file=../data/acc_GB_2015_Sant_EWtx_EWrx.txt
error_alert "error zipping together files"

python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/NS_transmitter_polarization/satpowerflight*' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/NS_transmitter_polarization/apm_*log' --lat0=38.4248532 --lon0=-79.8503723 --freq=137.5 --pol=NS_NS --rxant=N --acc_file=../data/acc_GB_2015_Nant_NStx_NSrx.txt
error_alert "error zipping together files"

python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/EW_transmitter_polarization/satpowerflight*' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/Green_bank_Aug_2015/North_dipole/EW_transmitter_polarization/apm_*log' --lat0=38.4248532 --lon0=-79.8503723 --freq=137.5 --pol=EW_EW --rxant=N --acc_file=../data/acc_GB_2015_Nant_EWtx_EWrx.txt
error_alert "error zipping together files"

#make the beams
echo ------------------------------------
echo making beams
echo ------------------------------------
python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_GB_2015_Sant_NStx_NSrx.txt
error_alert "error making GB beams"
python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_GB_2015_Sant_EWtx_EWrx.txt
error_alert "error making GB beams"
python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_GB_2015_Nant_NStx_NSrx.txt
error_alert "error making GB beams"
python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_GB_2015_Nant_EWtx_EWrx.txt
error_alert "error making GB beams"


python plot_ECHO_GB_power_rms_counts.py \
   --savefig=../figures/GB_power_rms_count.png
error_alert "error plotting GB power rms etc"


#4 panel plot showing orbcomm maps
python plot_ECHO_GB_maps.py --min=-40 --max=0 \
   --savefig=../figures/GB_4_power_patterns_2.png
error_alert "error plotting GB power patterns"


#E/H slices for 4 dipoles
python plot_GB_slices.py --savefig=../figures/GB_slices_quad.pdf
error_alert "error plotting GB slices"



#GB COMPARED TO MODEL
#orbcomm with tx model subtracted
python ECHO_sub_tx_beam.py ../data/acc_GB_2015_Nant_EWtx_EWrx_8_beam.fits \
 ../data/bicolog_legs_360.fits
error_alert "error subtracting TX beam"

python ECHO_sub_tx_beam.py ../data/acc_GB_2015_Sant_EWtx_EWrx_8_beam.fits \
 ../data/bicolog_legs_360.fits
error_alert "error subtracting TX beam"

#python ECHO_sub_tx_beam.py ../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits \
# ../data/bicolog_legs_360.fits --theta=-10
python ECHO_sub_tx_beam.py ../data/acc_GB_2015_Nant_NStx_NSrx_8_beam.fits \
../data/137MHz_tx_farfield_rot10deg_dist9cm.fits
#../data/137MHz_farfield_10deg.fits
#the North antenna (Ant A) gets fit for a transmitter offset against antenna B
#python ECHO_tx_fit.py
error_alert "error subtracting TX beam"

python ECHO_sub_tx_beam.py ../data/acc_GB_2015_Sant_NStx_NSrx_8_beam.fits \
 ../data/bicolog_legs_360.fits
error_alert "error subtracting TX beam"


# #orbcomm with tx model removed
# python plot_ECHO_GB_maps.py --min=-40 --max=0 --savefig=../figures/GB_4_power_patterns_txcal.png --mode=usecal
# error_alert "error plotting calibrated GB maps"
#
#
# #error maps
# python plot_ECHO_GB_maps.py --min=-40 --max=0 --savefig=../figures/GB_4_power_patterns_txcal.png --mode=err
# error_alert "error plotting error"

#GB COMPARED TO ORBCOMM DATA
#ratio maps and slices with comparison to orbcomm
# make 6 ratio maps
python plot_ECHO_GB_ratios.py
error_alert "error plottin GB ratios"


#BUILD a best model of the GB dipoles

#average the calibrated beams together
python plot_GB_avg_slices.py
error_alert "error plotting avg GB slices"

return


#MWA DATA
#mwa tile, rms and counts

#make the accumulation files and then fits maps

#the first day of mapping
python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/MWA_Tile_July2016/22_June_2016/spec*txt' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/MWA_Tile_July2016/22_June_2016/apm*log' --lat0=34.6205940 --lon0=-112.4479370 --freq=137.5 --pol=EW_EW --rxant=MWA --acc_file=../data/acc_MWAtile_2016_run1_EWtx_EWrx.txt

python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_MWAtile_2016_run1_EWtx_EWrx.txt


#the second day of mapping to get the top part
python ~/src/ECHO/scripts/ECHO_zipper.py --rx_files='/Users/djacobs/Google_Drive/ECHO/Experiments/MWA_Tile_July2016/21_July_2016/spec*txt' --apm_files='/Users/djacobs/Google_Drive/ECHO/Experiments/MWA_Tile_July2016/21_July_2016/apm*log' --lat0=34.6205940 --lon0=-112.4479370 --freq=137.5 --pol=EW_EW --rxant=MWA  --acc_file=../data/acc_MWAtile_2016_run2_EWtx_EWrx.txt

python ~/src/ECHO/scripts/ECHO_mk_beam.py ../data/acc_MWAtile_2016_run2_EWtx_EWrx.txt


#Add the two maps together with an offset for the different attenuators used on the different days
python combine_MWAtile_maps.py
error_alert "error combining MWA tile maps"

#subtract tx model
python ECHO_sub_tx_beam.py ../data/acc_MWAtile_2016_EWtx_EWrx_8_beam.fits \
 ../data/bicolog_legs_360.fits
error_alert "error subtracting TX beam"

#make a healpix map of the mwa tile beam
cd ../data/
source activate MWA
python ../scripts/MWATilemodel2hpm.py --nside=32 --freq=137.5 --outname=MWA_Tile_nside32
error_alert "error making MWA model"
cd -
source activate ECHO

#mwa tile slices
python plot_MWAtile_slices.py --savefig=../figures/MWATile_slice.png
error_alert "error plotting mwa tile slices"
#MWA COMPARED TO ORBCOMM DATA
#comparison to mwa model
