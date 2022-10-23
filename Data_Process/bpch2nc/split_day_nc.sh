#!/bin/bash
# Lixu, 20200825
# this script is used for split netcdf data into daily data, it can be called from run_split_day_nc PBS shell script

sub1="201807"
sub2="201808"
sub3="201809"

#chagne the file name my_bpch_data later by modified in nc_example.py

for entry in /glade/scratch/lixujin/nested_bpch/trac_avg.geosfp_025x03125_tropchem_na.*_diagns.nc
do
  #echo $entry
  if [[ "$entry" == *"$sub1"* ]]; then
  #there will be error info if code has space around '='
    #output='my_bpch_data_ok.201807'
    output='trac_avg.geosfp_025x03125_tropchem_na.201807'
  fi
  
  if [[ "$entry" == *"$sub2"* ]]; then
    #output='my_bpch_data_ok.201808'
    output='trac_avg.geosfp_025x03125_tropchem_na.201808'
  fi

  if [[ "$entry" == *"$sub3"* ]]; then
    #output='my_bpch_data_ok.201809'
    output='trac_avg.geosfp_025x03125_tropchem_na.201809'
  fi
  
  echo "$entry"
  echo "$output"
  cdo -splitday $entry $output
done


#try to change all file names into a standard format
##my_bpch_data_ok.2018 -> 2+4+4+2+4+4x1=6+6+8=20; 20+4=24
##file: trac_avg.geosfp_025x03125_tropchem_na.20180901.nc -> 4+3+6+3+5+8+2+6+1x7=44
##file: split_trac_avg.geosfp_025x03125_tropchem_na.201809 -> 5+4+3+6+3+5+8+2+6+1x8=
#surffix=".nc"
#for file in split*.nc
#do
##  #20 may be modified later
##  mv "$file" "trac_avg.geosfp_025x03125_tropchem_na.${file:36:44}0000.nc"
#  
#  echo ${file:44:45}
#
#done
