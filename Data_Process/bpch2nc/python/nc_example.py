# This python script is used to change bpch files into netcdf and extract variables
import xbpch
from glob import glob

#read a single file
#fns = "trac_avg.geosfp_025x03125_tropchem_na.201807190000"
#fns = 'ND49_20060101_ref_e2006_m2010.bpch'
#ds = xbpch.open_bpchdataset(fns)

#read multiple files
# List all the bpch files in the current directory
#fns = "trac_avg.geosfp_025x03125_tropchem_na.201809130000" #['ND49_1.bpch', 'ND49_2.bpch',''''''] 

fns = glob("trac_avg.geosfp_025x03125_tropchem_na.20180*0000") #['ND49_1.bpch', 'ND49_2.bpch',''''''] 

#data = xbpch.open_bpchdataset(fns)
#data = xbpch.common.fix_attr_encoding(data)
#data.to_netcdf('test.nc')

for index in range(len(fns)):
###for test
    time = fns[index].replace('trac_avg.geosfp_025x03125_tropchem_na', '')
#    #time_ok = time.replace('ND49_','')
    time_ok = time
    print('process '+time_ok)
#    print(fns[index])
    data = xbpch.open_bpchdataset(fns[index],categories=['IJ-AVG-$','BIOBSRCE','ANTHSRCE','BIOGSRCE'])
    data = xbpch.common.fix_attr_encoding(data)
    data.to_netcdf("trac_avg.geosfp_025x03125_tropchem_na" + time_ok+"_diagns"+".nc")

    
#change the input files, time, diagn, tracer information later!!

#time = ['201807180000', '201807190000', '201807200000', '201807210000', '201807220000', '201807230000', '201807240000',$
#            '201807250000', '201807260000', '201807270000', '201807280000', '201807290000', '201807300000', '201807310000',$
#            '201808010000', '201808020000', '201808030000', '201808040000', '201808050000', '201808060000', '201808070000',$
#            '201808080000', '201808090000', '201808100000', '201808110000', '201808120000', '201808130000', '201808140000',$
#            '201808150000', '201808160000', '201808170000', '201808180000', '201808190000', '201808200000', '201808210000',$
#            '201808220000', '201808230000', '201808240000', '201808250000', '201808260000', '201808270000', '201808280000',$
#            '201808290000', '201808300000', '201809010000', '201809020000', '201809030000', '201809040000', '201809050000',$
#            '201809060000', '201809070000', '201809070000', '201809080000', '201809090000', '201809100000', '201809110000',$
#            '201809120000', '201809130000', '201809140000', '201809150000', '201809160000']
