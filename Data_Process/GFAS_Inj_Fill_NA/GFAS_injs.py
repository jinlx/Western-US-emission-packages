# Lixu, 09/10/2021
# Raw script is from Jeff Pierce 
# This scirpt is used to create emission-flux weighted injections at corresponding GEOS-Chem resolution.
# The emission-flux weighted injections will help us tackle the issue the pixels without emissions would reduce the heights, if they are just simply averaged into a coarse grid.
from pylab import *
from netCDF4 import Dataset
import datetime

# grid setting
#grid='4x5'
#grid='2x2.5'
#grid='0.5x0.667'
#grid='0.5x0.625_NA'
grid='0.25x0.3125_NA'



#### Processing data from GC-ready data
# Directory setting
direc_input='/home/lj152920/project/geos-chem/ExtData/HEMCO/GFAS/v2018-09/2021/'
direc_output='/home/lj152920/project/geos-chem/ExtData/HEMCO/GFAS/v2021-09/2021/'

# Date setting
#yrmo=['201807','201808','201809','201810','201811','201812']
#yrmo=['201907','201908','201909']
#yrmo=['201808']
#yrmo=['202001','202002','202003']
#yrmo=['202001','202002','202003','202004','202005','202006','202007','202008','202009','202010','202011','202012']
#yrmo=['201801','201802','201803','201804','201805','201806','201807','201808','201809','201810','201811','201812']

#yrmo=['201901','201902','201903','201904',
#      '201905','201906','201907','201908',
#      '201909','201910','201911','201912']


yrmo=['202101','202102','202103','202104','202105','202106','202107','202108','202109','202110','202111','202112']



if grid == '4x5':
   late=empty((47))
   lone=empty((73))
   late[0]=-90.
   late[-1]=90.
   late[1:-1]=arange(-88.,90.,4.)
   lone[0]=0.
   lone[-1]=360.
   lone[1:-1]=arange(5.,360.,5.)
   switch_lon = False
elif grid == '0.5x0.667': # maybe has alignment issue
   late=empty((361))
   lone=empty((541))
   late[0]=-90.
   late[-1]=90.
   late[1:-1]=arange(-89.5,90.,0.5)
   lone[0]=0.
   lone[-1]=360.
   lone[1:-1]=arange(2./3.,360.,2./3.)
   switch_lon = False
elif grid == '2x2.5':
   late=empty((92))
   lone=empty((145))
   late[0]=-90.
   late[-1]=90.
   late[1:-1]=arange(-89.,90.,2.)
   lone[0]=0.
   lone[-1]=360.
   lone[1:-1]=arange(2.5,360.,2.5)
   switch_lon = False
elif grid == '0.5x0.625_NA':
   late=empty((122))
   lone=empty((162))
   late[0]=9.75
   late[-1]=70.25
   late[1:-1]=arange(10.25,70.,0.5)
   lone[0]=219.6875
   lone[-1]=320.3125
   lone[1:-1]=arange(219.6875+0.625,320.,0.625)
   switch_lon = True
elif grid == '0.25x0.3125_NA':
   late=empty((203))
   lone=empty((226))
   late[0]=9.75-0.25/2.
   late[-1]=60.+0.25/2.
   late[1:-1]=arange(9.75+0.25/2,60.,0.25)
   lone[0]=360.-130.-0.3125/2.
   lone[-1]=360.-60.+0.3125/2.
   lone[1:-1]=arange(360.-130.+0.3125/2.,360.-60.,0.3125)
   switch_lon = True
latc=(late[0:-1]+late[1:])/2.
lonc=(lone[0:-1]+lone[1:])/2.



for yr in yrmo:
   fname=direc_input + 'GFAS_' + yr  + '.nc'
   fid=Dataset(fname,'r')
   lata=fid.variables['lat']
   lat=lata[:]
   lona=fid.variables['lon']
   lon=lona[:]
   timea=fid.variables['time']
   time=timea[:]
   nt = len(time)

   mamia=fid.variables['mami']
   mami=mamia[:]

   ocfirea=fid.variables['ocfire']
   ocfire=ocfirea[:]

   mami_avg=empty((nt,len(latc),len(lonc)))
   mami_avg_weight=empty((nt,len(latc),len(lonc)))
   ocfire_avg=empty((nt,len(latc),len(lonc)))

   for y in range(0,len(latc)):
     print(y)
     for x in range(0,len(lonc)):
       islat=where(lat[:]>late[y])[0][0]
       inlat=where(lat[:]<late[y+1])[0][-1]+1
       iwlon=where(lon[:]>lone[x])[0][0]
       ielon=where(lon[:]<lone[x+1])[0][-1]+1
       #for t in range(0,nt):
       mamibox=mami[:,islat:inlat,iwlon:ielon]
       ocfirebox=ocfire[:,islat:inlat,iwlon:ielon]
    
       mami_avg_weight[:,y,x] = (mamibox*ocfirebox).sum(2).sum(1)/ocfirebox.sum(2).sum(1)
       mami_avg[:,y,x] = mamibox.mean(2).mean(1)
       ocfire_avg[:,y,x] = ocfirebox.mean(2).mean(1)
        

   mami_avg_weight=ma.masked_where(ocfire_avg == 0., mami_avg_weight)
   mami_avg=ma.masked_where(ocfire_avg == 0., mami_avg)
 

   #----------------------------------------------------------------
   # Save data to files
   #----------------------------------------------------------------
   # Make netCDF file
   newfile = direc_output+'GFAS_'+yr+'_hgts_'+grid+'.nc'
   nc_w_fid = Dataset(newfile, 'w', clobber=True,  format='NETCDF4')
   nc_w_fid.description = 'Emission fluxes weighted injection heights'
   nc_w_fid.history = 'Created ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
   # define file dimensions
   nc_w_fid.createDimension('time', None) # unlimited dimension
   nc_w_fid.createDimension('lon', len(lonc))
   nc_w_fid.createDimension('lat', len(latc))
   
   # create identity variables
   time_w = nc_w_fid.createVariable('time', np.float32, ('time',))
   lat_w = nc_w_fid.createVariable('lat', np.float32, ('lat',))
   lon_w = nc_w_fid.createVariable('lon', np.float32, ('lon',))

   if switch_lon:
      lon_w[:] = lonc-360.
   else:
      lon_w[:] = lonc
   lat_w[:] = latc
   time_w[:] = time
   
   # Create height variables
   mami_avg_weight_w = nc_w_fid.createVariable('mami_avg_weight', 
                                               np.float32,
                                               ('time', 'lat','lon',),
                                               fill_value=-1.e-31) 

   mami_avg_weight_w[:,:,:] = mami_avg_weight

   
   # Set attributes
   for a in lata.ncattrs():
      lat_w.setncattr(a,lata.getncattr(a))
   for a in lona.ncattrs():
      lon_w.setncattr(a,lona.getncattr(a))
   for a in timea.ncattrs():
      time_w.setncattr(a,timea.getncattr(a))
   for a in mamia.ncattrs()[1:]:
      mami_avg_weight_w.setncattr(a,mamia.getncattr(a))
   
   nc_w_fid.close()
