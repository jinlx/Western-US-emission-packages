# prerequisite: load cdo package, it could be installed by spack install or other methods. 
# This is used to make monthly average data to speed up calculation
cdo monmean GFAS_201901.nc GFAS_201901_monthly.nc
cdo monmean GFAS_201902.nc GFAS_201902_monthly.nc
cdo monmean GFAS_201903.nc GFAS_201903_monthly.nc
cdo monmean GFAS_201904.nc GFAS_201904_monthly.nc
cdo monmean GFAS_201905.nc GFAS_201905_monthly.nc
cdo monmean GFAS_201906.nc GFAS_201906_monthly.nc
cdo monmean GFAS_201907.nc GFAS_201907_monthly.nc
cdo monmean GFAS_201908.nc GFAS_201908_monthly.nc
cdo monmean GFAS_201909.nc GFAS_201909_monthly.nc
cdo monmean GFAS_201910.nc GFAS_201910_monthly.nc
cdo monmean GFAS_201911.nc GFAS_201911_monthly.nc
cdo monmean GFAS_201912.nc GFAS_201912_monthly.nc

# This is used to fix GFAS montly data lon issue
cdo sellonlatbox,-180,180,-90,90 GFAS_201901_monthly.nc GFAS_201901_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201902_monthly.nc GFAS_201902_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201903_monthly.nc GFAS_201903_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201904_monthly.nc GFAS_201904_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201905_monthly.nc GFAS_201905_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201906_monthly.nc GFAS_201906_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201907_monthly.nc GFAS_201907_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201908_monthly.nc GFAS_201908_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201909_monthly.nc GFAS_201909_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201910_monthly.nc GFAS_201910_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201911_monthly.nc GFAS_201911_monthly_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201912_monthly.nc GFAS_201912_monthly_grid.nc

# This used to fix raw GFAS daily data (optional)
# Only used it if we want to check specific days
cdo sellonlatbox,-180,180,-90,90 GFAS_201901.nc GFAS_201901_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201902.nc GFAS_201902_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201903.nc GFAS_201903_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201904.nc GFAS_201904_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201905.nc GFAS_201905_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201906.nc GFAS_201906_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201907.nc GFAS_201907_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201908.nc GFAS_201908_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201909.nc GFAS_201909_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201910.nc GFAS_201910_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201911.nc GFAS_201911_grid.nc
cdo sellonlatbox,-180,180,-90,90 GFAS_201912.nc GFAS_201912_grid.nc
