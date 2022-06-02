# this is a manucript to calculate the ice cover duration for the recreation activities based on their thresholds
# by Lei-Huang 2021-11-10
import numpy as np
import xarray as xr
from tqdm import tqdm
import glob
from dask.distributed import Client
import pandas as pd
################################################################################################################################
# define functions to read data
''' remove the variables that is nothing to do with lake data'''
def def_process_coords(exceptcv=[]):
    def process_coords(ds, except_coord_vars=exceptcv):
        coord_vars = []
        for v in np.array(ds.coords):
            if not v in except_coord_vars:
                coord_vars += [v]
        for v in np.array(ds.data_vars):
            if not v in except_coord_vars:
                coord_vars += [v]
        return ds.drop(coord_vars)
    return process_coords

'''The large ensemble data are stored in one member by one member, and in the subdirectory of each member, one netcdf file has daily data for ten years
    it is a little tricky to read in the large ensemble data, this function literally only works on our HPC, your have to change your path accordingly
    after read in data, the lake ice array is supposed to have a dimension of [ensemble,time,lat,lon]
'''
def read_in(var, exceptcv, domain='lnd/', freq='day_1/', stream='h6', chunks=dict(time=365), ens_s=-20, ens_e=-10):
    ens_dir = "large_ensemble_data_directory"
    histens_names = [member.split('archive/')[1][:-1]
                     for member in sorted(glob.glob(ens_dir + "b.e21.BHIST*LE2*[!old][!tmp]/"))]  # [10:]
    projens_names = [member.split('archive/')[1][:-1] for member in sorted(
        glob.glob(ens_dir + "b.e21.BSSP370*.f09_g17*[!old][!tmp]/"))]  # [10:]
    hist_ncfiles = []
    proj_ncfiles = []
    for i in np.arange(ens_s, ens_e):
        hist_fnames = sorted(glob.glob(
            ens_dir + histens_names[i] + "/" + domain + "proc/tseries/" + freq + "*" + stream + var + "*"))
        proj_fnames = sorted(glob.glob(
            ens_dir + projens_names[i] + "/" + domain + "proc/tseries/" + freq + "*" + stream + var + "*"))
        hist_ncfiles.append(hist_fnames)
        proj_ncfiles.append(proj_fnames)
    ens_numbers = [members.split('LE2-')[1]
                   for members in histens_names][ens_s:ens_e]
    hist_ds = xr.open_mfdataset(hist_ncfiles,
                                chunks=chunks,
                                preprocess=def_process_coords(exceptcv),
                                combine='nested',
                                concat_dim=[[*ens_numbers], 'time'],
                                parallel=True,
                                decode_cf=False,
                                decode_times=True)
    proj_ds = xr.open_mfdataset(proj_ncfiles,
                                chunks=chunks,
                                preprocess=def_process_coords(exceptcv),
                                combine='nested',
                                concat_dim=[[*ens_numbers], 'time'],
                                parallel=True,
                                decode_cf=False,
                                decode_times=True)
    if freq == 'day_1/':
        hist_ds = hist_ds.isel(time=np.arange(1, hist_ds.time.shape[0]))
        proj_ds = proj_ds.isel(time=np.arange(1, proj_ds.time.shape[0]))
    #    hist_ds['time'] = hist_ds.time.get_index('time').shift(-1, 'D')
    #    proj_ds['time'] = proj_ds.time.get_index('time').shift(-1, 'D')
    #if freq == 'month_1/':
    #    hist_ds['time'] = hist_ds.time.get_index('time').shift(-1, 'M')
    #    proj_ds['time'] = proj_ds.time.get_index('time').shift(-1, 'M')
    ens_ds = xr.concat((hist_ds, proj_ds), 'time')
    ens_ds = ens_ds.rename({'concat_dim': 'ensemble'})
    return ens_ds


################################################################################################################################
## build a mask array transfer grid array to landunit array, thus save memory
domain = 'lnd/'
freq = 'month_1/'
var = '.EFLX_LH_TOT.'

ens_dir = "large_ensemble_data_directory"
histens_names = [member[-37:]
                 for member in sorted(glob.glob(ens_dir + "b.e21.BHISTcmip6*LE2*"))]
projens_names = [member[-39:]
                 for member in sorted(glob.glob(ens_dir + "b.e21.BSSP370cmip6.f09_g17*"))]
hist_ncfiles = []
proj_ncfiles = []
for i in np.arange(10, 30):
    hist_fnames = sorted(glob.glob(
        ens_dir + histens_names[i] + "/" + domain + "proc/tseries/" + freq + "*" + 'h2' + var + "*"))
    proj_fnames = sorted(glob.glob(
        ens_dir + projens_names[i] + "/" + domain + "proc/tseries/" + freq + "*" + 'h2' + var + "*"))
    hist_ncfiles.append(hist_fnames)
    proj_ncfiles.append(proj_fnames)
ens_numbers = [members[-8:] for members in histens_names][10:30]

mask_ds = xr.open_mfdataset(hist_ncfiles[0][0],
                            combine='by_coords',
                            parallel=True)

ixy = mask_ds.land1d_ixy - 1
ixy['landunit'] = mask_ds.land1d_ityplunit
lake_ixy = ixy.sel(landunit = 5).astype(int).compute()
jxy = mask_ds.land1d_jxy - 1
jxy['landunit'] = mask_ds.land1d_ityplunit
lake_jxy = jxy.sel(landunit = 5).astype(int).compute()
## change lake_ixy lake_jxy coordinate to make it consistent with computation afterwards
lake_ixy['landunit'] = np.arange(lake_ixy.shape[0])
lake_jxy['landunit'] = np.arange(lake_jxy.shape[0])

landunit_idx = mask_ds.landunit
landunit_idx['landunit'] = mask_ds.land1d_ityplunit
lake_idx = landunit_idx.sel(landunit=5).values

## calculate the lake area and make mask to remove antarctic and greenland
pct_lake_ds = xr.open_dataset(
    '/share/CESM/cesm_input/lnd/clm2/surfdata_map/surfdata_0.9x1.25_78pfts_CMIP6_simyr1850_c170824.nc')
pct_lake = pct_lake_ds.PCT_LAKE / 100
pct_lake = pct_lake.rename({'lsmlat': 'lat', 'lsmlon': 'lon'})
pct_lake = pct_lake.where(pct_lake > 0)
pct_lake['lat'] = mask_ds.lat
pct_lake['lon'] = mask_ds.lon

landfrac = pct_lake_ds.LANDFRAC_PFT
landfrac = landfrac.rename({'lsmlat': 'lat', 'lsmlon': 'lon'})
landfrac['lat'] = mask_ds.lat
landfrac['lon'] = mask_ds.lon

cell_area = pct_lake_ds.AREA
cell_area = cell_area.rename({'lsmlat': 'lat', 'lsmlon': 'lon'})
cell_area['lat'] = mask_ds.lat
cell_area['lon'] = mask_ds.lon

## setup mask to remove ocean and ice sheet area
glacier_region = pct_lake_ds.GLACIER_REGION
glacier_region = glacier_region.rename(dict(lsmlat='lat', lsmlon='lon'))
glacier_region['lat'] = mask_ds.lat
glacier_region['lon'] = mask_ds.lon
landfrac = landfrac.where(glacier_region == 0, 0)
pct_lake = pct_lake.where(landfrac > 0)
lake_area = cell_area * landfrac * pct_lake

lake_mask = (pct_lake_ds.PCT_LAKE > 0)
lake_mask = lake_mask.rename({'lsmlat': 'lat', 'lsmlon': 'lon'})
lake_mask['lat'] = mask_ds.lat
lake_mask['lon'] = mask_ds.lon
## remove lakes in antarctic and greenland
lake_mask = lake_mask.where(landfrac > 0, False)

lake_area_lunit = lake_area[lake_jxy, lake_ixy]

lake_mask_lunit = lake_mask[lake_jxy, lake_ixy]

lake_ixy = lake_ixy.sel(landunit=lake_mask_lunit)
lake_jxy = lake_jxy.sel(landunit=lake_mask_lunit)
#lake_idx = lake_idx[lake_mask_lunit]
#lake_area_lunit = lake_area_lunit.sel(landunit=lake_mask_lunit)

lake_area_lunit = lake_area[lake_jxy, lake_ixy]
lake_mask_lunit = lake_mask[lake_jxy, lake_ixy]
################################################################################################################################
# calculation starts
if __name__ == "__main__":
    client = Client(
        scheduler_file='your_scheduler_directory/scheduler.json')
    variables = ['LAKEICETHICK']
    exceptcv = ['time', 'lat', 'lon', *variables]
    ice_day_ds = read_in(var='LAKEICETHICK',
                         chunks={'time': 365},
                         exceptcv=exceptcv,
                         domain='lnd/',
                         freq='day_1/',
                         stream='*',
                         ens_s=10,
                         ens_e=100)
    time_index = xr.cftime_range(
        start='1850', end='2101', freq='1D', calendar='noleap')[:-1]
    ice_day_ds = ice_day_ds.assign_coords(time=time_index)
    # northern Hemisphere
    daily_ice = ice_day_ds.LAKEICETHICK
    # instead of calculating the phenology data in a complete canlendar year, here by shifting the time index, we are going to calculate the ice phenology for a complete winter season instead of a complete calendar year
    daily_ice['time'] = daily_ice.time.get_index('time').shift(-232, "D")
    daily_ice = daily_ice.sel(time=slice('1850-01-01', '2099-12-31'))
    daily_ice_lunit = daily_ice[:, :, lake_jxy, lake_ixy]
    years = np.unique(daily_ice_lunit.time.dt.year)
    dayofyear = np.unique(daily_ice_lunit.time.dt.dayofyear)
    ind = pd.MultiIndex.from_product(
        (years, dayofyear), names=('year', 'dayofyear'))
    daily_ice_res = daily_ice_lunit.assign_coords(time=ind).unstack('time')
    # southern Hemisphere
    lake_ixy_SH = lake_ixy[lake_ixy.lat < 0]
    lake_jxy_SH = lake_jxy[lake_jxy.lat < 0]
    daily_ice_SH = ice_day_ds.LAKEICETHICK.sel(
        time=slice('1850-01-01', '2099-12-31'), lat=slice(-90, 0))
    daily_ice_SH_lunit = daily_ice_SH[:,:, lake_jxy_SH, lake_ixy_SH]
    daily_ice_SH_res = daily_ice_SH_lunit.assign_coords(time=ind).unstack('time')
    '''
    #==========================================================================================
    # calculate 2 inch threshold
    dayofyear_mask_lunit_2inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 0.0508).transpose('ensemble', ...)
    iceon_mod_lunit_2inch = xr.DataArray(np.nan,
                                   attrs=daily_ice.attrs,
                                   name='LAKEICETHICK',
                                   dims=('ensemble', 'landunit', 'year'),
                                   coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_2inch = xr.DataArray(np.nan,
                                    attrs=daily_ice.attrs,
                                    name='LAKEICETHICK',
                                    dims=('ensemble', 'landunit', 'year'),
                                    coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_2inch = xr.DataArray(np.nan,
                                         attrs=daily_ice.attrs,
                                         name='LAKEICETHICK',
                                         dims=('ensemble', 'landunit', 'year'),
                                         coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_2inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_2inch[i *
                                                                          10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_2inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_2inch[i *
                                                                    10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_2inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_2inch[i *
                                                                 10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 2 inch threshold for southern hemisphere
    dayofyear_mask_lunit_2inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 0.0508).transpose('ensemble', ...)
    iceon_mod_lunit_2inch_SH = xr.DataArray(np.nan,
                                   attrs=daily_ice.attrs,
                                   name='LAKEICETHICK',
                                   dims=('ensemble', 'landunit', 'year'),
                                   coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_2inch_SH = xr.DataArray(np.nan,
                                    attrs=daily_ice.attrs,
                                    name='LAKEICETHICK',
                                    dims=('ensemble', 'landunit', 'year'),
        coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_2inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_2inch_SH[i *
                                                                                10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_2inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_2inch_SH[i *
                                                                                 10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_2inch[:, :644, :] = iceon_mod_lunit_2inch_SH.values
    iceoff_mod_lunit_2inch[:, :644, :] = iceoff_mod_lunit_2inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_2inch.transpose(..., 'landunit').rename('iceduration_2inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_2inch.nc')
    iceon_mod_lunit_2inch.transpose(..., 'landunit').rename('iceon_2inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_2inch.nc')
    iceoff_mod_lunit_2inch.transpose(..., 'landunit').rename('iceoff_2inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_2inch.nc')
    #==========================================================================================
    '''
    # calculate 4 inch threshold
    dayofyear_mask_lunit_4inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 0.1016).transpose('ensemble', ...) # unit of lakeicethick is meter
    iceon_mod_lunit_4inch = xr.DataArray(np.nan,
                                         attrs=daily_ice.attrs,
                                         name='LAKEICETHICK',
                                         dims=('ensemble', 'landunit', 'year'),
                                         coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_4inch = xr.DataArray(np.nan,
                                          attrs=daily_ice.attrs,
                                          name='LAKEICETHICK',
                                          dims=('ensemble', 'landunit', 'year'),
                                          coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_4inch = xr.DataArray(np.nan,
                                               attrs=daily_ice.attrs,
                                               name='LAKEICETHICK',
                                               dims=('ensemble',
                                                     'landunit', 'year'),
                                               coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    # since the data set is large, so we load the data into memory ten member by ten member
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_4inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_4inch[i *
                                                                                      10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_4inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_4inch[i *
                                                                                10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_4inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_4inch[i *
                                                                                 10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 4 inch threshold for southern hemisphere
    dayofyear_mask_lunit_4inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 0.1016).transpose('ensemble', ...)
    iceon_mod_lunit_4inch_SH = xr.DataArray(np.nan,
                                            attrs=daily_ice.attrs,
                                            name='LAKEICETHICK',
                                            dims=('ensemble',
                                                  'landunit', 'year'),
                                            coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_4inch_SH = xr.DataArray(np.nan,
                                             attrs=daily_ice.attrs,
                                             name='LAKEICETHICK',
                                             dims=('ensemble',
                                                   'landunit', 'year'),
                                             coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})

    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_4inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_4inch_SH[i *
                                                                                      10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_4inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_4inch_SH[i *
                                                                                       10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_4inch[:, :644, :] = iceon_mod_lunit_4inch_SH.values
    iceoff_mod_lunit_4inch[:, :644, :] = iceoff_mod_lunit_4inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_4inch.transpose(..., 'landunit').rename('iceduration_4inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_4inch.nc')
    iceon_mod_lunit_4inch.transpose(..., 'landunit').rename('iceon_4inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_4inch.nc')
    iceoff_mod_lunit_4inch.transpose(..., 'landunit').rename('iceoff_4inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_4inch.nc')
    #==========================================================================================
    # calculate 5 inch threshold
    dayofyear_mask_lunit_5inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 0.127).transpose('ensemble', ...)
    iceon_mod_lunit_5inch = xr.DataArray(np.nan,
                                         attrs=daily_ice.attrs,
                                         name='LAKEICETHICK',
                                         dims=('ensemble', 'landunit', 'year'),
                                         coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_5inch = xr.DataArray(np.nan,
                                          attrs=daily_ice.attrs,
                                          name='LAKEICETHICK',
                                          dims=('ensemble', 'landunit', 'year'),
                                          coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_5inch = xr.DataArray(np.nan,
                                               attrs=daily_ice.attrs,
                                               name='LAKEICETHICK',
                                               dims=('ensemble',
                                                     'landunit', 'year'),
                                               coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_5inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_5inch[i *
                                                                                      10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_5inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_5inch[i *
                                                                                10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_5inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_5inch[i *
                                                                                 10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 5 inch threshold for southern hemisphere
    dayofyear_mask_lunit_5inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 0.127).transpose('ensemble', ...)
    iceon_mod_lunit_5inch_SH = xr.DataArray(np.nan,
                                            attrs=daily_ice.attrs,
                                            name='LAKEICETHICK',
                                            dims=('ensemble',
                                                  'landunit', 'year'),
                                            coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_5inch_SH = xr.DataArray(np.nan,
                                             attrs=daily_ice.attrs,
                                             name='LAKEICETHICK',
                                             dims=('ensemble',
                                                   'landunit', 'year'),
                                             coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})

    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_5inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_5inch_SH[i *
                                                                                      10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_5inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_5inch_SH[i *
                                                                                       10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_5inch[:, :644, :] = iceon_mod_lunit_5inch_SH.values
    iceoff_mod_lunit_5inch[:, :644, :] = iceoff_mod_lunit_5inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_5inch.transpose(..., 'landunit').rename('iceduration_5inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_5inch.nc')
    iceon_mod_lunit_5inch.transpose(..., 'landunit').rename('iceon_5inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_5inch.nc')
    iceoff_mod_lunit_5inch.transpose(..., 'landunit').rename('iceoff_5inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_5inch.nc')
    #==========================================================================================
    '''
    # calculate 8 inch threshold
    dayofyear_mask_lunit_8inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 0.2032).transpose('ensemble', ...)
    iceon_mod_lunit_8inch = xr.DataArray(np.nan,
                                         attrs=daily_ice.attrs,
                                         name='LAKEICETHICK',
                                         dims=('ensemble', 'landunit', 'year'),
                                         coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_8inch = xr.DataArray(np.nan,
                                          attrs=daily_ice.attrs,
                                          name='LAKEICETHICK',
                                          dims=('ensemble', 'landunit', 'year'),
                                          coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_8inch = xr.DataArray(np.nan,
                                               attrs=daily_ice.attrs,
                                               name='LAKEICETHICK',
                                               dims=('ensemble',
                                                     'landunit', 'year'),
                                               coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_8inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_8inch[i *
                                                                                      10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_8inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_8inch[i *
                                                                                10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_8inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_8inch[i *
                                                                                 10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 8 inch threshold for southern hemisphere
    dayofyear_mask_lunit_8inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 0.2032).transpose('ensemble', ...)
    iceon_mod_lunit_8inch_SH = xr.DataArray(np.nan,
                                            attrs=daily_ice.attrs,
                                            name='LAKEICETHICK',
                                            dims=('ensemble',
                                                  'landunit', 'year'),
                                            coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_8inch_SH = xr.DataArray(np.nan,
                                             attrs=daily_ice.attrs,
                                             name='LAKEICETHICK',
                                             dims=('ensemble',
                                                   'landunit', 'year'),
                                             coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})

    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_8inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_8inch_SH[i *
                                                                                      10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_8inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_8inch_SH[i *
                                                                                       10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_8inch[:, :644, :] = iceon_mod_lunit_8inch_SH.values
    iceoff_mod_lunit_8inch[:, :644, :] = iceoff_mod_lunit_8inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_8inch.transpose(..., 'landunit').rename('iceduration_8inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_8inch.nc')
    iceon_mod_lunit_8inch.transpose(..., 'landunit').rename('iceon_8inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_8inch.nc')
    iceoff_mod_lunit_8inch.transpose(..., 'landunit').rename('iceoff_8inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_8inch.nc')
    '''
    #==========================================================================================
    # calculate 12 inch threshold
    dayofyear_mask_lunit_12inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 0.3048).transpose('ensemble', ...)
    iceon_mod_lunit_12inch = xr.DataArray(np.nan,
                                          attrs=daily_ice.attrs,
                                          name='LAKEICETHICK',
                                          dims=('ensemble', 'landunit', 'year'),
                                          coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_12inch = xr.DataArray(np.nan,
                                           attrs=daily_ice.attrs,
                                           name='LAKEICETHICK',
                                           dims=('ensemble',
                                                 'landunit', 'year'),
                                           coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_12inch = xr.DataArray(np.nan,
                                                attrs=daily_ice.attrs,
                                                name='LAKEICETHICK',
                                                dims=('ensemble',
                                                      'landunit', 'year'),
                                                coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_12inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_12inch[i *
                                                                                        10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_12inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_12inch[i *
                                                                                  10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_12inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_12inch[i *
                                                                                   10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 12 inch threshold for southern hemisphere
    dayofyear_mask_lunit_12inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 0.3048).transpose('ensemble', ...)
    iceon_mod_lunit_12inch_SH = xr.DataArray(np.nan,
                                             attrs=daily_ice.attrs,
                                             name='LAKEICETHICK',
                                             dims=('ensemble',
                                                   'landunit', 'year'),
                                             coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_12inch_SH = xr.DataArray(np.nan,
                                              attrs=daily_ice.attrs,
                                              name='LAKEICETHICK',
                                              dims=('ensemble',
                                                    'landunit', 'year'),
                                              coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})

    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_12inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_12inch_SH[i *
                                                                                        10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_12inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_12inch_SH[i *
                                                                                         10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_12inch[:, :644, :] = iceon_mod_lunit_12inch_SH.values
    iceoff_mod_lunit_12inch[:, :644, :] = iceoff_mod_lunit_12inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_12inch.transpose(..., 'landunit').rename('iceduration_12inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_12inch.nc')
    iceon_mod_lunit_12inch.transpose(..., 'landunit').rename('iceon_12inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_12inch.nc')
    iceoff_mod_lunit_12inch.transpose(..., 'landunit').rename('iceoff_12inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_12inch.nc')
    #==========================================================================================
    # calculate 42 inch threshold
    dayofyear_mask_lunit_42inch = daily_ice_res.dayofyear.where(
        daily_ice_res >= 1.0668).transpose('ensemble', ...)
    iceon_mod_lunit_42inch = xr.DataArray(np.nan,
                                          attrs=daily_ice.attrs,
                                          name='LAKEICETHICK',
                                          dims=('ensemble', 'landunit', 'year'),
                                          coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_42inch = xr.DataArray(np.nan,
                                           attrs=daily_ice.attrs,
                                           name='LAKEICETHICK',
                                           dims=('ensemble',
                                                 'landunit', 'year'),
                                           coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceduration_mod_lunit_42inch = xr.DataArray(np.nan,
                                                attrs=daily_ice.attrs,
                                                name='LAKEICETHICK',
                                                dims=('ensemble',
                                                      'landunit', 'year'),
                                                coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_lunit.landunit, 'year': np.arange(1850, 2100)})
    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceduration_mod_lunit_42inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_42inch[i *
                                                                                        10:(i+1)*10, ...].count('dayofyear')
        iceon_mod_lunit_42inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_42inch[i *
                                                                                  10:(i+1)*10, ...].min('dayofyear') + 232
        iceoff_mod_lunit_42inch[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_42inch[i *
                                                                                   10:(i+1)*10, ...].max('dayofyear') + 232
    # calculate 42 inch threshold for southern hemisphere
    dayofyear_mask_lunit_42inch_SH = daily_ice_SH_res.dayofyear.where(
        daily_ice_SH_res >= 1.0668).transpose('ensemble', ...)
    iceon_mod_lunit_42inch_SH = xr.DataArray(np.nan,
                                             attrs=daily_ice.attrs,
                                             name='LAKEICETHICK',
                                             dims=('ensemble',
                                                   'landunit', 'year'),
                                             coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})
    iceoff_mod_lunit_42inch_SH = xr.DataArray(np.nan,
                                              attrs=daily_ice.attrs,
                                              name='LAKEICETHICK',
                                              dims=('ensemble',
                                                    'landunit', 'year'),
                                              coords={'ensemble': daily_ice.ensemble, 'landunit': daily_ice_SH_lunit.landunit, 'year': np.arange(1850, 2100)})

    for i in tqdm(np.arange(0, 9), desc='1st loop'):
        iceon_mod_lunit_42inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_42inch_SH[i *
                                                                                        10:(i+1)*10, ...].min('dayofyear')
        iceoff_mod_lunit_42inch_SH[i*10:(i+1)*10, :, :] = dayofyear_mask_lunit_42inch_SH[i *
                                                                                         10:(i+1)*10, ...].max('dayofyear')
    # combine southern Hemisphere and northern Hemisphere
    iceon_mod_lunit_42inch[:, :644, :] = iceon_mod_lunit_42inch_SH.values
    iceoff_mod_lunit_42inch[:, :644, :] = iceoff_mod_lunit_42inch_SH.values
    # write to netcdf
    iceduration_mod_lunit_42inch.transpose(..., 'landunit').rename('iceduration_42inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_42inch.nc')
    iceon_mod_lunit_42inch.transpose(..., 'landunit').rename('iceon_42inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_42inch.nc')
    iceoff_mod_lunit_42inch.transpose(..., 'landunit').rename('iceoff_42inch').to_netcdf(
        '/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_42inch.nc')

    '''
    # put all the inches together
    xr.concat((iceduration_mod_lunit_2inch,iceduration_mod_lunit_4inch,iceduration_mod_lunit_5inch,
                iceduration_mod_lunit_8inch, iceduration_mod_lunit_12inch),
              dim='thickness').transpose(..., 'landunit').rename('iceduration').to_netcdf('/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceduration_mod_lunit_allinch.nc')
    xr.concat((iceon_mod_lunit_2inch, iceon_mod_lunit_4inch, iceon_mod_lunit_5inch,
               iceon_mod_lunit_8inch, iceon_mod_lunit_12inch),
              dim='thickness').transpose(..., 'landunit').rename('iceon').to_netcdf('/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceon_mod_lunit_allinch.nc')
    xr.concat((iceoff_mod_lunit_2inch, iceoff_mod_lunit_4inch, iceoff_mod_lunit_5inch,
               iceoff_mod_lunit_8inch, iceoff_mod_lunit_12inch),
              dim='thickness').transpose(..., 'landunit').rename('iceoff').to_netcdf('/proj/lhuang/scripts/CESM2-LENS/LAKE_ICE_Recreation/iceoff_mod_lunit_allinch.nc')
    '''


    # transfer array from landunit coordinates to geographical grid, take 4inch as an example, so was the same for other thresholds
    iceduration_grid_4inch = np.empty((90,250,192,288))*np.nan
    for i in np.arange(iceduration_mod_lunit_4inch.landunit.size):
        iceduration_grid_4inch[:,:,lake_jxy[i],lake_ixy[i]] = iceduration_mod_lunit_4inch[:,:,i]
    iceduration_grid_4inch = xr.DataArray(
                                iceduration_grid_4inch,
                                dims = ('ensemble','year','lat','lon'),
                                coords = {
                                    'ensemble':iceduration_mod_lunit_4inch.ensemble,
                                    'year':iceduration_mod_lunit_4inch.year,
                                    'lat':mask_ds.lat,
                                    'lon':mask_ds.lon
                                },
                                attrs={
                                    'long_name':'ice duration for ice thickness at least 4inch',
                                    'units':'day'
                                },
                                name='iceduration'
                                )
    iceduration_grid_4inch_ensmean = iceduration_grid_4inch.mean('ensemble').to_netcdf('output_file_path')
print('End')
