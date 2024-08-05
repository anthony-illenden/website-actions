from datetime import datetime,timedelta
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from xarray import open_dataset
from xarray.backends import NetCDF4DataStore
import xarray as xr
from siphon.catalog import TDSCatalog
from netCDF4 import Dataset
import metpy
from metpy.plots import USCOUNTIES
import matplotlib.colors as mcolors
import metpy.calc as mpcalc

def wf_850_temp(temp):
    temp_value = temp.item() if hasattr(temp, 'item') else temp
    if temp_value <= -25:
        return 0.6
    elif -25 < temp_value <= -20:
        return 0.8
    elif -20 < temp_value <= -16:
        return 1.0
    elif -16 < temp_value <= -14:
        return 1.6
    elif -14 < temp_value <= -12:
        return 1.4
    elif -12 < temp_value <= -10:
        return 1.2
    elif -10 < temp_value <= -9:
        return 1.0
    elif -9 < temp_value <= -8:
        return 0.9
    elif -8 < temp_value <= -7:
        return 0.8
    elif -7 < temp_value <= -6:
        return 0.7
    elif -6 < temp_value <= -5:
        return 0.6
    elif -5 < temp_value <= 0:
        return 0.5
    else:
        return 0

def wf_850_700_rh(rh):
    rh_value = rh.item() if hasattr(rh, 'item') else rh
    if rh_value <= 10:
        return 0.2
    elif 10 < rh_value <= 20:
        return 0.4
    elif 20 < rh_value <= 30:
        return 0.6
    elif 30 < rh_value <= 40:
        return 0.8
    elif 40 < rh_value <= 50:
        return 1.0
    elif 50 < rh_value <= 60:
        return 1.2
    elif 60 < rh_value <= 70:
        return 1.4
    elif 70 < rh_value <= 80:
        return 1.6
    elif 80 < rh_value <= 90:
        return 1.8
    else:
        return 2.0

def wf_1000_850_ws(ws):
    ws_value = ws.item() if hasattr(ws, 'item') else ws
    if ws_value <= 3:
        return 0.8
    elif 3 < ws_value <= 9:
        return 0.9
    elif 9 < ws_value <= 12:
        return 1.0
    elif 12 < ws_value <= 18:
        return 1.1
    elif 18 < ws_value <= 21:
        return 1.2
    elif 21 < ws_value <= 24:
        return 1.3
    elif 24 < ws_value <= 27:
        return 1.4
    elif 27 < ws_value <= 30:
        return 1.5
    else:
        return 1.6
    
def find_time_dim(ds, var_name):
    possible_time_dims = ['time', 'time1', 'time2', 'time3']
    time_dim = None
    for dim in possible_time_dims:
        if dim in ds[var_name].dims:
            time_dim = dim
            break
    if time_dim is None:
        raise ValueError('Could not find the time dimension')
    return time_dim
    
def find_press_dim(ds, var_name):
    possible_iso_dims = ['isobaric', 'isobaric1', 'isobaric2', 'isobaric3']
    iso_dim = None
    for dim in possible_iso_dims:
        if dim in ds[var_name].dims:
            iso_dim = dim
            break
    if iso_dim is None:
        raise ValueError('Could not find the iso dimension')
    return iso_dim

tds_gfs = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/GFS/Global_0p25deg/latest.html')
gfs_ds = tds_gfs.datasets[0]
ds = xr.open_dataset(gfs_ds.access_urls['OPENDAP'])
ds = ds.metpy.parse_cf()
ds_latlon = ds.sel(lat=slice(50, 20), lon=slice(360-130, 360-60))
print(ds_latlon.dims)
iso_dim = find_press_dim(ds_latlon, 'Temperature_isobaric')
time_dim = find_time_dim(ds_latlon, 'Temperature_isobaric')  
ds_latlon = ds_latlon.isel({time_dim: 0})
temp = ds_latlon['Temperature_isobaric'].sel({iso_dim: 85000}).metpy.convert_units('degC')
rh = ds_latlon['Relative_humidity_isobaric']
u = ds_latlon['u-component_of_wind_isobaric']
v = ds_latlon['v-component_of_wind_isobaric']

rh_850_700 = (rh.sel({iso_dim: 85000}) + rh.sel({iso_dim: 70000})) / 2
u_850 = u.sel({iso_dim: 85000})
v_850 = v.sel({iso_dim: 85000})
u_1000 = u.sel({iso_dim: 100000})
v_1000 = v.sel({iso_dim: 100000})
wspd_850 = mpcalc.wind_speed(u_850, v_850)
wspd_1000 = mpcalc.wind_speed(u_1000, v_1000)
wspd = ((wspd_850 + wspd_1000) / 2) * 1.94384

count = 0
for i in range(0, 49):
    wf_850_temp_vec = np.vectorize(wf_850_temp)
    wf_850_700_rh_vec = np.vectorize(wf_850_700_rh)
    wf_1000_850_ws_vec = np.vectorize(wf_1000_850_ws)

    temp_wf = xr.apply_ufunc(wf_850_temp_vec, temp[:,:])
    rh_wf = xr.apply_ufunc(wf_850_700_rh_vec, rh_850_700[:,:])
    wspd_wf = xr.apply_ufunc(wf_1000_850_ws_vec, wspd[:,:])

    lsp = temp_wf * rh_wf * wspd_wf

    cmap = mcolors.ListedColormap(['none', 'yellow', 'red'])
    bounds = [0, 1, 2]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=(12, 9), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-91, -81, 40.5, 47.75])
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.STATES.with_scale('50m'))
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)
    ax.add_feature(cfeature.LAKES, zorder=1, color='white')
    ax.add_feature(cfeature.LAND, color='white')
    cf = plt.contourf(lsp['lon'], lsp['lat'], lsp, cmap=cmap, levels=np.arange(0,3.1, 1), extend='max')
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical',fraction=0.046, pad=0.04)
    cbar.set_ticks([0.5, 1.5, 2.5])
    cbar.set_ticklabels(['Low', 'Moderate', 'High'], fontsize=12)
    plt.title('{} GFS: LES Parameter for High Snowfall Rates | {} | FH: {}'.format(ds[time_dim][0].dt.strftime('%H00 UTC').item(), ds[time_dim][i].dt.strftime('%Y-%m-%d %H00 UTC').item(), count*3), fontsize=14)
    plt.savefig('plots/models/gfs/les/les_{}.png'.format(i), dpi=450, bbox_inches='tight')
    count += 1
