import xarray as xr
import metpy
import datetime as datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from metpy.plots import USCOUNTIES
from scipy.ndimage import gaussian_filter

start_date = datetime.datetime.now().replace(hour=12, minute=0, second=0, microsecond=0)

url = "https://thredds.ucar.edu/thredds/dodsC/grib/NCEP/NDFD/NWS/CONUS/NOAAPORT/Best"
ds = xr.open_dataset(url, engine='netcdf4')
ds = ds.metpy.parse_cf()
ds_latlon = ds.metpy.assign_latitude_longitude()

hi_temp = ds_latlon['Maximum_temperature_height_above_ground_12_Hour_Maximum'].sel(height_above_ground1=2)

possible_time_dims = ['time', 'time1', 'time2', 'time3']

time_dim = None
for dim in possible_time_dims:
    if dim in hi_temp.dims:
        time_dim = dim
        break
if time_dim is None:
    print("No valid time dimension found in hi_temp.")

start_date_formatted = np.datetime64(start_date)

hi_temp = hi_temp.sel(**{time_dim: slice(start_date_formatted, None)})

if time_dim:
    for i in range(min(9, len(hi_temp[time_dim]))):
        datetime_str = np.datetime_as_string(hi_temp[time_dim].values[i], unit='m')
        temp_f = (hi_temp[i,:,:] - 273.15) * 9/5 + 32
        temp_f_smoothed = gaussian_filter(temp_f, sigma=2)
        fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
        ax.set_extent([-87.5, -82, 41, 46.5])
        ax.add_feature(cfeature.STATES, linewidth=0.5)
        ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)
        plt.contourf(ds_latlon['longitude'], ds_latlon['latitude'], temp_f_smoothed, cmap='jet', levels=np.arange(30, 111, 1))
        cbar = plt.colorbar(label='Temperature (F)', fraction=0.046, pad=0.04)
        cbar.set_ticks(np.arange(30, 111, 10))
        cbar.set_ticklabels(['30', '40', '50', '60', '70', '80', '90', '100', '110'])
        contour = plt.contour(ds_latlon['longitude'], ds_latlon['latitude'], temp_f_smoothed, colors='black', levels=np.arange(30, 111, 1), linewidths=0.5)
        plt.clabel(contour, inline=True, fontsize=8, fmt='%1.0f')
        plt.title('NDFD High Temperature for {}'.format(datetime_str))
        plt.savefig('plots/temps/high/NDFD_High_Temp_{}.png'.format(i), dpi=450, bbox_inches='tight')
