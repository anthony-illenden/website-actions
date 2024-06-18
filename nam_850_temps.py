import numpy as np
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import xarray as xr
from scipy.ndimage import gaussian_filter

tds_nam = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/NAM/CONUS_12km/latest.html')
nam_ds = tds_nam.datasets[0]
ds = xr.open_dataset(nam_ds.access_urls['OPENDAP'])
ds = ds.metpy.parse_cf()
ds =  ds.metpy.assign_latitude_longitude()

possible_time_dims = ['time', 'time1', 'time2', 'time3']

time_dim = None
for dim in possible_time_dims:
    if dim in ds.dims:
        time_dim = dim
        break
if time_dim is None:
    raise ValueError('Could not find the time dimension')

count = 0
for i in range(0, len(ds['time'])):
    temp_850mb = ds['Temperature_isobaric'].sel(isobaric=85000)
    temp_850mb = temp_850mb.isel(**{time_dim: i})
    temp_850mb_smoothed = gaussian_filter(temp_850mb, sigma=2) 
    gph_850mb = ds['Geopotential_height_isobaric'].sel(isobaric=85000)
    gph_850mb = gph_850mb.isel(**{time_dim: i})
    gph_850mb_smoothed = gaussian_filter(gph_850mb, sigma=2)
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-125, -66.9, 23, 49.4])
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    cf = plt.contourf(temp_850mb['longitude'], temp_850mb['latitude'], temp_850mb_smoothed-273.15, cmap='hsv_r', levels=np.arange(-20, 31, 1),extend ='both')
    cbar = plt.colorbar(cf, orientation='horizontal', label='Temperature (C)', fraction=0.046, pad=0.04)
    isohypses = plt.contour(gph_850mb['longitude'], gph_850mb['latitude'], gph_850mb / 10, colors='k', levels=np.arange(90, 181, 3))
    plt.clabel(isohypses, inline=True, fontsize=12, fmt='%1.0f')
    isotherms = plt.contour(temp_850mb['longitude'], temp_850mb['latitude'], temp_850mb_smoothed-273.15, colors='white', levels=np.arange(-20, 31, 1), linestyles='dashed', linewidths=0.50)
    a = ds[time_dim][0].dt.strftime('%H00 UTC').item()
    b = temp_850mb[time_dim].dt.strftime('%Y-%m-%d %H00 UTC').item()
    c = count
    plt.title(f'{a} NAM 12KM: 850 hPa Temperature and Geopotential Heights | {b} | FH: {c*3}')
    plt.savefig('plots/models/nam/850/temps/temps_{}.png'.format(i), dpi=450, bbox_inches='tight')
    count +=1
