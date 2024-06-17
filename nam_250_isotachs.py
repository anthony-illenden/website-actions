import numpy as np
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import xarray as xr
from scipy.ndimage import gaussian_filter
from metpy.units import units

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

for i in range(0, len(ds['time'])):
    u_250mb = ds['u-component_of_wind_isobaric'].sel(isobaric=25000)
    v_250mb = ds['v-component_of_wind_isobaric'].sel(isobaric=25000)
    u_250mb = u_250mb.isel(**{time_dim: i})
    v_250mb = v_250mb.isel(**{time_dim: i})
    wspd = mpcalc.wind_speed(u_250mb, v_250mb)*1.94384
    div = mpcalc.divergence(u_250mb, v_250mb)
    wspd_smoothed = gaussian_filter(wspd, sigma=2)
    div_smoothed = gaussian_filter(div, sigma=4)

    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-125, -66.9, 23, 49.4])
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    cf = plt.contourf(wspd['longitude'], wspd['latitude'], wspd_smoothed, cmap='BuPu', levels=np.arange(60, 201, 20),extend ='both')
    cbar = plt.colorbar(cf, orientation='horizontal', label='Wind Speed (knots)', fraction=0.046, pad=0.04)
    isotachs = plt.contour(wspd['longitude'], wspd['latitude'], wspd_smoothed, colors='k', levels=np.arange(60, 201, 20))
    div_contour = plt.contour(wspd['longitude'], wspd['latitude'], div_smoothed*1e5, colors='magenta', levels=np.arange(-12, 0, 2))
    plt.clabel(div_contour, inline=True, fontsize=12, fmt='%1.0f')
    plt.title('NAM 12KM: 250 hPa Isotachs and Divergence {}'.format(u_250mb[time_dim].dt.strftime('%Y-%m-%d %H00 UTC').item()))
    plt.savefig('plots/models/nam/250/isotachs/isotachs_{}.png'.format(i), dpi=450, bbox_inches='tight')
    plt.show()
