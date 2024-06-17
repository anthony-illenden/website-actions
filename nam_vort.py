import numpy as np
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import xarray as xr
from scipy.ndimage import gaussian_filter
from matplotlib.colors import TwoSlopeNorm

tds_nam = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/NAM/CONUS_12km/latest.html')
nam_ds = tds_nam.datasets[0]
ds = xr.open_dataset(nam_ds.access_urls['OPENDAP'])
ds = ds.metpy.parse_cf()
ds =  ds.metpy.assign_latitude_longitude()

norm = TwoSlopeNorm(vmin=-10, vcenter=0, vmax=60)

possible_time_dims = ['time', 'time1', 'time2', 'time3']

time_dim = None
for dim in possible_time_dims:
    if dim in ds.dims:
        time_dim = dim
        break
if time_dim is None:
    raise ValueError('Could not find the time dimension')

for i in range(0, len(ds['time'])):
    abs_vort = ds['Absolute_vorticity_isobaric'].sel(isobaric1=50000)
    abs_vort = abs_vort.isel(**{time_dim: i})
    gph_500 = ds['Geopotential_height_isobaric'].sel(isobaric=50000)
    gph_500 = gph_500.isel(**{time_dim: i})
    abs_vort_smoothed = gaussian_filter(abs_vort, sigma=2) 
    gph_500_smoothed = gaussian_filter(gph_500, sigma=2) 
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-125, -66.9, 23, 49.4])
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    cf = plt.contourf(abs_vort['longitude'], abs_vort['latitude'], abs_vort_smoothed*1e5, cmap='RdGy_r', levels=np.arange(-10, 61, 5), norm=norm,extend ='both')
    cbar = plt.colorbar(cf, orientation='horizontal', label='Vorticity (1/s)', fraction=0.046, pad=0.04)
    isohypses = plt.contour(gph_500['longitude'], gph_500['latitude'], gph_500_smoothed / 10, colors='k', levels=np.arange(480, 620, 4))
    plt.clabel(isohypses, inline=True, fontsize=12, fmt='%1.0f')
    plt.title('NAM 12KM: 500 hPa Absolute Vorticity and Geopotential Height {}'.format(abs_vort['time'].dt.strftime('%Y-%m-%d %H UTC').item()))
    plt.savefig('plots/models/nam/vort/vort_{}.png'.format(i), dpi=450, bbox_inches='tight')
    #plt.show()
