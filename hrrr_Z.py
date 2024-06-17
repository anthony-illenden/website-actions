import numpy as np
import metpy.calc as mpcalc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from siphon.catalog import TDSCatalog
import xarray as xr
from scipy.ndimage import gaussian_filter
import metpy
from metpy.plots import USCOUNTIES

tds_hrrr = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/HRRR/CONUS_2p5km/latest.html')
hrrr_ds = tds_hrrr.datasets[0]
ds = xr.open_dataset(hrrr_ds.access_urls['OPENDAP'])
ds = ds.metpy.parse_cf()
ds =  ds.metpy.assign_latitude_longitude()

possible_time_dims = ['time', 'time1', 'time2', 'time3']

time_dim = None
for dim in possible_time_dims:
    if dim in ds['Composite_reflectivity_entire_atmosphere'].dims:
        time_dim = dim
        break
if time_dim is None:
    raise ValueError('Could not find the time dimension')

for i in range(0, len(ds['time'])):
    dbz = ds['Composite_reflectivity_entire_atmosphere'].isel(**{time_dim: i})
    dbz_masked = np.ma.masked_array(dbz, dbz < 1)
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([-91, -81, 40.5, 47.75])
    ax.add_feature(cfeature.STATES, linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25)
    img = plt.pcolormesh(dbz['longitude'], dbz['latitude'], dbz_masked, cmap=metpy.plots.ctables.registry.get_colortable('NWSStormClearReflectivity'), vmin=-30, vmax=80)
    cbar = fig.colorbar(img, ax=ax, orientation='vertical', label='Reflectivity (dBZ)', fraction=0.046, pad=0.04)
    plt.title('{} HRRR: Composite Reflectivity | {} | FH: {}'.format(ds[time_dim][0].dt.strftime('%H00 UTC').item(), ds[time_dim][i].dt.strftime('%Y-%m-%d %H00 UTC').item(), i))
    plt.savefig('plots/models/hrrr/reflectivity/reflectivity_{}.png'.format(i), dpi=450, bbox_inches='tight')
    plt.show()
