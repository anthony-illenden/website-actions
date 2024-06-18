from siphon.catalog import TDSCatalog
import xarray as xr
from datetime import datetime, timedelta
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import numpy as np
import metpy.units as units
from metpy.plots import USCOUNTIES
from scipy.ndimage import gaussian_filter

dt = datetime.utcnow()
rtma_cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RTMA/CONUS_2p5km/catalog.xml')
rtma_data = rtma_cat.datasets['Latest Collection for Real Time Mesoscale Analysis 2.5 km'].remote_access(use_xarray=True)
rtma_data = rtma_data.metpy.parse_cf()

t = rtma_data['Temperature_Analysis_height_above_ground'].sel(time=dt, method='nearest').squeeze()
tf = t.metpy.convert_units('degF')
tf_smoothed = gaussian_filter(tf, sigma=3) 

fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
isotherms_shaded = ax.contourf(t.metpy.x, t.metpy.y, tf_smoothed, transform=t.metpy.cartopy_crs, cmap='jet', levels=np.arange(30, 111, 1))
cbar = fig.colorbar(isotherms_shaded, label='Degrees Fahrenheit', fraction=0.046, pad=0.04)
isotherms_outlined = ax.contour(t.metpy.x, t.metpy.y, tf_smoothed, transform=t.metpy.cartopy_crs, colors='k', levels=np.arange(30, 111, 1), linewidths=0.5)
plt.clabel(isotherms_outlined, inline=1, fontsize=12, fmt="%i")
cbar.set_ticks(np.arange(30, 111, 10))
cbar.set_ticklabels(['30', '40', '50', '60', '70', '80', '90', '100', '110'])
ax.set_extent([-91, -81, 40.5, 47.75])
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25, alpha=0.25)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.title('RTMA Air Temperature Analysis {} UTC'.format(tf.reftime.dt.strftime("%Y-%m-%d %H:%M").item()))
plt.savefig('plots/rtma/temp/latest_temp.png')
plt.show()

td = rtma_data['Dewpoint_temperature_Analysis_height_above_ground'].sel(time=dt, method='nearest').squeeze()
tdf = td.metpy.convert_units('degF')
tdf_smoothed = gaussian_filter(tdf, sigma=2) 

fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
isodrosotherms_shaded = ax.contourf(tdf.metpy.x, tdf.metpy.y, tdf_smoothed, transform=tdf.metpy.cartopy_crs, cmap='YlGn', levels=np.arange(45, 86, 1), extend='both')
fig.colorbar(isodrosotherms_shaded, label='Degrees Fahrenheit', fraction=0.046, pad=0.04)
isodrosotherms_outlined = ax.contour(tdf.metpy.x, tdf.metpy.y, tdf_smoothed, transform=tdf.metpy.cartopy_crs, colors='Black', levels=np.arange(45, 86, 1))
plt.clabel(isodrosotherms_outlined, inline=1, fontsize=12, fmt="%i")
ax.set_extent([-91, -81, 40.5, 47.75])
ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.25, alpha=0.25)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.5)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.title('RTMA Dewpoint Temperature Analysis {} UTC'.format(tdf.reftime.dt.strftime("%Y-%m-%d %H:%M").item()))
plt.savefig('plots/rtma/dewp/latest_dewp.png')
plt.show()
