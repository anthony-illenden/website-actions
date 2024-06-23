import xarray as xr
import metpy
import datetime as datetime
import matplotlib.pyplot as plt
import numpy as np
from siphon.catalog import TDSCatalog
from datetime import timedelta
import matplotlib.dates as mdates
import pandas as pd

start_date = datetime.datetime.now().replace(hour=13, minute=0, second=0, microsecond=0)
end_date = start_date + timedelta(days=7)

url = "https://thredds.ucar.edu/thredds/dodsC/grib/NCEP/NDFD/NWS/CONUS/NOAAPORT/Best"
ds = xr.open_dataset(url, engine='netcdf4')
ds = ds.metpy.parse_cf()
ds_latlon = ds.metpy.assign_latitude_longitude()

hi_temp = ds_latlon['Maximum_temperature_height_above_ground_12_Hour_Maximum'].sel(height_above_ground1=2)
low_temp = ds_latlon['Minimum_temperature_height_above_ground_12_Hour_Minimum'].sel(height_above_ground1=2)
dewp = ds_latlon['Dewpoint_temperature_height_above_ground'].sel(height_above_ground1=2)

def find_time_dim_and_select(data, start_date):
    possible_time_dims = ['time', 'time1', 'time2', 'time3']
    time_dim = None
    for dim in possible_time_dims:
        if dim in data.dims:
            time_dim = dim
            break
    if time_dim is None:
        print("No valid time dimension found.")
        return data 
    
    start_date_formatted = np.datetime64(start_date)
    return data.sel(**{time_dim: slice(start_date_formatted, end_date)}), time_dim

def fetch_data(ds_latlon, hi_temp, low_temp, dewp, target_lat, target_lon):
    latitudes = ds_latlon['latitude'].values
    longitudes = ds_latlon['longitude'].values
    distances = np.sqrt((latitudes - target_lat)**2 + (longitudes - target_lon)**2)
    min_distance_idx = np.unravel_index(np.argmin(distances), distances.shape)
    y_idx, x_idx = min_distance_idx
    high = hi_temp.isel(x=x_idx, y=y_idx)
    low = low_temp.isel(x=x_idx, y=y_idx)
    dew = dewp.isel(x=x_idx, y=y_idx)

    return high, low, dew

def plot_temperature_forecast(city_name, lat, lon):
    city_abbr = city_abbreviations.get(city_name, city_name)
    highs, time_dim_high = find_time_dim_and_select(hi_temp, start_date)
    lows, time_dim_low = find_time_dim_and_select(low_temp, start_date)
    dews, time_dim_dew = find_time_dim_and_select(dewp, start_date)
    high, low, dew = fetch_data(ds_latlon, highs, lows, dews, lat, lon)

    high_dates = (pd.Series(high[time_dim_high]) - timedelta(hours=1)).dt.normalize()
    low_dates = (pd.Series(low[time_dim_low]) - timedelta(hours=13)).dt.normalize()
  
    highs_bar = plt.bar(high_dates, ((high - 273.15) * 1.8 + 32), color='red', label='High Temp', zorder=2)
    lows_bars = plt.bar(low_dates, ((low - 273.15) * 1.8 + 32), color='blue', label='Low Temp', width=0.65, zorder=2)
    plt.title('7-Day NDFD Temperature Forecast for {}'.format(city_name))

    ax = plt.gca()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%A,\n%m-%d'))
    plt.ylim(40, 100)
    ax.set_yticks(np.arange(40, 101, 5))
    ax.yaxis.grid(True, which='major', linestyle='--', linewidth=0.5, zorder=0, alpha=1, color='black')
    plt.ylabel('Temperature (F)')
    for bar in highs_bar:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom')

    for bar in lows_bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2., height, f'{int(height)}', ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig('plots/temps/ndfd/cities/{}.png'.format(city_abbr), dpi=450)
    plt.close()

cities = [("Rochester Hills, MI", 42.66, -83.41), ("Royal Oak, MI", 42.49, -83.14), ("Ferndale, MI", 42.46, -83.13), 
          ("Livonia, MI", 42.37, -83.35), ("Detroit, MI", 42.33, -83.05), ("Grand Rapids, MI", 42.96, -85.66)]
city_abbreviations = {
    "Rochester Hills, MI": "RH",
    "Royal Oak, MI": "RO",
    "Ferndale, MI": "FRN",
    "Livonia, MI": "LIV",
    "Detroit, MI": "DET",
    "Grand Rapids, MI": "GRR",
}
for city_name, lat, lon in cities:
    plot_temperature_forecast(city_name, lat, lon)
