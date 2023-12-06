# %%
import numpy as np
import netCDF4 as nc

#%%
# Edit these variables with needed information about data being processed
fn = '/home/maureenjcohen/lmd_data/standard_xtra.nc' # Path to file with model output
time_len = 100870 # Seconds passed between outputs from the simulation (interval_write)
heights = [0.00, 0.03, 0.12, 0.32, 0.68, 1.23, 2.03, 3.10, 4.50, 6.23, 8.35,
               10.8, 13.7, 17.0, 20.7, 24.6, 28.3, 31.9, 35.2, 38.4, 41.4, 44.2,
               46.9, 49.5, 51.9, 54.1, 56.2, 58.1, 60.1, 61.9, 63.7, 65.5, 67.2,
               68.8, 70.5, 72.2, 73.8, 75.5, 77.1, 78.7, 80.2, 81.8, 83.3, 84.8,
               86.2, 87.8, 90.1, 92.9, 94.9, 101.] # Heights of model output in km

### Class to hold data which will be put into its own netCDF file


# %%
def make_file(ncout, times, hghts, lats, lons, windtype='U'):
    """ Make an individual netCDF file from an empty Dataset"""
    # define axis size
    ncout.createDimension('time', len(times))
    ncout.createDimension('height', len(hghts))  
    ncout.createDimension('lat', len(lats))
    ncout.createDimension('lon', len(lons))
#    ncout.createDimension('lev', len(levs))

    # create longitude axis
    longitude = ncout.createVariable('Longitude', 'float32', ('lon',))
    longitude.units = 'degrees_east'
    longitude.axis = 'X'
    
    # create latitude axis
    latitude = ncout.createVariable('Latitude', 'float32', ('lat',))
    latitude.units = 'degrees_north'
    latitude.axis = 'Y'

    # create height axis
    height = ncout.createVariable('Height', 'float32', ('height',))
    height.units = 'm'
    height.axis = 'Z'
    height.positive = 'up'
    
    # create time axis
    time = ncout.createVariable('Time', 'float32', ('time',))
    time.units = 'seconds since 0000-01-01 00:00:00'

    uout = ncout.createVariable(f'{windtype}', 'float32', ('time', 'height', 'lat', 'lon'))
    uout.units = 'm/s'
    uout.interval_write = str(time_len)
