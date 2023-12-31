# %%
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

# %%
# Edit these variables with needed information about data being processed
fn = '/home/maureenjcohen/lmd_data/standard_xtra.nc' # Path to file with model output
experiment_name = 'VenusTest' # For labelling new files
outputdir = '/home/maureenjcohen/misc_data/' + experiment_name
t_select = (0,5) # Range of times to be included
h_select = (0,-1) # Range of heights to be included
rho = 65 # Density of atmosphere in kg/m3
g_constant = 8.87 # Gravitational constant of planet in m/s2
# If your atmospheric density varies significantly within the model domain,
# you will have to get a density cube
heights = np.array([0.00, 0.03, 0.12, 0.32, 0.68, 1.23, 2.03, 3.10, 4.50, 6.23, 8.35,
               10.8, 13.7, 17.0, 20.7, 24.6, 28.3, 31.9, 35.2, 38.4, 41.4, 44.2,
               46.9, 49.5, 51.9, 54.1, 56.2, 58.1, 60.1, 61.9, 63.7, 65.5, 67.2,
               68.8, 70.5, 72.2, 73.8, 75.5, 77.1, 78.7, 80.2, 81.8, 83.3, 84.8,
               86.2, 87.8, 90.1, 92.9, 94.9, 101.])*1e3 
# Heights of Venus model output in m

### Functions for reorganising and reformatting LMD Planets simulation output
# %%
def make_file(ncout, step, udata, vdata, wdata, hghts, lats, lons, 
              time_len):
    """ Make an individual netCDF file from an empty Dataset"""
    # Create the dimensions of the new file, same as the old file
    ncout.createDimension('time', 1)
    ncout.createDimension('height', len(hghts))  
    ncout.createDimension('lat', len(lats))
    ncout.createDimension('lon', len(lons))

    # Create variable to store longitudes
    longitude = ncout.createVariable('Longitude', 'float32', ('lon',))
    longitude.units = 'degrees_east'
    longitude.axis = 'X'
    
    # Create variable to store latitudes
    latitude = ncout.createVariable('Latitude', 'float32', ('lat',))
    latitude.units = 'degrees_north'
    latitude.axis = 'Y'

    # Create variable to store heights (not pressure levels)
    height = ncout.createVariable('Height', 'float32', ('height',))
    height.units = 'm'
    height.axis = 'Z'
    height.positive = 'up'
    
    # Create variable to hold timestamps
    time = ncout.createVariable('Time', 'float32', ('time',))
    time.units = 'seconds since 1987-03-30 00:00:00'

    # Now create the variables that will hold your wind data
    # Note: W input data must be in m/s (model output is in Pa/s)
    uout = ncout.createVariable('U', 'float32', ('time', 'height', 'lat', 'lon'))
    uout.units = 'm/s'
    uout.interval_write = str(time_len)
    uout[:,:,:,:] = udata

    vout = ncout.createVariable('V', 'float32', ('time', 'height', 'lat', 'lon'))
    vout.units = 'm/s'
    vout.interval_write = str(time_len)
    vout[:,:,:,:] = vdata

    wout = ncout.createVariable('W', 'float32', ('time', 'height', 'lat', 'lon'))
    wout.units = 'm/s'
    wout.interval_write = str(time_len)
    wout[:,:,:,:] = wdata

    # Fill in the dimensions with the arrays from the original sim files
    latitude[:] = lats
    longitude[:] = lons
    height[:] = hghts

    # Now do some funky time stuff
    secs_passed = step*time_len # Number of secs passed since start of sim
    secs = timedelta(seconds=secs_passed)
    date = datetime(1987,3,30) + secs # Add time passed to start date
    time[:] = nc.date2num(date, time.units)
    print('File written for:', time[:], time.units)

 # %%
def extract_metadata(ncfile):
    """ Input a netcdf4 file and extract the metadata that will be used
    to create a new, reformatted netcdf4 file 
    
    Outputs: arrays of longitudes, latitudes, timestamps, and scalar value
             of the time interval between each output cube (in seconds)  """

    lons = ncfile['lon'][:]
    lats = ncfile['lat'][:]
    times = ncfile['time_counter'][:]
    t_interval = np.diff(ncfile['time_counter'][:])[0]

    return lons, lats, times, t_interval
# %%
def process_data(windcube, windunits, windtype):
    """ The LMD output data needs to be made compatible with Parcels requirements.
    
    Steps: 1. Parcels accepts Arakawa C-grids (LMD grid type), but the grid indexing
              is different. The data cube needs to be shifted to match what Parcels
              expects. We shift up one row and discard the pole to achieve this.
           2. Parcels accepts vertical wind in m/s, not pressure vertical velocity
              in Pa/s, so this cube needs to be converted (based on hydrostatic
              relationship)"""
    
    if windtype=='W':
        print('Processing upward wind')
        if windunits == 'Pa/s':
            print('Wind units are Pa/s, converting to m/s')
            w_wind = windcube[:]
            winddata = -1*w_wind/(rho*g_constant)
        elif windunits == 'm/s':
            print('Wind units are m/s, good to go')
        else:
            print('Cannot parse wind units')
    elif windtype=='V':
        print('Processing northward wind')
        winddata = windcube[:]
    elif windtype=='U':
        print('Processing eastward wind')
        winddata = windcube[:]
    else:
        print('Invalid wind type, must be U, V, or W')

    windout = np.roll(winddata, 1, axis=2) # Shift cube 1 row up Y-axis (latitudes)
    windout[:,:,0,:] = 0 # Fill the emptied latitude row with 0s 

    return windout
# %%
def selector(ncfile, inputtimes, inputheights, trange=(0,-1), hrange=(0,-1)):
    """ Function that selects subsets of the data to include in the Parcels
        input files
        e.g. a subset of the time range, or a subset of the height levels """
    
    ucube = ncfile['vitu'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    vcube = ncfile['vitv'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    wcube = ncfile['vitw'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]

    select_times = inputtimes[trange[0]:trange[-1]]
    select_heights = inputheights[hrange[0]:hrange[-1]]

    return ucube, vcube, wcube, select_times, select_heights

# %%
def run_preprocess(inputfile, savedir, testname):

    """ Input the path to a file containing LMD Planets simulation data
        Input directory into which to save preprocessed files for Parcels
        Output files preprocessed to be compatible with Parcels     """
    
    lons, lats, times, t_interval = extract_metadata(inputfile)

    ucube, vcube, wcube, selected_times, selected_heights = selector(inputfile, times, heights, trange=t_select, hrange=h_select)

    u_data = process_data(ucube, 'm/s', 'U')
    v_data = process_data(vcube, 'm/s', 'V')
    w_data = process_data(wcube, 'Pa/s', 'W')

    for i in range(0,len(selected_times)):
        ncout = nc.Dataset(savedir + '/' + testname + f'_{i}.nc',
                           'w', format='NETCDF4')
        make_file(ncout, i, u_data[i,:,:,:], v_data[i,:,:,:],
                   w_data[i,:,:,:], selected_heights, 
                  lats, lons, t_interval)
    ncout.close(); del ncout

# %%
if __name__ == "__main__":

    ds = nc.Dataset(fn)

    run_preprocess(ds, outputdir, experiment_name)
    

# %%
