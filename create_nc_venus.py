"""" Script to extract data and metadata from Venus PCM outputs, reformat them as needed by
    the Parcels package, and write them to a new batch of netcdf files for trajectory analysis.
    
    Notes FOR VENUS:
        1. VENUS ROTATES BACKWARDS. The VPCM outputs have a flipped longitude axis that runs from
        +180 to -180. The U-wind direction is defined as positive when flowing towards -180. 
        When plotting winds with matplotlib, the package will reverse the x-axis (longitude) coordinates 
        from small to large and flip the data being plotted with it. You can then multiple the u-wind
        by -1 to get the true orientation. HOWEVER, Parcels doesn't do any such reorientation. Keeping
        the longitude axis orientation as-is leads to OutOfBounds errors when running the trajectory 
        analysis. Accordingly, this script reverses the longitude axis to run from -180 to +180 and flips
        all data arrays along the lon axis to match. It also still multiplies the U-wind by -1 to get
        the desired orientation (where negative is towards -180).
        2. The VPCM runs on an Arakawa C grid with different grid indexing than expected by Parcels.
        However, the outputs have the U/V/W fields on the same coords, so apparently the PCM or XIOS
        regrids the data before writing to netcdf. Hence the data should be treated as an A grid.

    Notes GENERAL:
        1. Parcels is mostly used for ocean data, so the vertical coordinate is depth in meters.
        2. Parcels only parses certain date formats. This script rewrites the time coordinate to use
        a format the package accepts.
        """
# %%
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import config  # per-run settings; copy config_example.py -> config.py

# %%
# Run settings come from config.py (paths, time/height selection, planet constants).
fn = config.INPUT_FILE            # Path to file with model output
experiment_name = config.EXPERIMENT_NAME  # For labelling new files
outputdir = config.OUTPUT_DIR
t_select = config.T_SELECT        # Range of times to be included
h_select = config.H_SELECT        # Range of heights to be included
rho = config.RHO                  # Density of atmosphere in kg/m3 (for Pa/s -> m/s vertical wind)
g_constant = config.G_CONSTANT    # Gravitational constant of planet in m/s2
# If your atmospheric density varies significantly within the model domain,
# you will have to get a density cube.
# The model level heights below are a fixed property of the Venus PCM output,
# not a per-run setting, so they stay here rather than in config.py.
heights = np.array([0.,  0.05,  0.2,  0.4,  0.8,  1.3,  2.2,  3.3,  4.7,  6.5,  8.6,
       11.1, 14., 17.3, 20.9, 24.7, 28.5, 32.1, 35.4, 38.6, 41.6, 44.4,
       47.1, 49.7, 52.1, 54.3, 56.4, 58.4, 60.3, 62.1, 63.9, 65.6, 67.4,
       69., 70.7, 72.3, 73.9, 75.4, 76.9, 78.4, 79.8, 81.2, 82.6, 84.,
       85.3, 86.8, 88.7, 91.2, 94.1, 97.])*1e3 
# Heights of Venus model output in m

### Functions for reorganising and reformatting LMD Planets simulation output
# %%
def make_file(ncout, udata, vdata, wdata, hghts, lats, lons,
              n_times, time_len):
    """ Fill an empty Dataset with the full run: all timesteps written to a
        single netCDF file (time dimension = n_times) rather than one file per
        step. udata/vdata/wdata are the full 4D (time, height, lat, lon) cubes."""
    # Create the dimensions of the new file, same as the old file
    ncout.createDimension('time', n_times)
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
    # Note: W input data must be in m/s
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

    # Now do some funky time stuff. Build one timestamp per step, evenly spaced
    # by time_len seconds and counted from the start date (index-based, so the
    # series always begins at zero regardless of where t_select starts).
    dates = [datetime(1987, 3, 30) + timedelta(seconds=int(step * time_len))
             for step in range(n_times)]
    time[:] = nc.date2num(dates, time.units)
    print('File written for', n_times, 'timesteps:',
          time[0], '->', time[-1], time.units)

 # %%
def extract_metadata(ncfile):
    """ Input a netcdf4 file and extract the metadata that will be used
    to create a new, reformatted netcdf4 file 
    
    Outputs: arrays of longitudes, latitudes, timestamps, and scalar value
             of the time interval between each output cube (in seconds)  """

    lons = ncfile['lon'][::-1] # Reverse longitude axis to run -180 to +180
    lats = ncfile['lat'][:]
    times = ncfile['time_counter'][:]
    t_interval = np.diff(ncfile['time_counter'][:])[0]

    return lons, lats, times, t_interval
# %%
def process_data(windcube, windunits, windtype):
    """ The LMD output data needs to be made compatible with Parcels requirements.
    
    Steps:  1. Extract the relevant wind cube
            2. If the vertical wind is in Pa/s, convert to m/s
                                                                        """
    
    if windtype=='W':
        print('Processing upward wind')
        if windunits == 'Pa/s':
            print('Wind units are Pa/s, converting to m/s')
            w_wind = windcube[:,:,:,::-1]
            winddata = -1*w_wind/(rho*g_constant)
        elif windunits == 'm/s':
            print('Wind units are m/s, good to go')
            winddata = windcube[:,:,:,::-1]
        else:
            print('Cannot parse wind units')
    elif windtype=='V':
        print('Processing northward wind')
        winddata = windcube[:,:,:,::-1]
    elif windtype=='U':
        print('Processing eastward wind')
        winddata = -windcube[:,:,:,::-1]
    else:
        print('Invalid wind type, must be U, V, or W')

    return winddata
# %%
def selector(ncfile, inputtimes, inputheights, trange=(0,None), hrange=(0,None)):
    """ Function that selects subsets of the data to include in the Parcels
        input files
        e.g. a subset of the time range, or a subset of the height levels """
    
    ucube = ncfile['vitu'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    vcube = ncfile['vitv'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    wcube = ncfile['vitwz'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]

    select_times = inputtimes[trange[0]:trange[1]]
    select_heights = inputheights[hrange[0]:hrange[1]]

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
    w_data = process_data(wcube, 'm/s', 'W')

    # Write every timestep into a single netCDF file (one file per run),
    # rather than one file per timestep.
    ncout = nc.Dataset(savedir + '/' + testname + '.nc', 'w', format='NETCDF4')
    make_file(ncout, u_data, v_data, w_data, selected_heights,
              lats, lons, len(selected_times), t_interval)
    ncout.close(); del ncout

# %%
if __name__ == "__main__":

    ds = nc.Dataset(fn)

    run_preprocess(ds, outputdir, experiment_name)
    

# %%
