"""" Script to extract data and metadata from the Mars PCM, OU version outputs, reformat them as needed by
    the Parcels package, and write them to a new batch of netcdf files for trajectory analysis.
    
    Notes FOR MARS:
        1. Time axis is in Martian sols (88775 seconds each).
        2. Since Mars has seasons (unlike Venus) and a long year (compared to Earth), a long time sample of at least one Mars year is necessary to capture dynamical changes over the year.

    Notes GENERAL:
        1. Parcels is mostly used for ocean data, so the vertical coordinate is depth in meters.
        2. Parcels only parses certain date formats. This script rewrites the time coordinate to use
        a format the package accepts.
        """
# %%
import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
import importlib
import config  # per-run settings; copy config_example.py -> config.py

# %%
importlib.reload(config) # Run if config settings have been changed

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
# The model level heights below are a fixed property of the Mars PCM output,
# not a per-run setting, so they stay here rather than in config.py.
# The Mars PCM 'level' dimension holds sigma coordinates (pressure / surface
# pressure). We convert them to log-pressure heights once, using the standard
# scale-height relation h = -H * ln(sigma) with H = 11 km. The sigma values
# below were read from a control output (level dimension) and hardcoded, the
# same way the Venus heights are hardcoded above in create_nc_venus.py.
scale_height = 11e3  # Mars atmospheric scale height H in m
MARS_SOL_S = 88775.0  # length of a Mars sol in seconds (PCM 'time' is in sols)
sigma = np.array([9.99500036e-01, 9.98395443e-01, 9.97060180e-01, 9.95446920e-01,
       9.93498921e-01, 9.91148591e-01, 9.88315403e-01, 9.84903872e-01,
       9.80801344e-01, 9.75875735e-01, 9.69973147e-01, 9.62916076e-01,
       9.54501629e-01, 9.44501340e-01, 9.32662129e-01, 9.18709874e-01,
       9.02355611e-01, 8.83306265e-01, 8.61279786e-01, 8.36025596e-01,
       8.07349801e-01, 7.75144577e-01, 7.39418328e-01, 7.00322866e-01,
       6.58173323e-01, 6.13454401e-01, 5.66809058e-01, 5.19008815e-01,
       4.70905840e-01, 4.23373938e-01, 3.77246320e-01, 3.33259106e-01,
       2.92008787e-01, 2.53927916e-01, 2.19279855e-01, 1.88169613e-01,
       1.60566464e-01, 1.36332795e-01, 1.15254290e-01, 9.70681310e-02,
       8.14869031e-02, 6.82174116e-02, 5.69743924e-02, 4.74896356e-02,
       3.95174474e-02, 3.28372344e-02, 2.72540580e-02, 2.25978158e-02,
       1.87215824e-02, 1.54994903e-02, 1.28244320e-02, 1.06057525e-02,
       8.76705907e-03, 7.24420371e-03, 5.98346628e-03, 4.93995240e-03,
       4.07619169e-03, 3.36092920e-03, 2.76808185e-03, 2.27584876e-03,
       1.86595216e-03, 1.52299611e-03, 1.23394514e-03, 9.87736043e-04,
       7.75077438e-04, 5.88519149e-04, 4.22874815e-04, 2.75911880e-04,
       1.48790888e-04, 4.52757049e-05])
heights = -scale_height * np.log(sigma)
# Heights of Mars model output in m (log-pressure, increasing = up)

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

    lons = ncfile['longitude'][:]
    lats = ncfile['latitude'][:]
    # The Mars PCM 'time' coordinate is in sols, but the reformatted Parcels
    # files (and the simulation dt/runtime) work in seconds, so convert here.
    times = ncfile['time'][:] * MARS_SOL_S
    t_interval = np.diff(ncfile['time'][:])[0] * MARS_SOL_S

    return lons, lats, times, t_interval
# %%
def process_data(windcube, windunits, windtype):
    """ The LMD output data needs to be made compatible with Parcels requirements.
    
    Steps:  1. Extract the relevant wind cube
            2. If the vertical wind is in Pa/s, convert to m/s
                                                                        """
    
    if windtype=='W':
        print('Processing upward wind')
        if windunits == 'Pa/s' or windunits == "Pa s-1":
            print('Wind units are Pa/s, converting to m/s')
            w_wind = windcube[...]
            winddata = -1*w_wind/(rho*g_constant)
        elif windunits == 'm/s' or windunits == "m s-1":
            print('Wind units are m/s, good to go')
            winddata = windcube[...]
        else:
            print('Cannot parse wind units')
    elif windtype=='V':
        print('Processing northward wind')
        winddata = windcube[...]
    elif windtype=='U':
        print('Processing eastward wind')
        winddata = windcube[...]  # Mars is prograde: no U negation, no lon flip
    else:
        print('Invalid wind type, must be U, V, or W')

    return winddata
# %%
def selector(ncfile, inputtimes, inputheights, trange=(0,None), hrange=(0,None)):
    """ Function that selects subsets of the data to include in the Parcels
        input files
        e.g. a subset of the time range, or a subset of the height levels """
    
    ucube = ncfile['U'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    vcube = ncfile['V'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    wcube = ncfile['W'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]

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
