"""" Script to extract data and metadata from the Mars PCM, OU version outputs, reformat them as needed by
    the Parcels V4 package, and write them to a new batch of netcdf files for trajectory analysis.
    
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
import xarray as xr
import importlib
import parcels._sgrid as sgrid
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
# The Mars PCM 'level' dimension holds sigma coordinates (pressure / surface
# pressure). We convert them to log-pressure heights once, using the standard
# scale-height relation h = -H * ln(sigma) with H = 11 km. The sigma values
# below were read from a control output (level dimension) and hardcoded.
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
        winddata = windcube[...] 
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
def add_zonal_halo(ds, n_cells):
    """ Add halo in longitude so that advection can wrap around.
        Parcels4 currently doesn't wrap.
        Hopefully temporary and v4 won't need this in its final release. """
    lon = ds.lon.values
    west = ds.isel(lon=slice(-n_cells, None)).assign_coords(lon=lon[-n_cells:] - 360.0)
    east = ds.isel(lon=slice(0, n_cells)).assign_coords(lon=lon[:n_cells] + 360.0)
    return xr.concat([west, ds, east], dim="lon")

# %%
def build_dataset(udata, vdata, wdata, hghts, lats, lons, n_times, time_len, n_cells):
    """ Create xarray Dataset from winds and coords. Attach SGRID metadata for v4's 
        FieldSet.from_sgrid_conventions.
        
        Coord names must be equal to their dims. Vertical coord must be named depth 
        (hardcoded in v4) even if it indexes the height dim.                    """
    
    base = np.datetime64('1987-03-30T00:00:00')
    times = base + (np.arange(n_times) * time_len).astype('int64').astype('timedelta64[s]')

    ds = xr.Dataset(
        data_vars={
            "U": (["time","height","lat","lon"], np.asarray(udata, dtype="float32"), {"units": "m/s"}),
            "V": (["time","height","lat","lon"], np.asarray(vdata, dtype="float32"), {"units": "m/s"}),
            "W": (["time","height","lat","lon"], np.asarray(wdata, dtype="float32"), {"units": "m/s"}),
        },
        coords={
            "lon": ("lon", np.asarray(lons, dtype="float32"), {"axis": "X", "units": "degrees_east"}),
            "lat": ("lat", np.asarray(lats, dtype="float32"), {"axis": "Y", "units": "degrees_north"}),
            "depth": ("height", np.asarray(hghts, dtype="float32"), {"axis": "Z", "units": "m", "positive":"up"}),
            "time": ("time", times, {"axis": "T"})
        }
    )
    ds = ds.sortby("lat")
    ds = add_zonal_halo(ds, n_cells)
    meta = sgrid.SGrid2DMetadata(
        cf_role="grid_topology",
        topology_dimension=2,
        node_dimensions=("lon", "lat"),
        face_dimensions=(
            sgrid.FaceNodePadding("lon_c", "lon", sgrid.Padding.LOW),
            sgrid.FaceNodePadding("lat_c", "lat", sgrid.Padding.LOW),
            ),
        vertical_dimensions=(sgrid.FaceNodePadding("height_c", "height", sgrid.Padding.LOW),),    
        )
    ds = sgrid._attach_sgrid_metadata(ds, meta) # Parcels4 helper function that attached metadata to grid
    return ds 

# %%
def run_preprocess(inputfile, savedir, testname, halo_cells):

    """ Input the path to a file containing LMD Planets simulation data
        Input directory into which to save preprocessed files for Parcels
        Output files preprocessed to be compatible with Parcels     """
    
    lons, lats, times, t_interval = extract_metadata(inputfile)

    ucube, vcube, wcube, selected_times, selected_heights = selector(inputfile, times, heights, trange=t_select, hrange=h_select)

    u_data = process_data(ucube, 'm/s', 'U')
    v_data = process_data(vcube, 'm/s', 'V')
    w_data = process_data(wcube, 'm/s', 'W')

    ds = build_dataset(u_data, v_data, w_data, selected_heights,
                       lats, lons, len(selected_times), t_interval, halo_cells)
    
    outpath = savedir + '/' + testname + '.nc'
    ds.to_netcdf(outpath,
                 encoding={"time": {"units": "seconds since 1987-03-30 00:00:00"}})
    print('File written: ', outpath, '| times', str(ds.time.values[0]),
        '->', str(ds.time.values[-1]))

# %%
if __name__ == "__main__":

    ds = nc.Dataset(fn)

    run_preprocess(inputfile=ds, savedir=outputdir, testname=experiment_name, halo_cells=3)
    

# %%
