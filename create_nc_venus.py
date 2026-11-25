"""" Script to extract data and metadata from Venus PCM outputs, reformat them as needed by
    the Parcels v4 package, and write them to a new batch of netcdf files for trajectory analysis.
    
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
import xarray as xr
import importlib
import parcels._sgrid as sgrid
from fieldset_tools import agrid_sgrid_metadata, add_zonal_halo
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
rho = 65                # Density of atmosphere in kg/m3 (for Pa/s -> m/s vertical wind)
g_constant = 8.87   # Gravitational constant of planet in m/s2
# If your atmospheric density varies significantly within the model domain,
# you will have to get a density cube.
# The model level heights below are a fixed property of the Venus PCM output,
# not a per-run setting, so they stay here rather than in config.py. VPCM runs
# come on either a 50-level or a 78-level vertical grid; both sets are kept (in
# METRES -- venuslab/venusdata.py keeps the same values in km) and the right one
# is chosen per input file by its level count, mirroring venusdata.py.

heights50 = np.array([9.4504710e-03, 4.8338074e-02, 1.5300426e-01, 3.7324557e-01,
       7.5411582e-01, 1.3392169e+00, 2.1695750e+00, 3.2827427e+00,
       4.7116103e+00, 6.4832373e+00, 8.6177559e+00, 1.1128335e+01,
       1.4020479e+01, 1.7288918e+01, 2.0912828e+01, 2.4761955e+01,
       2.8558361e+01, 3.2130924e+01, 3.5487930e+01, 3.8643879e+01,
       4.1625511e+01, 4.4466003e+01, 4.7189983e+01, 4.9781677e+01,
       5.2210793e+01, 5.4474300e+01, 5.6580246e+01, 5.8551716e+01,
       6.0436192e+01, 6.2270138e+01, 6.4059181e+01, 6.5802261e+01,
       6.7509895e+01, 6.9198883e+01, 7.0857117e+01, 7.2468269e+01,
       7.4031113e+01, 7.5544807e+01, 7.7022598e+01, 7.8483002e+01,
       7.9931213e+01, 8.1360779e+01, 8.2763893e+01, 8.4131378e+01,
       8.5458458e+01, 8.6925591e+01, 8.8879539e+01, 9.1365417e+01,
       9.4282623e+01, 9.7128662e+01])*1e3
# Heights of 50-level Venus model output in m

heights78 = np.array([9.45306290e-03, 4.83510271e-02, 1.53046221e-01, 3.73352647e-01,
       7.54337728e-01, 1.33960199e+00, 2.17016840e+00, 3.28359985e+00,
       4.71279430e+00, 6.48478937e+00, 8.61961460e+00, 1.11300869e+01,
       1.40210838e+01, 1.72868519e+01, 2.09060841e+01, 2.47483273e+01,
       2.85367489e+01, 3.21012650e+01, 3.54501762e+01, 3.85974350e+01,
       4.15684090e+01, 4.43950195e+01, 4.71000137e+01, 4.96712570e+01,
       5.20832405e+01, 5.43336220e+01, 5.64306679e+01, 5.83991661e+01,
       6.02850800e+01, 6.21196480e+01, 6.39083824e+01, 6.56521988e+01,
       6.73622131e+01, 6.90550003e+01, 7.07164383e+01, 7.23295288e+01,
       7.38933105e+01, 7.54046097e+01, 7.68734818e+01, 7.83174973e+01,
       7.97435608e+01, 8.11459122e+01, 8.25166473e+01, 8.38473129e+01,
       8.51307755e+01, 8.65320282e+01, 8.83726730e+01, 9.07071533e+01,
       9.34497375e+01, 9.65266037e+01, 9.90068893e+01, 1.00908829e+02,
       1.02824753e+02, 1.04765404e+02, 1.06736778e+02, 1.08729904e+02,
       1.10737305e+02, 1.12745575e+02, 1.14689415e+02, 1.16525398e+02,
       1.18274094e+02, 1.19964310e+02, 1.21610733e+02, 1.23218185e+02,
       1.24797585e+02, 1.26358299e+02, 1.27900208e+02, 1.29420120e+02,
       1.30924683e+02, 1.32434402e+02, 1.33977646e+02, 1.35579971e+02,
       1.37256210e+02, 1.39007202e+02, 1.40822159e+02, 1.42683762e+02,
       1.44574203e+02, 1.46479889e+02])*1e3
# Heights of 78-level Venus model output in m

LEVEL_HEIGHTS = {50: heights50, 78: heights78}

def heights_for(n_levels):
    """ Return the model level heights [m] for a VPCM run with n_levels
        vertical levels. Raises rather than guessing for an unknown grid. """
    try:
        return LEVEL_HEIGHTS[n_levels]
    except KeyError:
        raise ValueError(
            f"Input file has {n_levels} vertical levels; known VPCM grids: "
            f"{sorted(LEVEL_HEIGHTS)}. Add its level heights to LEVEL_HEIGHTS "
            "in create_nc_venus.py."
        ) from None

# Scalar (non-wind) fields to carry through to the Parcels input file, mapping
# the output name Parcels will see to the VPCM source variable name and units.
# Each is written only if present in the input file, so older wind-only outputs
# still preprocess without modification.
SCALAR_FIELDS = {
    'TEMP': ('temp', 'K'),
    'PRES': ('pres', 'Pa'),
    'RHO':  ('rho',  'kg/m3'),
}

### Functions for reorganising and reformatting LMD Planets simulation output
# %%
def build_dataset(udata, vdata, wdata, hghts, lats, lons,
                  n_times, time_len, n_cells, scalars=None):
    """ Create an xarray Dataset holding the full run (all timesteps) with the
        SGRID metadata needed by Parcels v4's FieldSet.from_sgrid_conventions.
        udata/vdata/wdata are the full 4D (time, height, lat, lon) cubes.

        scalars is an optional {name: (data, units)} dict of extra 4D fields
        (temperature, pressure, density) written on the same dims as the winds.

        Coord names must equal their dims, except the vertical coordinate, which
        Parcels v4 requires to be named depth even though it indexes the height
        dim. Timestamps are evenly spaced by time_len seconds from the start date
        (index-based, so the series always begins at zero regardless of where
        t_select starts). A zonal halo of n_cells is added on each side so that
        particles can be wrapped in longitude. """
    base = np.datetime64('1987-03-30T00:00:00')
    times = base + (np.arange(n_times) * time_len).astype('int64').astype('timedelta64[s]')

    data_vars = {
        "U": (["time","height","lat","lon"], np.asarray(udata, dtype="float32"), {"units": "m/s"}),
        "V": (["time","height","lat","lon"], np.asarray(vdata, dtype="float32"), {"units": "m/s"}),
        "W": (["time","height","lat","lon"], np.asarray(wdata, dtype="float32"), {"units": "m/s"}),
    }
    # Extra scalar fields, on the same dims and the same reversed lon ordering as
    # the winds so Parcels reads them onto an identical grid.
    if scalars is None:
        scalars = {}
    for name, (data, units) in scalars.items():
        data_vars[name] = (["time","height","lat","lon"], np.asarray(data, dtype="float32"), {"units": units})

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "lon": ("lon", np.asarray(lons, dtype="float32"), {"axis": "X", "units": "degrees_east"}),
            "lat": ("lat", np.asarray(lats, dtype="float32"), {"axis": "Y", "units": "degrees_north"}),
            "depth": ("height", np.asarray(hghts, dtype="float32"), {"axis": "Z", "units": "m", "positive": "up"}),
            "time": ("time", times, {"axis": "T"})
        }
    )
    ds = ds.sortby("lat")
    ds = add_zonal_halo(ds, n_cells)
    ds = sgrid._attach_sgrid_metadata(ds, agrid_sgrid_metadata())
    return ds

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
            raise ValueError(f"Cannot parse wind units {windunits!r}; expected 'Pa/s' or 'm/s'")
    elif windtype=='V':
        print('Processing northward wind')
        winddata = windcube[:,:,:,::-1]
    elif windtype=='U':
        print('Processing eastward wind')
        winddata = -windcube[:,:,:,::-1]
    else:
        raise ValueError(f"Invalid wind type {windtype!r}, must be U, V, or W")

    return winddata
# %%
def process_scalar(datacube, fieldname):
    """ Reformat a scalar (non-wind) field for Parcels.

    Scalars get the same longitude reversal as the winds (see module docstring
    note 1) but no sign flip. Reversing the axis is a reindexing of the data;
    only vector components pointing along that axis change sign, which is why
    process_data negates U but nothing here is negated. """

    print('Processing scalar field', fieldname)

    return datacube[:,:,:,::-1]
# %%
def selector(ncfile, inputtimes, inputheights, trange=(0,None), hrange=(0,None)):
    """ Function that selects subsets of the data to include in the Parcels
        input files
        e.g. a subset of the time range, or a subset of the height levels

        Scalar fields in SCALAR_FIELDS are picked up when the input file has
        them and skipped when it does not, so wind-only outputs still work. """

    ucube = ncfile['vitu'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    vcube = ncfile['vitv'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
    wcube = ncfile['vitwz'][trange[0]:trange[1],hrange[0]:hrange[1],:,:]

    scalar_cubes = {}
    for name, (srcname, _units) in SCALAR_FIELDS.items():
        if srcname in ncfile.variables:
            scalar_cubes[name] = ncfile[srcname][trange[0]:trange[1],hrange[0]:hrange[1],:,:]
        else:
            print('Field', srcname, 'not found in input file, skipping')

    select_times = inputtimes[trange[0]:trange[1]]
    select_heights = inputheights[hrange[0]:hrange[1]]

    return ucube, vcube, wcube, scalar_cubes, select_times, select_heights

# %%
def run_preprocess(inputfile, savedir, testname, halo_cells=3):

    """ Input the path to a file containing LMD Planets simulation data
        Input directory into which to save preprocessed files for Parcels
        Output files preprocessed to be compatible with Parcels     """
    
    lons, lats, times, t_interval = extract_metadata(inputfile)

    # Pick the 50- or 78-level height array to match this file's vertical grid
    # (mirrors venuslab/venusdata.py, which keys off the presnivs length).
    heights = heights_for(inputfile['vitu'].shape[1])

    ucube, vcube, wcube, scalar_cubes, selected_times, selected_heights = selector(inputfile, times, heights, trange=t_select, hrange=h_select)

    u_data = process_data(ucube, 'm/s', 'U')
    v_data = process_data(vcube, 'm/s', 'V')
    w_data = process_data(wcube, 'm/s', 'W')

    # Named scalar_data, not rho/temp/pres: `rho` at module scope is the scalar
    # background density from config.py that process_data uses for the Pa/s
    # conversion, and must not be shadowed here.
    scalar_data = {name: (process_scalar(cube, name), SCALAR_FIELDS[name][1])
                   for name, cube in scalar_cubes.items()}

    # Write every timestep into a single netCDF file (one file per run).
    ds = build_dataset(u_data, v_data, w_data, selected_heights,
                       lats, lons, len(selected_times), t_interval, halo_cells,
                       scalars=scalar_data)
    outpath = savedir + '/' + testname + '.nc'
    ds.to_netcdf(outpath,
                 encoding={"time": {"units": "seconds since 1987-03-30 00:00:00"}})
    print('File written: ', outpath, '| times', str(ds.time.values[0]),
          '->', str(ds.time.values[-1]))

# %%
if __name__ == "__main__":

    ds = nc.Dataset(fn)

    run_preprocess(ds, outputdir, experiment_name)
    

# %%
