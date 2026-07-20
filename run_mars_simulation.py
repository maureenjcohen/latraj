# %%
import netCDF4 as nc
import xarray as xr
import os, re
from datetime import timedelta
import parcels
from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, AdvectionRK4_3D
import custom_kernels
import config  # per-run settings; copy config_example.py -> config.py
import importlib

# %%
importlib.reload(config) # Run if config settings have been changed
importlib.reload(custom_kernels)
from custom_kernels import *

# Get filepaths of all files in datadir
# %%
def alphanumeric_sort(lst):

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(lst, key=alphanum_key)


# %%
def main():
    """ Build the FieldSet, seed the cloud particles, and run the simulation.

    Wrapped in a function (rather than running at import time) so that
    custom_kernels and helpers can be imported and tested without kicking off a
    multi-day integration. Run the simulation with:  python run_simulation.py
    """
    datadir = config.DATA_DIR  # Where the preprocessed input files live
    savedir = config.SAVE_DIR  # Where to save output zarrs

    paths = [datadir + '/' + f for f in alphanumeric_sort(os.listdir(datadir))]

    # Set up Parcels inputs
    filenames = {'U': paths,
                 'V': paths,
                 'W': paths}

    variables = {'U': 'U',
                 'V': 'V',
                 'W': 'W'}

    dimensions = {'time': 'Time',
                  'depth': 'Height',
                  'lat': 'Latitude',
                  'lon': 'Longitude'}

    # Create the FieldSet
    fieldset = FieldSet.from_netcdf(filenames,
                                    variables,
                                    dimensions,
                                    allow_time_extrapolation=False)

    fieldset.add_constant("halo_west", fieldset.U.grid.lon[0])
    fieldset.add_constant("halo_east", fieldset.U.grid.lon[-1])
    fieldset.add_periodic_halo(zonal=True, meridional=False)
    x = fieldset.U.grid.lon
    y = fieldset.U.grid.lat

    # Create particle set (initial positions come from config.py)
    pset_clouds = ParticleSet.from_list(
        fieldset=fieldset,
        pclass=JITParticle,
        lon=config.PARTICLE_LON,
        lat=config.PARTICLE_LAT,
        depth=config.PARTICLE_DEPTH,)

    output_file = pset_clouds.ParticleFile(
                  name=savedir + config.OUTPUT_NAME,
                  outputdt=timedelta(minutes=config.OUTPUT_MINUTES),
    )

    pset_clouds.execute([AdvectionRK4_3D, surface_bounce, periodicBC],
                 runtime=timedelta(days=config.RUNTIME_DAYS),
                 dt=timedelta(minutes=config.DT_MINUTES),
                 output_file=output_file,
                 verbose_progress=True,
    )


# %%
if __name__ == "__main__":
    main()

# %%
