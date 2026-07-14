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
    print(paths)
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

    cell_areas = parcels.Field(
        name="cell_areas", data=fieldset.U.cell_areas(), lon=x, lat=y)
    fieldset.add_field(cell_areas)
    fieldset.add_constant("Cs", 0.1)

    # Convective vertical-wind (OU/AR(1)) parameters, Vega-1 calibration.
    # See vega_convection.py for their derivation from the Vega balloon W_a data.
    fieldset.add_constant("conv_sigma", 0.6)     # convective w std [m/s], de-biased W_a
    fieldset.add_constant("conv_tau", 1200.0)    # correlation time [s] (~20 min)
    fieldset.add_constant("conv_z_lo", 48000.0)  # convective layer bottom [m]
    fieldset.add_constant("conv_z_hi", 55000.0)  # convective layer top [m]
    fieldset.add_constant("conv_edge", 2000.0)   # taper half-width at each edge [m]

    # Create particle set (initial positions come from config.py)
    pset_clouds = ParticleSet.from_list(
        fieldset=fieldset,
        pclass=VenusParticle,
        lon=config.PARTICLE_LON,
        lat=config.PARTICLE_LAT,
        depth=config.PARTICLE_DEPTH,)

    output_file = pset_clouds.ParticleFile(
                  name=savedir + config.OUTPUT_NAME,
                  outputdt=timedelta(minutes=config.OUTPUT_MINUTES),
    )

    pset_clouds.execute([AdvectionRK4_3D, smagdiff, convection_ou, periodicBC],
                 runtime=timedelta(days=config.RUNTIME_DAYS),
                 dt=timedelta(minutes=config.DT_MINUTES),
                 output_file=output_file,
                 verbose_progress=True,
    )


# %%
if __name__ == "__main__":
    main()
