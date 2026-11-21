# %%
import xarray as xr, numpy as np
import os
from datetime import timedelta
from parcels import FieldSet, SphericalMesh, ParticleSet, ParticleFile
from parcels.kernels import AdvectionRK4_3D
import custom_kernels
import config  # per-run settings; copy config_example.py -> config.py
import importlib

# %%
importlib.reload(config) # Run if config settings have been changed
importlib.reload(custom_kernels)
from custom_kernels import *

# %%
PLANET_RADIUS = 3.3895e6 # Mars radius in m

# %%
def main():
    """ Build the FieldSet, seed the cloud particles, and run the simulation.

    Wrapped in a function (rather than running at import time) so that
    custom_kernels and helpers can be imported and tested without kicking off a
    multi-day integration. Run the simulation with:  python run_mars_simulation.py
    """
    datadir = config.DATA_DIR  # Where the preprocessed input files live
    savedir = config.SAVE_DIR  # Where to save output parquet files

    input_path = os.path.join(datadir, config.EXPERIMENT_NAME + '.nc')
    ds = xr.open_dataset(input_path)

    fieldset = FieldSet.from_sgrid_conventions(ds, mesh=SphericalMesh(radius = PLANET_RADIUS))
    # periodicBC wrap domain
    fieldset.add_context("lon_min", 0.0)
    fieldset.add_context("lon_span", 360.0)
    lat = np.asarray(fieldset.U.grid.lat)
    fieldset.add_context("lat_min", float(lat.min()))   # -87.159096
    fieldset.add_context("lat_max", float(lat.max()))   #  87.159096

    # Create particle set (initial positions come from config.py)
    pset = ParticleSet(
        fieldset=fieldset,
        pclass=MarsParticle,
        x=config.PARTICLE_LON,
        y=config.PARTICLE_LAT,
        z=config.PARTICLE_DEPTH,)

    output_path = os.path.join(savedir, os.path.splitext(config.OUTPUT_NAME)[0] + '.parquet')
    output_file = ParticleFile(output_path, outputdt=timedelta(minutes=config.OUTPUT_MINUTES), mode='w')

    pset.execute([AdvectionRK4_3D, surface_bounce, periodicBC, boundary_stick],
                 runtime=timedelta(days=config.RUNTIME_DAYS),
                 dt=timedelta(minutes=config.DT_MINUTES),
                 output_file=output_file,
                 verbose_progress=True,
    )



# %%
if __name__ == "__main__":
    main()

# %%
