# %%
import xarray as xr, numpy as np
import os
from datetime import timedelta
from parcels import FieldSet, SphericalMesh, ParticleSet, ParticleFile
from parcels.kernels import AdvectionRK4
import custom_kernels
from fieldset_tools import add_cell_areas
import config  # per-run settings; copy config_example.py -> config.py
import importlib

# %%
importlib.reload(config) # Run if config settings have been changed
importlib.reload(custom_kernels)
from custom_kernels import *

# %%
PLANET_RADIUS = 6.0518e6 # Venus radius in m

# %%
def add_venus_context(fieldset):
    """ Add the run constants (diffusion, convection and Aerobot parameters)
        to the fieldset as Parcels v4 context values. """
    fieldset.add_context("Cs", 0.1)

    # Convective vertical-wind (OU/AR(1)) parameters, Vega-1 calibration.
    # See vega_convection.py for their derivation from the Vega balloon W_a data.
    fieldset.add_context("conv_sigma", 0.6)     # convective w std [m/s], de-biased W_a
    fieldset.add_context("conv_tau", 1200.0)    # correlation time [s] (~20 min)
    fieldset.add_context("conv_z_lo", 48000.0)  # convective layer bottom [m]
    fieldset.add_context("conv_z_hi", 55000.0)  # convective layer top [m]
    fieldset.add_context("conv_edge", 2000.0)   # taper half-width at each edge [m]

    # Aerobot parameters
    # Drag and virtual mass coefficients:
    fieldset.add_context("C_D_top", 0.8) # Drag coefficient for vertical motion
    fieldset.add_context("C_m", 0.2) # Virtual mass coefficient
    fieldset.add_context("A_top", 19.6) # Upper area [m^2]
    fieldset.add_context("V_infl", 72.6) # Enclosed volume of the profile [m^3]
    fieldset.add_context("g_Venus", 8.87) # Venus gravitational acceleration [m/s^2]
    fieldset.add_context("m_total", 62) # Total mass: helium + envelopes + payload [kg]
    fieldset.add_context("m_gas_ZP", 5.75) # Mass of helium in the ZP balloon [kg]
    #fieldset.add_context("m_gas_SP", ) # Mass of helium in the SP balloon [kg] - not needed until later

# %%
def check_seed_depths(fieldset, depths):
    """ Seed depths must lie inside the file's vertical axis.
        Catches unit mistakes (e.g. a Height axis accidentally written in km)
        up front, instead of as cryptic out-of-bounds errors mid-run. """
    depth_axis = np.asarray(fieldset.U.grid.depth)
    for d in depths:
        if not depth_axis.min() <= d <= depth_axis.max():
            raise ValueError(
                f"Particle seed depth {d} m is outside the input file's height axis "
                f"[{depth_axis.min():g}, {depth_axis.max():g}] m. Check PARTICLE_DEPTH "
                "in config.py and the Height axis units in the input files "
                "(Parcels expects metres)."
            )

# %%
def main():
    """ Build the FieldSet, seed the cloud particles, and run the simulation.

    Wrapped in a function (rather than running at import time) so that
    custom_kernels and helpers can be imported and tested without kicking off a
    multi-day integration. Run the simulation with:  python run_venus_simulation.py
    """
    datadir = config.DATA_DIR  # Where the preprocessed input file lives
    savedir = config.SAVE_DIR  # Where to save output parquet files

    input_path = os.path.join(datadir, config.EXPERIMENT_NAME + '.nc')
    ds = xr.open_dataset(input_path)

    fieldset = FieldSet.from_sgrid_conventions(ds, mesh=SphericalMesh(radius=PLANET_RADIUS))
    fieldset = add_cell_areas(fieldset, PLANET_RADIUS)
    # periodicBC wrap domain
    fieldset.add_context("lon_min", -180.0)
    fieldset.add_context("lon_span", 360.0)
    add_venus_context(fieldset)

    check_seed_depths(fieldset, config.PARTICLE_DEPTH)

    # Create particle set (initial positions come from config.py)
    pset_clouds = ParticleSet(
        fieldset=fieldset,
        pclass=BalloonParticle,
        x=config.PARTICLE_LON,
        y=config.PARTICLE_LAT,
        z=config.PARTICLE_DEPTH,)

    output_path = os.path.join(savedir, os.path.splitext(config.OUTPUT_NAME)[0] + '.parquet')
    output_file = ParticleFile(output_path, outputdt=timedelta(minutes=config.OUTPUT_MINUTES), mode='w')

    pset_clouds.execute([AdvectionRK4, balloon_vertical, smagdiff, periodicBC, boundary_stick],
                 runtime=timedelta(days=config.RUNTIME_DAYS),
                 dt=timedelta(minutes=config.DT_MINUTES),
                 output_file=output_file,
                 verbose_progress=True,
    )


# %%
if __name__ == "__main__":
    main()
