http://localhost:8888/lab?token=60a6d0017b83416cd86410ef6595266c96e96f494f321e09""" Per-run and machine-specific configuration for the latraj pipeline.

    HOW TO USE THIS FILE
    --------------------
    Copy it to config.py and edit the values for your machine and experiment:

        cp config_example.py config.py

    config.py is listed in .gitignore, so your local paths and particle choices
    stay OUT of version control. That means they never appear in your pull
    requests and never conflict with anyone else's config. This template
    (config_example.py) IS tracked, so it documents every setting the scripts
    expect — if you add a new knob, add it here too.

    create_nc.py and run_simulation.py both read their settings from config.py.
    Physics constants that rarely change (the model level heights, the Vega
    convection calibration) live in the scripts themselves, not here.
"""

# --- Machine-specific paths -------------------------------------------------
# create_nc.py: the raw Venus PCM output to preprocess, a label for this
# experiment, and where to write the Parcels-ready netCDF files.
INPUT_FILE = '/exomars/data/internal/working/mc5526/VPCM_age_of_air/aoa35_96x96x50/Xins_141to145.nc'
EXPERIMENT_NAME = 'VenusTest'
OUTPUT_DIR = '/exomars/projects/mc5526/lagrangian_trajectory/' + EXPERIMENT_NAME

# run_simulation.py: where the preprocessed input files live, and where to
# write the output trajectory zarr.
DATA_DIR = '/exomars/projects/mc5526/lagrangian_trajectory/VenusTest/'
SAVE_DIR = '/exomars/projects/mc5526/lagrangian_trajectory/outputs/'
OUTPUT_NAME = 'CloudTest2.zarr'   # written into SAVE_DIR

# --- Preprocessing selection (create_nc.py) ---------------------------------
T_SELECT = (0, 20)      # range of timesteps to include (start, stop)
H_SELECT = (0, None)    # range of height levels to include (start, stop)
RHO = 65                # atmospheric density [kg/m3], for Pa/s -> m/s vertical-wind conversion
G_CONSTANT = 8.87       # planetary surface gravity [m/s2]

# --- Simulation run settings (run_simulation.py) ----------------------------
# Initial particle positions: one entry per particle, so all three lists must
# be the same length. This is the block you will most often edit.
PARTICLE_LON = [90., 90.]           # degrees east, -180..180
PARTICLE_LAT = [40, 40]        # degrees north
PARTICLE_DEPTH = [43000.0, 43000.0]  # altitude in metres (positive = up)

RUNTIME_DAYS = 60       # total simulated time
DT_MINUTES = 5          # integration timestep
OUTPUT_MINUTES = 30     # how often particle positions are written out
