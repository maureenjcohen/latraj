# %%
import xarray as xr, numpy as np
import os
from datetime import timedelta
import parcels._sgrid as sgrid
from parcels import FieldSet, SphericalMesh, ParticleSet, ParticleFile, Particle
from parcels.kernels import AdvectionRK4_3D
import custom_kernels_v4
import config  # per-run settings; copy config_example.py -> config.py
import importlib

# %%
importlib.reload(config) # Run if config settings have been changed
importlib.reload(custom_kernels_v4)
from custom_kernels_v4 import *

# %%
PLANET_RADIUS = 3.3895e6 # Mars radius in m

# %%
def add_cell_areas(fieldset, radius):
    """ Add cell_areas field to fieldset (needed for smagdiff)
        Computed directly from lat/lon and deg2m. Not needed
        if there is a pre-existing areas cube (like for VPCM). """
    lat = np.asarray(fieldset.U.grid.lat)
    lon = np.asarray(fieldset.U.grid.lon)
    xlon, ylat = np.meshgrid(lon, lat)
    dlat = np.deg2rad(np.abs(np.gradient(ylat, axis=0)))
    dlon = np.deg2rad(np.abs(np.gradient(xlon, axis=1)))
    dy = dlat*radius
    dx = dlon*radius*np.abs(np.cos(np.deg2rad(ylat)))
    area = (dy*dx).astype("float32")

    ca = xr.Dataset(
        {"cell_areas": (["time", "height","lat","lon"], area[None, None])},
        coords={
            "lon": ("lon", lon.astype("float32"), {"axis": "X"}),
            "lat": ("lat", lat.astype("float32"), {"axis": "Y"}),
            "depth": ("height", np.array([0.0], dtype="float32"), {"axis": "Z"}),
            "time": ("time", np.array([fieldset.time_interval.left]), {"axis":"T"})
        },
    )
    meta = sgrid.SGrid2DMetadata(
        cf_role="grid_topology", topology_dimension=2,
        node_dimensions=("lon", "lat"),
        face_dimensions=(sgrid.FaceNodePadding("lon_c", "lon", sgrid.Padding.LOW),
                         sgrid.FaceNodePadding("lat_c", "lat", sgrid.Padding.LOW)),
        vertical_dimensions=(sgrid.FaceNodePadding("height_c", "height", sgrid.Padding.LOW),))
    ca = sgrid._attach_sgrid_metadata(ca, meta)
    ca_fs = FieldSet.from_sgrid_conventions(ca, mesh=SphericalMesh(radius=radius))
    fieldset.add_field(ca_fs.cell_areas, name="cell_areas")





# %%
def main():
    """ Build the FieldSet, seed the cloud particles, and run the simulation.

    Wrapped in a function (rather than running at import time) so that
    custom_kernels and helpers can be imported and tested without kicking off a
    multi-day integration. Run the simulation with:  python run_mars_simulation_v4.py
    """
    datadir = config.DATA_DIR  # Where the preprocessed input files live
    savedir = config.SAVE_DIR  # Where to save output zarrs

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
        pclass=Particle,
        x=config.PARTICLE_LON,
        y=config.PARTICLE_LAT,
        z=config.PARTICLE_DEPTH,)

    output_path = os.path.join(savedir, os.path.splitext(config.OUTPUT_NAME)[0] + '.parquet')
    output_file = ParticleFile(output_path, outputdt=timedelta(minutes=config.OUTPUT_MINUTES), mode='w')

    pset.execute([periodicBC, AdvectionRK4_3D, surface_bounce],
                 runtime=timedelta(days=config.RUNTIME_DAYS),
                 dt=timedelta(minutes=config.DT_MINUTES),
                 output_file=output_file,
                 verbose_progress=True,
    )
    print("t    :", pset._data["t"][:])
    print("dt   :", pset._data["dt"][:])
    print("state:", pset._data["state"][:])
    print("x,y,z:", pset._data["x"][:], pset._data["y"][:], pset._data["z"][:])
    print("nparts:", len(pset))



# %%
if __name__ == "__main__":
    main()

# %%
