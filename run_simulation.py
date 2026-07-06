# %%
import netCDF4 as nc
import xarray as xr
import os, re
from datetime import timedelta
import parcels
from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, AdvectionRK4_3D
from custom_kernels import periodicBC, convection_ou, smagdiff, VenusParticle

# Define any variables for the run
# %%
datadir = '/exomars/projects/mc5526/lagrangian_trajectory/VenusTest/' # Where my data is located
savedir = '/exomars/projects/mc5526/lagrangian_trajectory/outputs/' # Where to save output zarrs

# Get filepaths of all files in datadir
# %%
def alphanumeric_sort(lst):

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(lst, key=alphanum_key)

# %%
paths = []
for f in alphanumeric_sort(os.listdir(datadir)):
    paths.append(datadir+f)

# Set up Parcels inputs
# %%
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
# %%
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

# Create particle set
# %%
pset_clouds = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=VenusParticle,
    lon=[90., 90.],
    lat=[-0.0, 41.25],
    depth=[49000.0, 53000.0],)

# %%
pset_polar = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=JITParticle,
    lon=[10., -90.],
    lat=[65.0, 65.0],
    depth=[35400.0, 35400.0],)

# %%
output_file = pset_clouds.ParticleFile(
              name=savedir+'CloudTest2.zarr',
              outputdt=timedelta(minutes=30),
)

pset_clouds.execute([AdvectionRK4_3D, smagdiff, convection_ou, periodicBC],
             runtime=timedelta(days=60),
             dt=timedelta(minutes=5),
             output_file=output_file,
             verbose_progress=True, 
)


# %%
