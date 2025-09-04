# %%
import netCDF4 as nc
import xarray as xr
import os, re
from datetime import timedelta
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D, StatusCode

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
#fieldset.add_constant("halo_south", fieldset.U.grid.lat[0])
#fieldset.add_constant("halo_north", fieldset.U.grid.lat[-1])
fieldset.add_periodic_halo(zonal=True, meridional=False)

# Create particle set
# %%
pset = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=JITParticle,
    lon=[90., 90.],
    lat=[-0.0, 41.25],
    depth=[72300.0, 72300.0],)

# %%
pset_polar = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=JITParticle,
    lon=[10., -90.],
    lat=[65.0, 65.0],
    depth=[35400.0, 35400.0],)

# Define kernels
# %%
def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()

def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle_dlon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle_dlon -= fieldset.halo_east - fieldset.halo_west

# %%
output_file = pset_polar.ParticleFile(
              name=savedir+'PolarTest_long.zarr',
              outputdt=timedelta(hours=12),
)

pset_polar.execute([AdvectionRK4_3D, periodicBC],
             runtime=timedelta(days=60),
             dt=timedelta(minutes=30),
             output_file=output_file,
             verbose_progress=True, 
)


# %%
