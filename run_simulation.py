# %%
import netCDF4 as nc
import os
from datetime import timedelta
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D, AdvectionRK4, StatusCode

# Define any variables for the run
# %%
datadir = '/exomars/projects/mc5526/lagrangian_trajectory/VenusTest/' # Where my data is located

# Get filepaths of all files in datadir
# %%
paths = []
for f in sorted(os.listdir(datadir)):
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

# Create particle set
# %%
pset = ParticleSet.from_list(
    fieldset=fieldset,
    pclass=JITParticle,
    lon=[90.],
    lat=[0.0],
    depth=[72300.0],)

# Define kernels
# %%
def DeleteErrorParticle(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()

# %%
output_file = pset.ParticleFile(
              name='VenusTest.zarr',
              outputdt=timedelta(hours=1),
)

pset.execute([AdvectionRK4_3D],
             runtime=timedelta(days=5),
             dt=timedelta(minutes=30),
             output_file=output_file,
)

# %%
