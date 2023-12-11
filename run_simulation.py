# %%
import netCDF4 as nc
import os
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile

# Define any variables for the run
# %%
datadir = '/home/maureenjcohen/misc_data/VenusTest/' # Where my data is located

# Set up Parcels inputs
# %%
filenames = {'U': datadir + 'VenusTest_*.nc',
             'V': datadir + 'VenusTest_*.nc',
             'W': datadir + 'VenusTest_*.nc'}

#ds = nc.Dataset(datadir + os.listdir(datadir)[0])
#t_interval = float(ds['U'].interval_write)
#t_dim, h_dim = len(os.listdir(datadir)), ds['U'].shape[1]
#timestamps = [float(x*t_interval) for x in range(1,t_dim+1)]

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
                                allow_time_extrapolation=True)

fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
fieldset.add_periodic_halo(zonal=True)


# %%
