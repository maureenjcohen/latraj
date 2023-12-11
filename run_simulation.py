# %%
import netCDF4 as nc
from parcels import FieldSet, ParticleSet, Variable, JITParticle, AdvectionRK4, plotTrajectoriesFile

# Define some variables for the run
datadir = '/home/maureenjcohen/misc_data/VenusTest/' # Where my data is located

# Set up Parcels inputs
filenames = {'U': datadir + 'VenusTest_U.nc',
             'V': datadir + 'VenusTest_V.nc',
             'W': datadir + 'VenusTest_W.nc'}

variables = {'U': 'U',
             'V': 'V',
             'W': 'W'}

dimensions = {'time': 'Time',
              'height': 'Height',
              'lat': 'Latitude',
              'lon': 'Longitude'}

# Create the FieldSet
fieldset = FieldSet.from_netcdf(filenames,
                                variables,
                                dimensions,
                                allow_time_extrapolation=True)

fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
fieldset.add_periodic_halo(zonal=True)

