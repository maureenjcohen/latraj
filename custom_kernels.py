""" Custom kernels for Parcels trajectory simulations in the Venus atmosphere """
# %% Imports
import parcels
from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, StatusCode, Variable
import math
import numpy as np

# %%
class VenusParticle(JITParticle):
    """ Custom particle class for Venus simulations """
    next_convection = Variable('next_convection', dtype=np.float32, 
                               initial=0.0, to_write=True)

# %%
def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()

# %%
def periodicBC(particle, fieldset, time):
    if particle.lon < fieldset.halo_west:
        particle_dlon += fieldset.halo_east - fieldset.halo_west
    elif particle.lon > fieldset.halo_east:
        particle_dlon -= fieldset.halo_east - fieldset.halo_west

# %%
def smagdiff(particle, fieldset, time):
    """ Smagorinsky horizontal diffusion scheme adapted from
        https://oceanparcels.org/tutorials/advanced/smagorinsky_diffusion.html
    """
    dx = 0.01
    # gradients are computed by using a local central difference.
    updx, vpdx = fieldset.UV[time, particle.depth, particle.lat, particle.lon + dx]
    umdx, vmdx = fieldset.UV[time, particle.depth, particle.lat, particle.lon - dx]
    updy, vpdy = fieldset.UV[time, particle.depth, particle.lat + dx, particle.lon]
    umdy, vmdy = fieldset.UV[time, particle.depth, particle.lat - dx, particle.lon]

    dudx = (updx - umdx) / (2 * dx)
    dudy = (updy - umdy) / (2 * dx)

    dvdx = (vpdx - vmdx) / (2 * dx)
    dvdy = (vpdy - vmdy) / (2 * dx)

    A = fieldset.cell_areas[time, 0, particle.lat, particle.lon]
    sq_deg_to_sq_m = (1852 * 60) ** 2 * math.cos(particle.lat * math.pi / 180)
    A = A / sq_deg_to_sq_m
    Kh = fieldset.Cs * A * math.sqrt(dudx**2 + 0.5 * (dudy + dvdx) ** 2 + dvdy**2)

    dlat = parcels.ParcelsRandom.normalvariate(0.0, 1.0) * math.sqrt(
        2 * math.fabs(particle.dt) * Kh
    )
    dlon = parcels.ParcelsRandom.normalvariate(0.0, 1.0) * math.sqrt(
        2 * math.fabs(particle.dt) * Kh
    )

    particle_dlat += dlat
    particle_dlon += dlon

# %%
def convection(particle, fieldset, time):
    """ Simple convection kernel that moves particles up/down by a distance Z
        drawn from a normal distribution with mean 0 and stddev 500 m every X
        seconds, with X also drawn from a normal distribution with a mean of
        2 hours and a standard deviation of 1 hour.
    """
    if time >= particle.next_convection:
        if particle.depth <=55000.0 and particle.depth >= 48000.0: # convection layer
            Z = parcels.ParcelsRandom.normalvariate(0.0, 500.0) # in meters
            particle_ddepth = particle_ddepth + Z
            particle.next_convection = time + math.fabs(
                parcels.ParcelsRandom.normalvariate(2 * 3600, 1 * 3600))
# %%
