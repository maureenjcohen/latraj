""" Custom kernels for Parcels v4 trajectory simulations in the Venus atmosphere """
# %% Imports
from parcels import Particle, Variable, StatusCode
import numpy as np

# %%
## Custom particle class for Venus simulations.
VenusParticle = Particle.add_variable(
                Variable("u_conv", dtype=np.float32, initial=0.0, to_write=True)
)

# %%
def CheckOutOfBounds(particles, fieldset):
    state = np.asarray(particles.state)
    particles.state = np.where(state == StatusCode.ErrorOutOfBounds, 
                               StatusCode.Delete, state)

# %%
def periodicBC(particles, fieldset):
    lon_min = fieldset.lon_min
    particles.x = np.mod(np.asarray(particles.x) - lon_min, fieldset.lon_span) + lon_min

# %%
def smagdiff(particles, fieldset):
    """ Smagorinsky horizontal diffusion scheme adapted from
        https://oceanparcels.org/tutorials/advanced/smagorinsky_diffusion.html
        Now adapted for Parcels4
        Needs 'cell_areas' Field on the fieldset (area in m^2).
    """
    dx = 0.01
    # gradients are computed by using a local central difference.
    updx, vpdx = fieldset.UV[particles.time, particles.z, particles.y, particles.x + dx, particles]
    umdx, vmdx = fieldset.UV[particles.time, particles.z, particles.y, particles.x - dx, particles]
    updy, vpdy = fieldset.UV[particles.time, particles.z, particles.y + dx, particles.x, particles]
    umdy, vmdy = fieldset.UV[particles.time, particles.z, particles.y - dx, particles.x, particles]

    dudx = (updx - umdx) / (2 * dx)
    dudy = (updy - umdy) / (2 * dx)

    dvdx = (vpdx - vmdx) / (2 * dx)
    dvdy = (vpdy - vmdy) / (2 * dx)

    A = fieldset.cell_areas[particles.time, particles.z, particles.y, particles.x, particles]
    deg2m = fieldset.UV.grid.deg2m
    sq_deg_to_sq_m = (deg2m ** 2) * np.cos(np.asarray(particles.y) * np.pi / 180)
    A = A / sq_deg_to_sq_m
    Kh = fieldset.Cs * A * np.sqrt(dudx**2 + 0.5 * (dudy + dvdx) ** 2 + dvdy**2)

    dlat = np.random.normal(0.0, 1.0, size=np.asarray(particles.y).shape) * np.sqrt(2 * np.fabs(particles.dt) * Kh)
    dlon = np.random.normal(0.0, 1.0, size=np.asarray(particles.y).shape) * np.sqrt(2 * np.fabs(particles.dt) * Kh)

    particles.dy += dlat
    particles.dx += dlon

# %%
def convection_ou(particles, fieldset):
    """ Convective vertical-wind kernel: an Ornstein-Uhlenbeck / AR(1) red-noise
        process calibrated to the Vega balloon anemometer data (see
        vega_convection.py). Each parcel carries a dimensionless, zero-mean,
        unit-variance state u_conv, advanced by the exact OU update

            a      = exp(-|dt| / tau)
            u_conv = a * u_conv + sqrt(1 - a^2) * g,   g ~ N(0, 1)

        so u_conv stays stationary N(0,1) with correlation time fieldset.conv_tau.
        The physical convective velocity applied to the parcel is

            w_conv = fieldset.conv_sigma * envelope(depth) * u_conv

        where envelope(z) is 1 inside [conv_z_lo, conv_z_hi], tapers linearly to 0
        over conv_edge at each edge, and is 0 outside. This is layered on top of
        the large-scale W from the advection kernel; it does not replace it. The
        update is exact for any dt, so no sub-stepping is needed.
    """
    z = np.asarray(particles.z)
    lo, hi, edge = fieldset.conv_z_lo, fieldset.conv_z_hi, fieldset.conv_edge
    # trapezoidal envelope: 1 in [lo,hi], linear rampge over `edge` at each side, 0 outside
    env = np.clip(np.minimum((z - (lo - edge)) / edge, ((hi + edge) - z)/ edge), 0.0, 1.0)
    active = env > 0.0

    a = np.exp(-np.fabs(particles.dt) / fieldset.conv_tau)
    g = np.random.normal(0.0, 1.0, size=z.shape)
    u_new = a * np.asarray(particles.u_conv) + np.sqrt(1.0 - a * a) * g
    particles.u_conv = np.where(active, u_new, np.asarray(particles.u_conv))
    particles.dz += np.where(active, 
                             fieldset.conv_sigma * env * np.asarray(particles.u_conv) * particles.dt, 
                             0.0)
# %%
def surface_bounce(particles, fieldset):
    hit = np.asarray(particles.state) == StatusCode.ErrorThroughSurface
    particles.dz = np.where(hit, fieldset.W.grid.depth[0] - np.asarray(particles.z), 
                            np.asarray(particles.dz))
    particles.state = np.where(hit, StatusCode.Evaluate, np.asarray(particles.state))
# %%
