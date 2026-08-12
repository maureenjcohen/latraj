""" Custom kernels for Parcels trajectory simulations in the Venus atmosphere """
# %% Imports
import parcels
from parcels import FieldSet, ParticleSet, JITParticle, ScipyParticle, StatusCode, Variable
from datetime import timedelta
import math
import numpy as np
from mpmath import mpf

# %%
class VenusParticle(JITParticle):
    """ Custom particle class for Venus simulations.

    Carries u_conv, the dimensionless Ornstein-Uhlenbeck / AR(1) red-noise state
    used by the convection kernel. It starts at 0; the OU process relaxes to its
    stationary N(0,1) distribution within a few correlation times (~1 hr for
    tau = 20 min), which is negligible against the multi-day integration.
    """
    u_conv = Variable('u_conv', dtype=np.float32, initial=0.0, to_write=True)
    stuck = Variable('stuck', dtype=np.int32, initial=0.0, to_write=True)

# %%
class BalloonParticle(JITParticle):
    """ Custom particle class for Aerobot simulations.

    Carries w_bal, the vertical velocity of the balloon (m/s), starting at 0;
    and v_bal, the current balloon volume (m^3), starting at ....
    """
    w_bal = Variable('w_bal', dtype=np.float32, initial=0.0, to_write=True)
    v_bal = Variable('v_bal', dtype=np.float32, initial=0.0, to_write=True)
    u_conv = Variable('u_conv', dtype=np.float32, initial=0.0, to_write=True)
    stuck = Variable('stuck', dtype=np.int32, initial=0.0, to_write=True)

# %%
def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()

# %%
def boundary_stick(particle, fieldset, time):
    """ Freeze particles that leave the domain (e.g. over the pole) at their
        last in-bounds position instead of crashing the simulation. """
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle_dlon = 0.0
        particle_dlat = 0.0
        particle_ddepth = 0.0
        particle.stuck = 1
        particle.state = StatusCode.Success

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
def convection_ou(particle, fieldset, time):
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
    z = particle.depth
    env = 0.0
    if z > fieldset.conv_z_lo - fieldset.conv_edge and z < fieldset.conv_z_hi + fieldset.conv_edge:
        if z < fieldset.conv_z_lo:
            env = (z - (fieldset.conv_z_lo - fieldset.conv_edge)) / fieldset.conv_edge
        elif z > fieldset.conv_z_hi:
            env = ((fieldset.conv_z_hi + fieldset.conv_edge) - z) / fieldset.conv_edge
        else:
            env = 1.0

    if env > 0.0:
        a = math.exp(-math.fabs(particle.dt) / fieldset.conv_tau)
        g = parcels.ParcelsRandom.normalvariate(0.0, 1.0)
        particle.u_conv = a * particle.u_conv + math.sqrt(1.0 - a * a) * g
        particle_ddepth += fieldset.conv_sigma * env * particle.u_conv * particle.dt
        
# %%
def surface_bounce(particle, fieldset, time):
     if particle.state == StatusCode.ErrorThroughSurface:
          particle_ddepth = 0.0
          particle.state = StatusCode.Success

# %%
def balloon_vertical(particle, fieldset, time):
    """ write a detailed description here like in example kernels"""
    displacement, i = 0.0, 0.0
    dt_inner = particle.dt / 30 # Sub-timesteps [s]
    rho_atm = fieldset.RHO[time, particle.depth, particle.lat, particle.lon]
    w_atm = fieldset.W[time, particle.depth, particle.lat, particle.lon]

    rho1 = fieldset.RHO[time, particle.depth + 1000, particle.lat, particle.lon] 
    rho2 = fieldset.RHO[time, particle.depth - 1000, particle.lat, particle.lon]     
    z1 = particle.depth + 1000
    z2 = particle.depth - 1000
    slope = (rho1 - rho2)/(z1 - z2)
    #print(slope)
    # Restructure convection_ou to add directly to w_atm:
    z = particle.depth 
    env = 0.0
    if z > fieldset.conv_z_lo - fieldset.conv_edge and z < fieldset.conv_z_hi + fieldset.conv_edge:
        if z < fieldset.conv_z_lo:
            env = (z - (fieldset.conv_z_lo - fieldset.conv_edge)) / fieldset.conv_edge
        elif z > fieldset.conv_z_hi:
            env = ((fieldset.conv_z_hi + fieldset.conv_edge) - z) / fieldset.conv_edge
        else:
            env = 1.0

    if env > 0.0:
        a = math.exp(-math.fabs(particle.dt) / fieldset.conv_tau)
        g = parcels.ParcelsRandom.normalvariate(0.0, 1.0)
        particle.u_conv = a * particle.u_conv + math.sqrt(1.0 - a * a) * g
        w_atm += fieldset.conv_sigma * env * particle.u_conv 
    
    while i < 30:
        i += 1
        # Compute displaced volume & virtual mass:
        Vol = 10.86 * fieldset.m_gas_ZP / rho_atm # Displaced volume [m^3] uses Vol instead of V because of JITParticle errors
        if Vol > fieldset.V_infl:
            Vol = fieldset.V_infl # Caps volume at maximum inflation 
        m_virtual = fieldset.C_m * rho_atm * Vol # Apparent extra mass [kg]
        
        w_rel_old = particle.w_bal - w_atm
        tau_vertical = (fieldset.m_total + m_virtual) / (0.5 * rho_atm * fieldset.C_D_top * fieldset.A_top * math.fabs(w_rel_old)) # Drag relaxation [s]
        F_drag = 0.5 * rho_atm * fieldset.C_D_top * fieldset.A_top * w_rel_old * math.fabs(w_rel_old) # Drag force [N]
        F_net = rho_atm*Vol*fieldset.g_Venus - fieldset.m_total*fieldset.g_Venus - F_drag # Net force from Eq. (1) [N]
        
        # Update vertical velocities and compute displacement:
        a_buoy = (rho_atm*Vol*fieldset.g_Venus - fieldset.m_total*fieldset.g_Venus) / (fieldset.m_total + m_virtual) # Acceleration due to buoyancy [m/s^2]
        w_eq = a_buoy * tau_vertical # Equilibrium velocity [m/s]
        w_rel = w_eq + (w_rel_old - w_eq)*math.exp(-math.fabs(dt_inner) / tau_vertical) # Relative velocity [m/s]
        
        particle.w_bal = w_rel + w_atm # Update vertical velocity [m/s]
        displacement += particle.w_bal*dt_inner # Update displacement [m]  
        particle_ddepth += displacement # Add displacement to particle position 
        rho_atm = math.fabs(slope*particle.depth)
        #print(rho_atm)
        #w_term = math.sqrt(math.fabs((2*(rho_atm*Vol*fieldset.g_Venus - fieldset.m_total*fieldset.g_Venus))/(rho_atm * fieldset.C_D_top * fieldset.A_top) ))
        #print(w_term)
        