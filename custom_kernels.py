""" Custom kernels for Parcels v4 trajectory simulations in planetary atmospheres """
# %% Imports
from parcels import Particle, Variable, StatusCode
import numpy as np

# %%
## Custom particle class for Venus simulations.
## u_conv is the dimensionless Ornstein-Uhlenbeck / AR(1) red-noise state used by
## the convection kernel. It starts at 0; the OU process relaxes to its stationary
## N(0,1) distribution within a few correlation times (~1 hr for tau = 20 min),
## which is negligible against the multi-day integration.
## stuck is set to 1 by boundary_stick when a particle is frozen at the boundary.
VenusParticle = Particle.add_variable([
                Variable("u_conv", dtype=np.float32, initial=0.0, to_write=True),
                Variable("stuck", dtype=np.int32, initial=0, to_write=True),
])

# %%
## Custom particle class for Mars simulations.
## stuck is set to 1 by boundary_stick when a particle is frozen at the boundary.
MarsParticle = Particle.add_variable(
                Variable("stuck", dtype=np.int32, initial=0, to_write=True),
)

# %%
## Custom particle class for Venus Aerobot simulations.
## w_bal is the vertical velocity of the balloon (m/s) and v_bal the current
## balloon volume (m^3), both starting at 0. u_conv and stuck are as in VenusParticle.
BalloonParticle = Particle.add_variable([
                Variable("w_bal", dtype=np.float32, initial=0.0, to_write=True),
                Variable("v_bal", dtype=np.float32, initial=0.0, to_write=True),
                Variable("u_conv", dtype=np.float32, initial=0.0, to_write=True),
                Variable("stuck", dtype=np.int32, initial=0, to_write=True),
])

# %%
def CheckOutOfBounds(particles, fieldset):
    state = np.asarray(particles.state)
    particles.state = np.where(state == StatusCode.ErrorOutOfBounds,
                               StatusCode.Delete, state)

# %%
def boundary_stick(particles, fieldset):
    """ Freeze particles that leave the domain (e.g. over the pole) at their
        last position instead of crashing the simulation. """
    oob = np.asarray(particles.state) == StatusCode.ErrorOutOfBounds
    particles.dx = np.where(oob, 0.0, np.asarray(particles.dx))
    particles.dy = np.where(oob, 0.0, np.asarray(particles.dy))
    particles.dz = np.where(oob, 0.0, np.asarray(particles.dz))
    particles.stuck = np.where(oob, 1, np.asarray(particles.stuck))
    particles.state = np.where(oob, StatusCode.Success, np.asarray(particles.state))

# %%
def periodicBC(particles, fieldset):
    """ Wrap the end-of-step longitude x + dx into [lon_min, lon_min + lon_span).
        Parcels applies dx after the whole kernel chain, so the wrap is done on dx;
        this kernel must come after every kernel that adds to dx (advection,
        diffusion) so that output longitudes are always inside the domain. """
    x = np.asarray(particles.x)
    x_new = np.mod(x + np.asarray(particles.dx) - fieldset.lon_min, fieldset.lon_span) + fieldset.lon_min
    particles.dx = x_new - x

# %%
def polarBC(particles, fieldset):
    lat = np.asarray(particles.y)
    lon = np.asarray(particles.x)
    lon_min, lon_span = fieldset.lon_min, fieldset.lon_span
    lat_max, lat_min = fieldset.lat_max, fieldset.lat_min

    north = lat > lat_max
    south = lat < lat_min

    lat = np.where(north, 180.0 - lat, lat)
    lat = np.where(south, -180.0 - lat, lat)
    lon = np.where(north | south, np.mod(lon - lon_min + 180.0, lon_span) + lon_min, lon)

    particles.y = lat
    particles.x = lon

# %%
def smagdiff(particles, fieldset):
    """ Smagorinsky horizontal diffusion scheme adapted from
        https://oceanparcels.org/tutorials/advanced/smagorinsky_diffusion.html
        Needs 'cell_areas' Field on the fieldset (area in m^2).
    """
    dx = 0.01
    # gradients are computed by using a local central difference.
    updx, vpdx = fieldset.UV[particles.t, particles.z, particles.y, particles.x + dx, particles]
    umdx, vmdx = fieldset.UV[particles.t, particles.z, particles.y, particles.x - dx, particles]
    updy, vpdy = fieldset.UV[particles.t, particles.z, particles.y + dx, particles.x, particles]
    umdy, vmdy = fieldset.UV[particles.t, particles.z, particles.y - dx, particles.x, particles]

    dudx = (updx - umdx) / (2 * dx)
    dudy = (updy - umdy) / (2 * dx)

    dvdx = (vpdx - vmdx) / (2 * dx)
    dvdy = (vpdy - vmdy) / (2 * dx)

    A = fieldset.cell_areas[particles.t, particles.z, particles.y, particles.x, particles]
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
    # trapezoidal envelope: 1 in [lo,hi], linear ramp over `edge` at each side, 0 outside
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
    """ Particles pushed below the lowest model level are put back on it. """
    hit = np.asarray(particles.state) == StatusCode.ErrorThroughSurface
    particles.dz = np.where(hit, fieldset.W.grid.depth[0] - np.asarray(particles.z),
                            np.asarray(particles.dz))
    particles.state = np.where(hit, StatusCode.Success, np.asarray(particles.state))

# %%
def balloon_vertical(particles, fieldset):
    """ Balloon vertical-wind kernel: updates the balloon's vertical velocity
        and displacement using the equations for net force, drag and buoyancy
        from the prototype Venus Aerobot paper by Izraelevitz et al:

            (m_total + m_virtual) * d2z/dt2 = rho_atm * V * g  -  m_total * g  -  F_drag (Eq. 1)
            F_drag = 0.5 * rho_atm * C_D * A_ref * v^2 * v_hat (Eq. 11)
            m_virtual = C_m * rho_atm * V (Eq. 13)

        In the outer step, rho_atm and w_atm are interpolated at the particle position,
        and the Ornstein-Uhlenbeck process from convection_ou is used to directly advance
        w_atm, the vertical atmospheric wind.

        In the inner loop, the balloon's volume and virtual mass are calculated at each
        substep along with the forces on the balloon. The balloon's relative velocity is advanced
        by a semi-exponential update similar to the OU process in convection_ou:

            a_buoy = (rho*V*fieldset.g_Venus - fieldset.m_total*fieldset.g_Venus) / (fieldset.m_total + m_virtual)
            w_eq = a_buoy * tau_vertical
            w_rel = w_eq + (w_rel_old - w_eq)*exp(-|dt_inner| / tau_vertical)

        Additionally, a local gradient is used to extrapolate rho_atm at each substep.
        Since convection_ou is restructured in this kernel, it should not be added to the
        kernel list in run_venus_simulation.py when using BalloonParticle.

        All particles are updated together as arrays; the substep loop runs over time,
        not over particles.
        """
    n_inner = 60
    dt = np.asarray(particles.dt)
    dt_inner = dt / n_inner
    z = np.asarray(particles.z)
    rho0 = np.asarray(fieldset.RHO[particles.t, particles.z, particles.y, particles.x, particles])
    w_atm = np.asarray(fieldset.W[particles.t, particles.z, particles.y, particles.x, particles])
    rho_atm = rho0

    # Compute density gradient:
    rho1 = np.asarray(fieldset.RHO[particles.t, z + 200, particles.y, particles.x, particles])
    rho2 = np.asarray(fieldset.RHO[particles.t, z - 200, particles.y, particles.x, particles])
    slope = (rho1 - rho2) / 400.0 # d_rho/dz

    # Convective OU noise (as in convection_ou), added directly to w_atm:
    lo, hi, edge = fieldset.conv_z_lo, fieldset.conv_z_hi, fieldset.conv_edge
    env = np.clip(np.minimum((z - (lo - edge)) / edge, ((hi + edge) - z) / edge), 0.0, 1.0)
    active = env > 0.0
    a = np.exp(-np.fabs(dt) / fieldset.conv_tau)
    g = np.random.normal(0.0, 1.0, size=z.shape)
    u_conv = np.asarray(particles.u_conv)
    u_conv = np.where(active, a * u_conv + np.sqrt(1.0 - a * a) * g, u_conv)
    particles.u_conv = u_conv
    w_atm = w_atm + fieldset.conv_sigma * env * u_conv

    # Inner sub-step loop:
    displacement = np.zeros_like(rho0)
    w_bal = np.asarray(particles.w_bal, dtype=np.float64)
    for _ in range(n_inner):
        # Compute volume & virtual mass:
        Vol = np.minimum(10.86 * fieldset.m_gas_ZP / rho_atm, fieldset.V_infl) # Displaced volume [m^3], capped at full inflation
        m_virtual = fieldset.C_m * rho_atm * Vol # Apparent extra mass [kg] from Eq. (13)

        # Compute forces:
        w_rel_old = w_bal - w_atm # Relative velocity [m/s]
        tau_vertical = (fieldset.m_total + m_virtual) / (0.5 * rho_atm * fieldset.C_D_top * fieldset.A_top * np.fabs(w_rel_old)) # Drag relaxation [s]

        # Update relative velocity:
        a_buoy = (rho_atm*Vol*fieldset.g_Venus - fieldset.m_total*fieldset.g_Venus) / (fieldset.m_total + m_virtual) # Acceleration due to buoyancy [m/s^2]
        w_eq = a_buoy * tau_vertical # Equilibrium velocity [m/s]
        w_rel = w_eq + (w_rel_old - w_eq)*np.exp(-np.fabs(dt_inner) / tau_vertical) # Relative velocity [m/s]

        # Update vertical velocity, displacement and atmospheric density:
        w_bal = w_rel + w_atm # Balloon's vertical velocity [m/s]
        displacement += w_bal*dt_inner # Accumulates displacement [m]
        rho_atm = rho0 + slope*displacement # Extrapolates atmospheric density [kg/m^3]

    particles.w_bal = w_bal
    particles.v_bal = Vol # Volume at the start of the last substep
    # Update altitude position [m]:
    particles.dz += displacement
# %%
