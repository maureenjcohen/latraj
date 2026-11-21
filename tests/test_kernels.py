import custom_kernels
import pytest
import numpy as np
import math
import xarray as xr
import parcels._sgrid as sgrid
from custom_kernels import (convection_ou, balloon_vertical, boundary_stick, periodicBC,
                            surface_bounce, smagdiff, MarsParticle, VenusParticle, BalloonParticle)
from fieldset_tools import agrid_sgrid_metadata, add_cell_areas, check_timestep
from parcels import FieldSet, ParticleSet, SphericalMesh
from parcels.kernels import AdvectionRK4, AdvectionRK4_3D

def toy_fieldset(fields, lon, lat, alt, n_times=1):
    """ Build a time-constant Parcels v4 FieldSet from 3-D (alt, lat, lon) arrays,
        using the same SGRID A-grid layout as the preprocessed model files.
        n_times > 1 repeats the data at hourly timestamps, which gives the
        FieldSet a finite time interval (needed by add_cell_areas). """
    data_vars = {name: (["time", "height", "lat", "lon"],
                        np.repeat(np.asarray(arr, dtype=np.float32)[None], n_times, axis=0))
                 for name, arr in fields.items()}
    times = np.datetime64("1987-03-30T00:00:00", "ns") + np.arange(n_times) * np.timedelta64(1, "h")
    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "lon": ("lon", lon, {"axis": "X", "units": "degrees_east"}),
            "lat": ("lat", lat, {"axis": "Y", "units": "degrees_north"}),
            "depth": ("height", alt, {"axis": "Z", "units": "m", "positive": "up"}),
            "time": ("time", times, {"axis": "T"}),
        },
    )
    ds = sgrid._attach_sgrid_metadata(ds, agrid_sgrid_metadata())
    return FieldSet.from_sgrid_conventions(ds, mesh=SphericalMesh())

@pytest.fixture
def convection_fieldset():
    """ Toy fieldset for test sims with convection_ou kernel """
    lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    U = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    V = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    W = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    # 3-D wind field with all winds 0 m/s
    fieldset = toy_fieldset({"U": U, "V": V, "W": W}, lon, lat, alt)
    fieldset.add_context("conv_sigma", 0.6)
    fieldset.add_context("conv_tau", 1200.0)
    fieldset.add_context("conv_z_lo", 48000.0)
    fieldset.add_context("conv_z_hi", 55000.0)
    fieldset.add_context("conv_edge", 2000.0)
    # Convection parameters
    return fieldset

@pytest.fixture
def balloon_fieldset():
    """ Toy fieldset for test sims with balloon kernel """
    lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    U = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    V = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    W = np.full((alt.size, lat.size, lon.size), 0.5, dtype=np.float32) 
    RHO = np.full((alt.size, lat.size, lon.size), 2.0, dtype=np.float32)
    # 3-D wind field with most winds 0 m/s
    fieldset = toy_fieldset({"U": U, "V": V, "W": W, 'RHO': RHO}, lon, lat, alt)
    fieldset.add_context("conv_sigma", 0.6)
    fieldset.add_context("conv_tau", 1200.0)
    fieldset.add_context("conv_z_lo", 48000.0)
    fieldset.add_context("conv_z_hi", 55000.0)
    fieldset.add_context("conv_edge", 2000.0)
    # Convection parameters
    fieldset.add_context("C_D_top", 0.8) # Drag coefficient for vertical motion
    fieldset.add_context("C_D_side", 1.0) # Drag coefficient for horizontal motion
    fieldset.add_context("C_m", 0.2) # Virtual mass coefficient
    fieldset.add_context("A_top", 19.6) # Upper area [m^2]
    fieldset.add_context("A_side", 22.9) # Silhouette area of profile [m^2]
    fieldset.add_context("g_Venus", 8.87) # Venus gravitational acceleration [m/s^2]
    fieldset.add_context("m_total", 62) # Total mass: helium + envelopes + payload [kg]
    fieldset.add_context("V_infl", 72.6) # Maximum volume [m^3]
    # Balloon parameters
    return fieldset

@pytest.fixture
def balloon_fieldset_variable_rho():
    """ Toy fieldset for test sims with balloon kernel with 
        variable rho for ceiling behaviour test """
    lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    rho_decay = np.array([68.4787*math.exp(-i/12900) for i in alt])
    U = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    V = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    W = np.zeros((alt.size, lat.size, lon.size), dtype=np.float32)
    RHO = np.tile(rho_decay[:, np.newaxis, np.newaxis], (1, lat.size, lon.size))
    # Still air, so only buoyancy can lift the balloon to its ceiling
    fieldset = toy_fieldset({"U": U, "V": V, "W": W, 'RHO': RHO}, lon, lat, alt)
    fieldset.add_context("conv_sigma", 0)
    fieldset.add_context("conv_tau", 1200.0)
    fieldset.add_context("conv_z_lo", 48000.0)
    fieldset.add_context("conv_z_hi", 55000.0)
    fieldset.add_context("conv_edge", 2000.0)
    # Convection parameters
    fieldset.add_context("C_D_top", 0.8) # Drag coefficient for vertical motion
    fieldset.add_context("C_D_side", 1.0) # Drag coefficient for horizontal motion
    fieldset.add_context("C_m", 0.2) # Virtual mass coefficient
    fieldset.add_context("A_top", 19.6) # Upper area [m^2]
    fieldset.add_context("A_side", 22.9) # Silhouette area of profile [m^2]
    fieldset.add_context("g_Venus", 8.87) # Venus gravitational acceleration [m/s^2]
    fieldset.add_context("m_total", 62) # Total mass: helium + envelopes + payload [kg]
    fieldset.add_context("V_infl", 72.6) # Maximum volume [m^3]
    # Balloon parameters
    return fieldset

def test_kernels_importable():
    assert callable(custom_kernels.smagdiff)
    assert callable(custom_kernels.convection_ou)
    assert callable(custom_kernels.balloon_vertical)

@pytest.mark.parametrize("z", [35000, 60000, 0])
def test_no_convection_outside_envelope(convection_fieldset, z):
    """ Check that particle displacement through convection is 0
        if particle is outside the vertical range where convection
        occurs """
    x0, y0, z0 = 6.1, 6.2, z # Initial particle location - cycle through z values
    pset = ParticleSet(convection_fieldset, pclass=VenusParticle, x=x0, y=y0, z=z0)

    pset.execute(convection_ou, runtime=4, dt=1) # Execute only convection_ou kernel
    assert pset.z == pytest.approx(z0) # Check that final particle depth is same as initial

def test_conv_distribution_stationary(convection_fieldset):
    """ Test that the mean of u_conv is zero and the variance
        is one over a reasonable sample of particles """
    n_particles = 5000 # Run an ensemble of particles to get a statistical sample
    tau = 1200.0 # Decorrelation time for convection - sim length should be much longer than this
    np.random.seed(1234) # Set random seed for reproducibility
    pset = ParticleSet(convection_fieldset, pclass=VenusParticle, 
                       x=np.full(n_particles, 6.1),
                       y=np.full(n_particles, 6.2),
                       z=np.full(n_particles, 52000.0)) # Initial location inside envelope

    pset.execute(convection_ou, runtime=10 * tau, dt=600) # Execute only convection_ou kernel

    u_conv = np.asarray(pset.u_conv) # Final u_conv values for all particles
    assert np.mean(u_conv) == pytest.approx(0.0, abs=0.05)
    assert np.var(u_conv) == pytest.approx(1.0, abs=0.1)

def test_balloon_tracks_watm(balloon_fieldset):
    """ Tests that the balloon kernel tracks W_atm and does
        not double-count the vertical wind. Convective noise is switched off:
        the kernel adds it to w_atm, so with it on w_bal tracks a random wind
        rather than the field value this test compares against. """
    tau = 1200.0
    balloon_fieldset.context["conv_sigma"] = 0.0
    w_atm = balloon_fieldset.W[0, 49000, 20, 90]
    balloon_fieldset.add_context('m_gas_ZP', (balloon_fieldset.m_total / 10.86))
    pset = ParticleSet(balloon_fieldset, pclass=BalloonParticle, x=90, y=20, z=49000)

    pset.execute(balloon_vertical, runtime=10*tau, dt=600)
    assert pset.w_bal == pytest.approx(w_atm, rel=0.01)

@pytest.mark.parametrize("m_gas_ZP", [62 / 10.86, 5.75])
def test_balloon_zero_relative_velocity(balloon_fieldset, m_gas_ZP):
    """ A balloon moving exactly with the air (w_rel = 0, so no drag) must give a
        finite update, and accelerates in the direction of its net buoyancy. """
    balloon_fieldset.context["conv_sigma"] = 0.0
    balloon_fieldset.add_context('m_gas_ZP', m_gas_ZP)
    w_atm = float(np.asarray(balloon_fieldset.W[0, 49000, 20, 90]).squeeze())
    pset = ParticleSet(balloon_fieldset, pclass=BalloonParticle, x=90, y=20, z=49000, w_bal=w_atm)

    pset.execute(balloon_vertical, runtime=600, dt=600)
    assert np.isfinite(pset.w_bal[0]) and np.isfinite(pset.z[0])
    if m_gas_ZP > 62 / 10.86:
        assert pset.w_bal[0] > w_atm
    else:
        assert pset.w_bal[0] == pytest.approx(w_atm, abs=1e-3)

def test_ceiling_behaviour(balloon_fieldset_variable_rho):
    """ Tests that the balloon arrests only when
        V = V_infl at the altitude where 
        m_total = rho_atm*V_infl """
    tau = 1200.0
    balloon_fieldset_variable_rho.add_context('m_gas_ZP', 5.75)
    pset = ParticleSet(balloon_fieldset_variable_rho, pclass=BalloonParticle, x=90, y=20, z=49000)

    pset.execute(balloon_vertical, runtime=100*tau, dt=600)
    final_depth = pset.z[0]
    rho_final = balloon_fieldset_variable_rho.RHO[0, final_depth, 20, 90]
    assert rho_final*balloon_fieldset_variable_rho.V_infl == pytest.approx(balloon_fieldset_variable_rho.m_total, rel=0.01)
    assert pset.w_bal == pytest.approx(0, abs=0.01)
    assert pset.v_bal == pytest.approx(balloon_fieldset_variable_rho.V_infl, rel=0.01)

@pytest.mark.parametrize("dt, ok", [(300, True), (120, False)])
def test_check_timestep(balloon_fieldset, dt, ok):
    """ The horizontal drag relaxation time at the reference case is ~15.3 s, so
        the shortest valid timestep is ~153 s: 300 s passes and 120 s raises. """
    balloon_fieldset.add_context("m_gas_ZP", 5.75)
    balloon_fieldset.add_context("rho_ceiling", 0.85)
    balloon_fieldset.add_context("v_rel_h_ref", 0.5)
    if ok:
        check_timestep(balloon_fieldset, dt)
    else:
        with pytest.raises(ValueError, match="too short"):
            check_timestep(balloon_fieldset, dt)

def uniform_wind_fieldset(u=0.0, v=0.0, w=0.0, lon=None, n_times=1):
    """ Toy fieldset with spatially uniform winds on a 1-degree grid """
    if lon is None:
        lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    shape = (alt.size, lat.size, lon.size)
    return toy_fieldset({"U": np.full(shape, u), "V": np.full(shape, v), "W": np.full(shape, w)},
                        lon, lat, alt, n_times=n_times)

@pytest.mark.parametrize("pclass", [MarsParticle, VenusParticle, BalloonParticle])
def test_boundary_stick_freezes_particle(pclass):
    """ A particle advected out through the northern edge of the domain is frozen
        at its last position and flagged stuck, instead of crashing the run.
        Covers every particle class used with boundary_stick. """
    fieldset = uniform_wind_fieldset(v=100.0) # 100 m/s northward: 0.54 deg per 600 s step
    y0 = 88.8 # Within one step of the last latitude node (89 deg)
    pset = ParticleSet(fieldset, pclass=pclass, x=10.0, y=y0, z=30000.0)

    pset.execute([AdvectionRK4, boundary_stick], runtime=1800, dt=600, verbose_progress=False)
    assert pset.stuck[0] == 1
    assert pset.y[0] == pytest.approx(y0)

def test_boundary_stick_leaves_inbounds_particle():
    """ A particle that stays inside the domain is advected normally and not flagged """
    fieldset = uniform_wind_fieldset(v=10.0)
    pset = ParticleSet(fieldset, pclass=MarsParticle, x=10.0, y=0.0, z=30000.0)

    pset.execute([AdvectionRK4, boundary_stick], runtime=1800, dt=600, verbose_progress=False)
    assert pset.stuck[0] == 0
    assert pset.y[0] > 0.0

@pytest.mark.parametrize("u, x0", [(50.0, 359.9), (-50.0, 0.1)])
def test_periodicBC_wraps_every_output_step(u, x0):
    """ A particle crossing the 0/360 seam (eastward and westward) is wrapped on
        the step it crosses, so its longitude never leaves [lon_min, lon_min + lon_span).
        The toy grid carries a 3-degree halo on each side, as the preprocessed files do. """
    lon = np.arange(-3, 363, dtype=np.float32)
    fieldset = uniform_wind_fieldset(u=u, lon=lon)
    fieldset.add_context("lon_min", 0.0)
    fieldset.add_context("lon_span", 360.0)
    pset = ParticleSet(fieldset, pclass=MarsParticle, x=x0, y=0.0, z=30000.0)

    step_deg = abs(u) * 600 / fieldset.UV.grid.deg2m # Longitude change per step at the equator
    xs = []
    for _ in range(4):
        pset.execute([AdvectionRK4, periodicBC, boundary_stick], runtime=600, dt=600, verbose_progress=False)
        xs.append(float(pset.x[0]))
    xs = np.array(xs)
    assert np.all((xs >= 0.0) & (xs < 360.0)), f"unwrapped longitude in {xs}"
    # Unwrapped path is a straight line: each step moves step_deg in the wind direction
    unwrapped = np.unwrap(np.deg2rad(np.concatenate([[x0], xs]))) * 180 / np.pi
    np.testing.assert_allclose(np.diff(unwrapped), np.sign(u) * step_deg, rtol=1e-3)

def test_surface_bounce_returns_particle_to_lowest_level():
    """ A particle driven down through the lowest model level is put back on it
        instead of raising ErrorThroughSurface. """
    fieldset = uniform_wind_fieldset(w=-1.0) # 1 m/s downward: 600 m per step
    pset = ParticleSet(fieldset, pclass=MarsParticle, x=10.0, y=0.0, z=100.0)

    pset.execute([AdvectionRK4_3D, surface_bounce, boundary_stick], runtime=600, dt=600, verbose_progress=False)
    assert pset.z[0] == pytest.approx(float(fieldset.W.grid.depth[0]))

def smagdiff_fieldset(u):
    """ Toy fieldset with winds u (3-D array), cell_areas and the Smagorinsky constant """
    lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    zeros = np.zeros_like(u)
    fieldset = toy_fieldset({"U": u, "V": zeros, "W": zeros}, lon, lat, alt, n_times=2)
    fieldset = add_cell_areas(fieldset, fieldset.UV.grid.mesh.radius)
    fieldset.add_context("Cs", 0.1)
    return fieldset

def test_smagdiff_no_diffusion_without_wind():
    """ With no wind there are no velocity gradients, so the Smagorinsky
        diffusivity and the random displacement are both zero. """
    fieldset = smagdiff_fieldset(np.zeros((35, 180, 360), dtype=np.float32))
    pset = ParticleSet(fieldset, pclass=MarsParticle, x=np.full(10, 100.0),
                       y=np.full(10, 10.0), z=np.full(10, 30000.0))

    pset.execute(smagdiff, runtime=1800, dt=600, verbose_progress=False)
    np.testing.assert_array_equal(pset.x, 100.0)
    np.testing.assert_array_equal(pset.y, 10.0)

def test_smagdiff_displacement_scales_with_root_strain():
    """ Kh is proportional to the strain rate, and the displacement to sqrt(Kh).
        Doubling a linear shear doubles every velocity gradient, so with the same
        random draws each displacement grows by exactly sqrt(2). """
    lat = np.arange(-90, 90, dtype=np.float32)
    shear = np.broadcast_to((0.5 * lat)[None, :, None], (35, 180, 360)).astype(np.float32)
    displacements = []
    for factor in (1.0, 2.0):
        fieldset = smagdiff_fieldset(factor * shear)
        np.random.seed(7)
        pset = ParticleSet(fieldset, pclass=MarsParticle, x=np.full(50, 100.0),
                           y=np.full(50, 10.0), z=np.full(50, 30000.0))
        pset.execute(smagdiff, runtime=600, dt=600, verbose_progress=False)
        displacements.append(np.asarray(pset.y, dtype=np.float64) - 10.0)
    np.testing.assert_allclose(displacements[1], np.sqrt(2.0) * displacements[0], rtol=1e-3, atol=1e-6)

def test_smagdiff_diffuses_in_sheared_flow():
    """ Meridional shear in U gives a non-zero diffusivity: particles spread out,
        with finite displacements centred on the release point. """
    lat = np.arange(-90, 90, dtype=np.float32)
    u = np.broadcast_to((0.5 * lat)[None, :, None], (35, 180, 360)).astype(np.float32) # dU/dlat = 0.5 m/s per degree
    fieldset = smagdiff_fieldset(u)
    np.random.seed(42)
    n = 2000
    pset = ParticleSet(fieldset, pclass=MarsParticle, x=np.full(n, 100.0),
                       y=np.full(n, 10.0), z=np.full(n, 30000.0))

    pset.execute(smagdiff, runtime=1800, dt=600, verbose_progress=False)
    dy = np.asarray(pset.y) - 10.0
    assert np.all(np.isfinite(dy))
    assert np.std(dy) > 0.0
    assert abs(np.mean(dy)) < 3 * np.std(dy) / np.sqrt(n)
