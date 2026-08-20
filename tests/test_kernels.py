import custom_kernels, parcels
import pytest
import numpy as np
import xarray as xr
from custom_kernels import convection_ou, balloon_vertical, VenusParticle, BalloonParticle
from parcels import FieldSet, ParticleSet, ScipyParticle, Variable

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
    fieldset = FieldSet.from_data({"U": U, "V": V, "W": W},
                                {"lon": lon, "lat": lat, "depth": alt})
    fieldset.add_constant("conv_sigma", 0.6)
    fieldset.add_constant("conv_tau", 1200.0)
    fieldset.add_constant("conv_z_lo", 48000.0)
    fieldset.add_constant("conv_z_hi", 55000.0)
    fieldset.add_constant("conv_edge", 2000.0)
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
    fieldset = FieldSet.from_data({"U": U, "V": V, "W": W, 'RHO': RHO},
                                {"lon": lon, "lat": lat, "depth": alt})
    fieldset.add_constant("conv_sigma", 0.6)
    fieldset.add_constant("conv_tau", 1200.0)
    fieldset.add_constant("conv_z_lo", 48000.0)
    fieldset.add_constant("conv_z_hi", 55000.0)
    fieldset.add_constant("conv_edge", 2000.0)
    # Convection parameters
    fieldset.add_constant("C_D_top", 0.8) # Drag coefficient for vertical motion
    fieldset.add_constant("C_D_side", 1.0) # Drag coefficient for horizontal motion
    fieldset.add_constant("C_m", 0.2) # Virtual mass coefficient
    fieldset.add_constant("A_top", 19.6) # Upper area [m^2]
    fieldset.add_constant("A_side", 22.9) # Silhouette area of profile [m^2]
    fieldset.add_constant("g_Venus", 8.87) # Venus gravitational acceleration [m/s^2]
    fieldset.add_constant("m_total", 62) # Total mass: helium + envelopes + payload [kg]
    #fieldset.add_constant("m_gas_ZP", 5.75) # Mass of helium in the ZP balloon [kg]
    fieldset.add_constant("V_infl", 72.6) # Maximum volume [m^3]
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
    pset = ParticleSet(convection_fieldset, pclass=VenusParticle, lon=x0, lat=y0, depth=z0)

    pset.execute(convection_ou, runtime=4, dt=1) # Execute only convection_ou kernel
    assert pset.depth == pytest.approx(z0) # Check that final particle depth is same as initial

def test_conv_distribution_stationary(convection_fieldset):
    """ Test that the mean of u_conv is zero and the variance
        is one over a reasonable sample of particles """
    n_particles = 5000 # Run an ensemble of particles to get a statistical sample
    tau = 1200.0 # Decorrelation time for convection - sim length should be much longer than this
    parcels.ParcelsRandom.seed(1234) # Set random seed for reproducibility
    pset = ParticleSet(convection_fieldset, pclass=VenusParticle, 
                       lon=np.full(n_particles, 6.1), 
                       lat=np.full(n_particles, 6.2), 
                       depth=np.full(n_particles, 52000.0)) # Initial location inside envelope

    pset.execute(convection_ou, runtime=10 * tau, dt=600) # Execute only convection_ou kernel

    u_conv = np.asarray(pset.u_conv) # Final u_conv values for all particles
    assert np.mean(u_conv) == pytest.approx(0.0, abs=0.05)
    assert np.var(u_conv) == pytest.approx(1.0, abs=0.1)

def test_balloon_tracks_watm(balloon_fieldset):
    """ Tests that the balloon kernel tracks W_atm and does 
        not double-count the vertical wind """
    n_particles = 5000
    tau = 1200.0 
    #w_atm = balloon_fieldset.W
    w_atm = 0.5
    balloon_fieldset.add_constant('m_gas_ZP', (balloon_fieldset.m_total / 10.86))
    pset = ParticleSet(balloon_fieldset, pclass=BalloonParticle, lon=90, lat=20, depth=49000)

    pset.execute(balloon_vertical, runtime=10*tau, dt=600)
    w_bal = pset.w_bal
    assert w_bal == pytest.approx(w_atm, rel=(w_atm*0.01))

@pytest.mark.parametrize("v_rel", [0.5, 1, 2, 3])
def test_tau_horizontal(balloon_fieldset, v_rel):
    """ Tests that tau_horizontal is well below
        the outer dt and is around 10s """
    tau = 1200.0
    #rho_atm = balloon_fieldset.RHO # how to multiply with uniform fields?
    rho_atm = 2.0
    balloon_fieldset.add_constant('m_gas_ZP', 5.75)
    pset = ParticleSet(balloon_fieldset, pclass=BalloonParticle, lon=90, lat=20, depth=49000)

    pset.execute(balloon_vertical, runtime=10*tau, dt=600)
    m_virtual = balloon_fieldset.C_m*rho_atm*pset.v_bal
    tau_horizontal = (balloon_fieldset.m_total + m_virtual)/(0.5*rho_atm*balloon_fieldset.C_D_side*balloon_fieldset.A_side*v_rel)
    assert tau_horizontal <= 10

def test_ceiling_behaviour(balloon_fieldset):
    """ Tests that the balloon arrests only when
        V = V_infl at the altitude where 
        m_total = rho_atm*V_infl """
    tau = 1200.0
    balloon_fieldset.add_constant('m_gas_ZP', 5.75)