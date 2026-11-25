""" The minimal reference model reproduces the existing balloon_vertical kernel. """
import math
import numpy as np
import pytest
from parcels import ParticleSet
from custom_kernels import balloon_vertical, BalloonParticle
from tests.test_kernels import toy_fieldset
from balloon_model import state
from balloon_model.params import BalloonParams, ModelConfig
from balloon_model.environment import ProfileAtmosphere
from balloon_model.integrate import kernel_step, linearise

W_ATM = 0.3
CONTEXT = dict(conv_sigma=0.0, conv_tau=1200.0, conv_z_lo=48000.0, conv_z_hi=55000.0, conv_edge=2000.0,
               C_D_top=0.8, C_m=0.2, A_top=19.6, g_Venus=8.87, m_total=62.0, V_infl=72.6, m_gas_ZP=5.75)


def as_float32(a):
    """ The values the kernel sees: field data and particle variables are float32. """
    return np.asarray(a, dtype=np.float32).astype(np.float64)


@pytest.fixture(scope="module")
def setups():
    """ Exponential density (scale height 12.9 km) and a uniform 0.3 m/s
        updraft, as a Parcels FieldSet and as the equivalent column profile. """
    lon = np.arange(0, 360, dtype=np.float32)
    lat = np.arange(-90, 90, dtype=np.float32)
    alt = np.arange(0, 70000, 2000, dtype=np.float32)
    rho = np.array([68.4787 * math.exp(-z / 12900) for z in alt], dtype=np.float32)
    shape = (alt.size, lat.size, lon.size)
    fieldset = toy_fieldset({"U": np.zeros(shape), "V": np.zeros(shape), "W": np.full(shape, W_ATM),
                             "RHO": np.broadcast_to(rho[:, None, None], shape)}, lon, lat, alt)
    for name, value in CONTEXT.items():
        fieldset.add_context(name, value)
    one = np.ones(alt.size)
    env = ProfileAtmosphere(as_float32(alt), rho=as_float32(rho), T=300.0 * one, P=1e5 * one,
                            u=0.0 * one, v=0.0 * one, w=as_float32(W_ATM) * one)
    return fieldset, env


@pytest.mark.parametrize("z0, w0", [
    (49000.0, -0.4),    # free branch, rising after an initial descent
    (56300.0, 0.2),     # capped branch near equilibrium
    (57500.0, 1.0),     # capped branch above equilibrium, sinking
])
def test_matches_balloon_vertical(setups, z0, w0):
    """ One 600 s outer step with convective noise off. Both runs use the same
        scheme, so they differ only by float32 storage in the kernel: z to one
        float32 spacing at 50-60 km (3.9e-3 m; 1.5e-3 m found), w_bal to 1e-5
        m/s (4e-6 found), v_bal to 1e-4 m^3 (5e-6 found). """
    fieldset, env = setups
    pset = ParticleSet(fieldset, pclass=BalloonParticle, x=90, y=20, z=z0, w_bal=w0)
    pset.execute(balloon_vertical, runtime=600, dt=600, verbose_progress=False)

    p = BalloonParams.balloon_vertical(m_total=CONTEXT["m_total"], m_gas_ZP=CONTEXT["m_gas_ZP"],
                                       C_D_top=CONTEXT["C_D_top"], C_m=CONTEXT["C_m"],
                                       A_top=CONTEXT["A_top"], g=CONTEXT["g_Venus"],
                                       V_infl=CONTEXT["V_infl"])
    y0 = state.pack(z=as_float32(z0), w=as_float32(w0), m_ZP=p.m_gas_ZP)
    y, V = kernel_step(y0, linearise(env, 0.0, state.get(y0, "z")), p, ModelConfig(), 600.0, 60)

    assert abs(float(pset.z[0]) - y[state.IDX["z"], 0]) < 3.9e-3
    assert abs(float(pset.w_bal[0]) - y[state.IDX["w"], 0]) < 1e-5
    assert abs(float(pset.v_bal[0]) - V[0]) < 1e-4
    assert abs(y[state.IDX["z"], 0] - z0) > 100.0     # the step moves the balloon
