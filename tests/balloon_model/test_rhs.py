import numpy as np
import pytest
from balloon_model import state
from balloon_model.params import BalloonParams, ModelConfig
from balloon_model.environment import UniformAtmosphere
from balloon_model.rhs import displaced_volume, vertical_acceleration, total_mass, rhs

P = BalloonParams()
CFG = ModelConfig()


def test_volume_free_and_capped():
    """ Free branch: V = m (M_a / M_gas) / rho_a, so rho_a V = 10.855 m_gas
        whatever the density; capped at V_infl when that exceeds it. """
    rho = np.array([2.0, 0.5])
    V = displaced_volume(5.75, rho, P, CFG)
    assert V[0] == pytest.approx(5.75 * (43.45 / 4.0026) / 2.0, rel=1e-12)
    assert V[1] == P.V_infl


def test_total_mass_uses_state_gas():
    assert total_mass(5.0, 1.0, P) == pytest.approx(P.m_struct + 6.0)


def test_no_acceleration_when_neutral_and_moving_with_air():
    """ rho_a V = m_tot and w = w_a: no net force. """
    m_tot = 62.0
    rho = 0.9
    V = m_tot / rho
    assert vertical_acceleration(0.3, 0.3, rho, V, m_tot, P, CFG) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("w_r", [0.5, -0.5])
def test_drag_opposes_relative_motion(w_r):
    """ Neutral buoyancy: the acceleration has the sign opposite to w_r and
        magnitude 0.5 rho C_D A w_r^2 / (m_tot + C_m rho V). """
    m_tot, rho = 62.0, 0.9
    V = m_tot / rho
    a = vertical_acceleration(0.1 + w_r, 0.1, rho, V, m_tot, P, CFG)
    expected = -np.sign(w_r) * 0.5 * rho * P.C_D_top * P.A_top * w_r**2 / (m_tot + P.C_m * rho * V)
    assert a == pytest.approx(expected, rel=1e-12)


def test_rhs_tendencies():
    """ dz/dt = w; dw/dt from the momentum equation; every other variable has
        zero tendency in the minimal model. """
    env = UniformAtmosphere(rho=0.9, T=300.0, P=5e4, w=0.2)
    y = state.pack(z=np.array([50000.0, 51000.0]), w=np.array([0.0, 1.0]), m_ZP=5.75, T_ZP=300.0)
    d = rhs(0.0, y, P, CFG, env)
    np.testing.assert_array_equal(d[state.IDX["z"]], [0.0, 1.0])
    V = displaced_volume(5.75, 0.9, P, CFG)
    expected = vertical_acceleration(np.array([0.0, 1.0]), 0.2, 0.9, V, P.m_struct + 5.75, P, CFG)
    np.testing.assert_array_equal(d[state.IDX["w"]], expected)
    others = [i for n, i in state.IDX.items() if n not in ("z", "w")]
    assert np.all(d[others] == 0.0)


def test_rhs_counts_superpressure_gas_in_mass():
    """ Gas in the SP balloon adds to the mass but not to the ZP volume. """
    env = UniformAtmosphere(rho=0.9, T=300.0, P=5e4)
    y = state.pack(z=50000.0, m_ZP=5.75, m_SP=1.0)
    d = rhs(0.0, y, P, CFG, env)
    V = displaced_volume(5.75, 0.9, P, CFG)
    expected = vertical_acceleration(0.0, 0.0, 0.9, V, P.m_struct + 6.75, P, CFG)
    np.testing.assert_allclose(d[state.IDX["w"]], expected, rtol=1e-14)
