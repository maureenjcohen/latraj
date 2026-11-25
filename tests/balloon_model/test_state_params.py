import numpy as np
import pytest
from balloon_model import state
from balloon_model.params import BalloonParams, ModelConfig, M_HE


def test_state_layout_holds_full_model():
    """ The layout carries every prognostic variable of the full model: two gas
        nodes, four envelope nodes, horizontal velocity and the convective gust
        variable. """
    for name in ("z", "w", "u", "v", "m_ZP", "m_SP", "T_ZP", "T_SP", "T_e1", "T_e2", "T_e3", "T_e4",
                 "u_conv"):
        assert name in state.IDX
    assert state.N_STATE == len(state.STATE_VARS) == len(set(state.STATE_VARS))


def test_pack_and_get():
    y = state.pack(z=np.array([1.0, 2.0, 3.0]), m_ZP=5.75)
    assert y.shape == (state.N_STATE, 3) and y.dtype == np.float64
    np.testing.assert_array_equal(state.get(y, "z"), [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(state.get(y, "m_ZP"), 5.75)
    np.testing.assert_array_equal(state.get(y, "T_e4"), 0.0)
    with pytest.raises(KeyError):
        state.pack(altitude=1.0)


def test_default_config_is_minimal_model():
    cfg = ModelConfig()
    assert (cfg.thermal_model, cfg.drag_model, cfg.area_model) == ("isothermal", "constant", "fixed")


@pytest.mark.parametrize("switch, value", [("thermal_model", "one_node"), ("drag_model", "Re"),
                                           ("area_model", "shape"), ("horizontal", "drag_lag"),
                                           ("balloon_type", "superpressure")])
def test_planned_switch_values_raise(switch, value):
    with pytest.raises(NotImplementedError):
        ModelConfig(**{switch: value})


def test_unknown_switch_value_raises():
    with pytest.raises(ValueError):
        ModelConfig(drag_model="cubic")


def test_balloon_vertical_params():
    """ The equivalent of the balloon_vertical constants keeps the total mass
        and the 10.86 lift ratio of that kernel. """
    p = BalloonParams.balloon_vertical(m_total=62.0, m_gas_ZP=5.75)
    assert p.m_struct + p.m_gas_ZP + p.m_gas_SP == pytest.approx(62.0, abs=1e-12)
    assert p.M_atm / p.M_gas == pytest.approx(10.86, rel=1e-12)


def test_balloon_vertical_overrides():
    """ M_atm may be overridden; m_payload is derived from m_total, so passing
        it is an error rather than being silently replaced. """
    p = BalloonParams.balloon_vertical(m_total=62.0, m_gas_ZP=5.75, M_atm=43.45e-3)
    assert p.M_atm == 43.45e-3
    with pytest.raises(TypeError, match="m_payload"):
        BalloonParams.balloon_vertical(m_total=62.0, m_gas_ZP=5.75, m_payload=10.0)


def test_default_total_mass_and_gas_constants():
    p = BalloonParams()
    assert p.m_struct + p.m_gas_ZP + p.m_gas_SP == pytest.approx(62.0, abs=1e-12)
    assert p.M_gas == M_HE
    assert p.R_gas == pytest.approx(2077.3, rel=1e-4)
    assert p.cv_gas == pytest.approx(1.5 * 8.314462618 / 4.0026e-3, rel=1e-12)
    assert p.cp_gas / p.cv_gas == pytest.approx(5.0 / 3.0)
    assert p.R_atm == pytest.approx(8.314462618 / 43.45e-3, rel=1e-12)


def test_struct_mass_includes_both_envelopes():
    assert BalloonParams(m_env_SP=2.0).m_struct == pytest.approx(BalloonParams().m_struct + 2.0)


def test_params_accept_per_particle_arrays():
    """ Any parameter can be an array over particles. """
    p = BalloonParams(m_payload=np.array([30.0, 40.0]))
    np.testing.assert_allclose(p.m_struct, [51.3, 61.3])
