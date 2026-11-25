""" Integration tests of the isothermal model: neutral buoyancy, analytic
    terminal velocity, equilibrium altitude, zero relative velocity, kernel
    scheme vs truth, and input checks. """
import numpy as np
import pytest
from balloon_model import state
from balloon_model.params import BalloonParams, ModelConfig
from balloon_model.environment import IsothermalAtmosphere, UniformAtmosphere, ProfileAtmosphere
from balloon_model.integrate import integrate_truth, integrate_kernel_scheme, kernel_step, linearise, DZ_PROBE
from balloon_model import rhs as rhs_module

P = BalloonParams()
CFG = ModelConfig()
M_TOT = P.m_struct + P.m_gas_ZP     # 62 kg
DAY = 86400.0
Z, W = state.IDX["z"], state.IDX["w"]


def venus_isothermal(w=0.0):
    return IsothermalAtmosphere(T=300.0, P_ref=5.0e4, z_ref=55000.0, M_atm=P.M_atm, g=P.g, w=w)


def equilibrium_altitude(env):
    """ Altitude where rho_a V_infl = m_tot. The default craft has gas lift
        10.855 * 5.75 = 62.42 kg > 62 kg, so it floats on the capped branch. """
    return env.altitude_of_density(M_TOT / P.V_infl)


def test_neutral_buoyancy_holds_for_a_day():
    """ Released at rest at its equilibrium altitude, the balloon stays there.
        The only forcing is round-off in rho_a V - m_tot (~1e-14 kg), so both
        integrators must hold z to well under a millimetre. """
    env = venus_isothermal()
    z_eq = equilibrium_altitude(env)
    y0 = state.pack(z=z_eq, m_ZP=P.m_gas_ZP)
    _, Y = integrate_truth(env, y0, P, CFG, (0.0, DAY), t_eval=[DAY])
    assert abs(Y[-1, Z, 0] - z_eq) < 1e-3
    _, Yk = integrate_kernel_scheme(env, y0, P, CFG, 0.0, 600.0, 144, 60)
    assert np.max(np.abs(Yk[:, Z, 0] - z_eq)) < 1e-3


def terminal_solution(t, rho):
    """ Rise from rest in still uniform air with constant net force F:
        m_eff dw/dt = F - k w^2, k = 0.5 rho C_D A, so
        w = w_t tanh(t / tau), z - z0 = w_t tau ln cosh(t / tau),
        w_t = sqrt(F / k), tau = m_eff / sqrt(F k). """
    V = P.m_gas_ZP * (P.M_atm / P.M_gas) / rho
    F = (rho * V - M_TOT) * P.g
    k = 0.5 * rho * P.C_D_top * P.A_top
    m_eff = M_TOT + P.C_m * rho * V
    w_t, tau = np.sqrt(F / k), m_eff / np.sqrt(F * k)
    return w_t, w_t * np.tanh(t / tau), w_t * tau * np.log(np.cosh(t / tau))


def test_terminal_velocity_matches_analytic():
    """ Uniform density keeps the free-branch volume, and so the net force,
        constant. The truth integrator must follow the tanh solution to its
        tolerance (1e-10 relative, so 1e-7 is a safe bound), and reach w_t.
        Samples every 2 s resolve the approach, whose time scale tau is 14.6 s. """
    rho = 0.9    # Free branch: V = 69.4 m^3 < V_infl
    env = UniformAtmosphere(rho=rho, T=300.0, P=5e4)
    t = np.arange(0.0, 3601.0, 2.0)
    y0 = state.pack(z=50000.0, m_ZP=P.m_gas_ZP)
    _, Y = integrate_truth(env, y0, P, CFG, (0.0, t[-1]), t_eval=t)
    w_t, w, dz = terminal_solution(t, rho)
    np.testing.assert_allclose(Y[:, W, 0], w, rtol=1e-7, atol=1e-10)
    np.testing.assert_allclose(Y[:, Z, 0] - 50000.0, dz, rtol=1e-7, atol=1e-6)
    assert Y[-1, W, 0] == pytest.approx(w_t, rel=1e-9)


def test_settles_at_equilibrium_altitude():
    """ Released 2 km below its equilibrium altitude, the balloon rises,
        overshoots, oscillates and settles where rho_a V_infl = m_tot.
        Quadratic drag damps the oscillation ever more weakly as it shrinks
        (amplitude falls roughly as 1/t), so after one day a residual of a few
        millimetres remains; 0.05 m bounds it with margin. """
    env = venus_isothermal()
    z_eq = equilibrium_altitude(env)
    y0 = state.pack(z=z_eq - 2000.0, m_ZP=P.m_gas_ZP)
    t, Y = integrate_truth(env, y0, P, CFG, (0.0, DAY), t_eval=np.arange(0.0, DAY + 1, 60.0))
    z = Y[:, Z, 0]
    assert z.max() - z_eq > 1.0                        # overshoot
    assert abs(z[-1] - z_eq) < 0.05                    # settled
    assert np.max(np.abs(z[t > 0.5 * DAY] - z_eq)) < 0.1


@pytest.mark.parametrize("excess", [0.0, 0.05])
def test_zero_relative_velocity_is_finite(excess):
    """ A balloon moving exactly with the air has no drag. The kernel scheme
        floors |w_rel| in its relaxation time, so the update stays finite;
        neutral, the balloon keeps moving with the air, and with extra gas
        it accelerates upward relative to it. """
    env = UniformAtmosphere(rho=0.9, T=300.0, P=5e4, w=0.3)
    # Free branch: rho V = 10.855 m_gas, so neutral gas mass solves 10.855 m = m_struct + m
    m_ratio = P.M_atm / P.M_gas
    m_gas = P.m_struct / (m_ratio - 1.0) * (1.0 + excess)
    y0 = state.pack(z=50000.0, w=0.3, m_ZP=m_gas)
    y, _ = kernel_step(y0, linearise(env, 0.0, state.get(y0, "z")), P, CFG, 600.0, 60)
    w = y[W, 0]
    assert np.all(np.isfinite(y))
    if excess == 0.0:
        assert w == pytest.approx(0.3, abs=1e-9)
    else:
        assert w > 0.3


def linear_density_env(w=0.3):
    """ Density exactly linear in z (one profile segment), so the kernel
        scheme's start-of-step density gradient is exact and its error comes
        from substepping alone. Slope -rho/H with H = 6472 m, as near 55 km. """
    z = np.array([40000.0, 70000.0])
    slope = -0.854 / 6472.0
    one = np.ones(2)
    return ProfileAtmosphere(z, rho=0.854 + slope * (z - 55000.0), T=300.0 * one, P=5e4 * one,
                             u=0.0 * one, v=0.0 * one, w=w * one)


def test_kernel_scheme_converges_to_truth_at_first_order():
    """ Rise from 2 km below equilibrium through the volume cap, with a 0.3 m/s
        updraft, over 2 h of 600 s outer steps. The scheme freezes the drag
        coefficient and buoyancy over each substep and moves z with the
        end-of-substep velocity, so its error is first order: doubling the
        substep count halves the z and w errors (ratio 2 within 15% once in
        the asymptotic range, n_sub >= 80). At the production n_sub = 60 the z error is
        below 10 m on a 2 km rise and the w error below 1e-3 m/s. """
    env = linear_density_env()
    slope = -0.854 / 6472.0
    z_eq = 55000.0 + (M_TOT / P.V_infl - 0.854) / slope
    y0 = state.pack(z=z_eq - 2000.0, m_ZP=P.m_gas_ZP)
    T = 7200.0
    t_eval = np.arange(0.0, T + 1, 600.0)
    _, Y = integrate_truth(env, y0, P, CFG, (0.0, T), t_eval=t_eval)
    errors = {}
    for n_sub in (60, 80, 160, 320):
        _, Yk = integrate_kernel_scheme(env, y0, P, CFG, 0.0, 600.0, int(T / 600), n_sub)
        errors[n_sub] = (np.max(np.abs(Yk[:, Z, 0] - Y[:, Z, 0])), np.max(np.abs(Yk[:, W, 0] - Y[:, W, 0])))
    assert errors[60][0] < 10.0 and errors[60][1] < 1e-3
    for n in (80, 160):
        assert errors[n][0] / errors[2 * n][0] == pytest.approx(2.0, rel=0.15)
        assert errors[n][1] / errors[2 * n][1] == pytest.approx(2.0, rel=0.15)


def test_kernel_scheme_density_linearisation_error_is_bounded():
    """ In an exponential (isothermal) atmosphere the start-of-step density
        gradient misses the curvature of rho(z); over a 600 s step at ~0.7 m/s
        the density error reaches (dz/H)^2/2 ~ 0.2%, a sizeable share of the
        0.4 kg net lift. This error does not shrink with n_sub. Over a day
        from 2 km below equilibrium it stays below 20 m (14 m found), and both
        integrators settle at the same altitude to within 0.05 m. """
    env = venus_isothermal()
    z_eq = equilibrium_altitude(env)
    y0 = state.pack(z=z_eq - 2000.0, m_ZP=P.m_gas_ZP)
    t_eval = np.arange(0.0, DAY + 1, 600.0)
    _, Y = integrate_truth(env, y0, P, CFG, (0.0, DAY), t_eval=t_eval)
    _, Yk = integrate_kernel_scheme(env, y0, P, CFG, 0.0, 600.0, 144, 60)
    assert np.max(np.abs(Yk[:, Z, 0] - Y[:, Z, 0])) < 20.0
    assert abs(Yk[-1, Z, 0] - z_eq) < 0.05


def test_kernel_step_vectorised_over_particles():
    """ Particles with different states advance independently: a batch step
        equals single-particle steps. """
    env = venus_isothermal(w=0.1)
    z = np.array([50000.0, 55000.0, 58000.0])
    w = np.array([-0.5, 0.0, 0.8])
    y0 = state.pack(z=z, w=w, m_ZP=P.m_gas_ZP)
    y, _ = kernel_step(y0, linearise(env, 0.0, z), P, CFG, 600.0, 60)
    for i in range(3):
        yi0 = y0[:, i:i + 1]
        yi, _ = kernel_step(yi0, linearise(env, 0.0, z[i:i + 1]), P, CFG, 600.0, 60)
        np.testing.assert_array_equal(yi[:, 0], y[:, i])


def test_linearise_near_profile_edges():
    """ Within DZ_PROBE of either end of a profile the gradient is a one-sided
        difference; for a density linear in z it is still exact. """
    env = linear_density_env()
    slope = -0.854 / 6472.0
    z = np.array([40000.0, 40000.0 + 0.5 * DZ_PROBE, 55000.0, 70000.0 - 0.5 * DZ_PROBE, 70000.0])
    np.testing.assert_allclose(linearise(env, 0.0, z).drho_dz, slope, rtol=1e-9)


def test_truth_integrates_only_evolved_variables():
    """ The solver passes the right-hand side only finite states, and the
        variables the configuration does not evolve keep their initial values. """
    env = venus_isothermal()
    calls = {"nonfinite": 0}

    class Checked:
        def sample(self, t, z, lat=None, lon=None):
            calls["nonfinite"] += int(not np.all(np.isfinite(z)))
            return env.sample(t, z)

    y0 = state.pack(z=equilibrium_altitude(env) - 500.0, m_ZP=P.m_gas_ZP, T_ZP=300.0, T_e1=290.0)
    _, Y = integrate_truth(Checked(), y0, P, CFG, (0.0, 3600.0), t_eval=[0.0, 3600.0])
    assert calls["nonfinite"] == 0
    frozen = [i for name, i in state.IDX.items() if name not in rhs_module.evolved_variables(CFG)]
    np.testing.assert_array_equal(Y[-1, frozen], y0[frozen])
    assert Y[-1, Z, 0] != y0[Z, 0]


def test_kernel_scheme_resamples_time_varying_environment():
    """ Each outer step samples the environment at its own start time. """
    class Gust:
        def sample(self, t, z, lat=None, lon=None):
            return UniformAtmosphere(rho=0.9, T=300.0, P=5e4, w=0.0 if t < 600.0 else 0.5).sample(t, z)

    env = Gust()
    y0 = state.pack(z=50000.0, m_ZP=P.m_gas_ZP)
    _, Yk = integrate_kernel_scheme(env, y0, P, CFG, 0.0, 600.0, 2, 60)
    y1, _ = kernel_step(y0, linearise(env, 0.0, state.get(y0, "z")), P, CFG, 600.0, 60)
    y2, _ = kernel_step(y1, linearise(env, 600.0, state.get(y1, "z")), P, CFG, 600.0, 60)
    np.testing.assert_array_equal(Yk[2], y2)


@pytest.mark.parametrize("n_sub", [0, -1, 2.5])
def test_kernel_step_rejects_bad_substep_count(n_sub):
    env = UniformAtmosphere(rho=0.9, T=300.0, P=5e4)
    y0 = state.pack(z=50000.0, m_ZP=P.m_gas_ZP)
    with pytest.raises(ValueError, match="n_sub"):
        kernel_step(y0, linearise(env, 0.0, state.get(y0, "z")), P, CFG, 600.0, n_sub)


def test_kernel_step_rejects_zero_drag():
    env = UniformAtmosphere(rho=0.9, T=300.0, P=5e4)
    y0 = state.pack(z=50000.0, m_ZP=P.m_gas_ZP)
    with pytest.raises(ValueError, match="C_D_top"):
        kernel_step(y0, linearise(env, 0.0, state.get(y0, "z")), BalloonParams(C_D_top=0.0), CFG, 600.0, 60)
