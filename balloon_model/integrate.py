""" Integrators for the balloon model.

    integrate_truth: scipy solve_ivp with an implicit stiff method and tight
    tolerances; the accuracy reference.

    kernel_step: the fixed-substep scheme the Parcels kernel calls. It sees
    the environment only through values sampled at the start of the outer
    step (LocalEnvironment), as a kernel does, and loops over substeps in
    time, not over particles. integrate_kernel_scheme drives kernel_step
    through a column environment for testing.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from . import state
from .rhs import displaced_volume, total_mass, rhs, evolved_variables

DZ_PROBE = 200.0        # Half-separation of the density probes [m]
W_REL_FLOOR = 1e-6      # Lower bound on |w_rel| in the drag time scale [m/s]


@dataclass
class LocalEnvironment:
    """ Environment at the start of an outer step: density, its vertical
        gradient and the vertical wind, one value per particle. """
    rho0: np.ndarray
    drho_dz: np.ndarray
    w_a: np.ndarray


def linearise(env, t, z):
    """ Sample env at z and at z +/- DZ_PROBE for a centred density gradient.
        An environment with a finite altitude range (attribute z_range) clips
        the probes to it, so within DZ_PROBE of an edge the gradient is a
        one-sided difference over the span that remains. """
    z = np.asarray(z, dtype=np.float64)
    z_lo, z_hi = getattr(env, "z_range", (-np.inf, np.inf))
    a0 = env.sample(t, z)
    z_up = np.minimum(z + DZ_PROBE, z_hi)
    z_dn = np.maximum(z - DZ_PROBE, z_lo)
    rho_up = env.sample(t, z_up)["rho"]
    rho_dn = env.sample(t, z_dn)["rho"]
    return LocalEnvironment(rho0=a0["rho"], drho_dz=(rho_up - rho_dn) / (z_up - z_dn), w_a=a0["w"])


def kernel_step(y, local, p, cfg, dt, n_sub):
    """ Advance the state y by one outer step dt with n_sub substeps.

        In each substep the drag is linearised about the relative velocity at
        the substep start, which makes the momentum equation

            (m_tot + C_m rho V) dw/dt = (rho V - m_tot) g - 0.5 rho C_D A_top |w_r| w_r

        linear in w_r with relaxation time tau and equilibrium w_eq = a_buoy tau;
        that linear equation is advanced exactly. |w_r| is floored so tau stays
        finite when the balloon moves with the air. Density follows the
        balloon by the start-of-step gradient; the wind is held fixed.

        Returns the new state and the displaced volume at the start of the
        last substep. Needs n_sub >= 1 and non-zero drag, since tau is
        inversely proportional to C_D_top A_top. """
    if cfg.drag_model != "constant" or cfg.area_model != "fixed":
        raise NotImplementedError((cfg.drag_model, cfg.area_model))
    if int(n_sub) != n_sub or n_sub < 1:
        raise ValueError(f"n_sub must be a positive integer, got {n_sub!r}")
    if np.any(np.asarray(p.C_D_top * p.A_top) <= 0):
        raise ValueError("kernel_step needs C_D_top * A_top > 0")
    y = np.array(y, dtype=np.float64)
    dt_sub = dt / n_sub
    w = state.get(y, "w").copy()
    m_ZP = state.get(y, "m_ZP")
    m_tot = total_mass(m_ZP, state.get(y, "m_SP"), p)
    rho = local.rho0
    displacement = np.zeros_like(w)
    for _ in range(n_sub):
        V = displaced_volume(m_ZP, rho, p, cfg)
        m_eff = m_tot + p.C_m * rho * V
        w_rel_old = w - local.w_a
        tau = m_eff / (0.5 * rho * p.C_D_top * p.A_top * np.maximum(np.abs(w_rel_old), W_REL_FLOOR))
        a_buoy = (rho * V - m_tot) * p.g / m_eff
        w_eq = a_buoy * tau
        w_rel = w_eq + (w_rel_old - w_eq) * np.exp(-np.abs(dt_sub) / tau)
        w = w_rel + local.w_a
        displacement += w * dt_sub
        rho = local.rho0 + local.drho_dz * displacement
    y[state.IDX["z"]] += displacement
    y[state.IDX["w"]] = w
    return y, V


def integrate_kernel_scheme(env, y0, p, cfg, t0, dt, n_steps, n_sub):
    """ Run kernel_step for n_steps outer steps of length dt, re-sampling the
        environment at the start of each. Returns times (n_steps + 1,) and
        states (n_steps + 1, N_STATE, n). """
    y = np.array(y0, dtype=np.float64)
    times = t0 + dt * np.arange(n_steps + 1)
    out = np.empty((n_steps + 1,) + y.shape)
    out[0] = y
    for k in range(n_steps):
        local = linearise(env, times[k], state.get(y, "z"))
        y, _ = kernel_step(y, local, p, cfg, dt, n_sub)
        out[k + 1] = y
    return times, out


def integrate_truth(env, y0, p, cfg, t_span, t_eval=None, method="Radau", rtol=1e-10, atol=1e-10):
    """ Integrate the full right-hand side with solve_ivp. Returns times and
        states (len(t), N_STATE, n). Raises if the solver fails.

        Only the variables the configuration evolves are passed to the solver;
        the others are constant, and the finite-difference Jacobian of an
        implicit method perturbs a zero-valued constant row by an amount that
        overflows to inf. """
    y0 = np.array(y0, dtype=np.float64)
    shape = y0.shape
    rows = [state.IDX[name] for name in evolved_variables(cfg)]
    active_shape = (len(rows), shape[1])

    def full(yflat):
        y = y0.copy()
        y[rows] = yflat.reshape(active_shape)
        return y

    def f(t, yflat):
        return rhs(t, full(yflat), p, cfg, env)[rows].ravel()

    sol = solve_ivp(f, t_span, y0[rows].ravel(), method=method, t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    out = np.broadcast_to(y0, (sol.t.size,) + shape).copy()
    out[:, rows] = sol.y.T.reshape((sol.t.size,) + active_shape)
    return sol.t, out
