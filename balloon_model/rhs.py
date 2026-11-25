""" Right-hand side of the balloon equations.

    Functions take and return arrays over particles. p is a BalloonParams,
    cfg a ModelConfig.
"""
import numpy as np
from . import state


def total_mass(m_ZP, m_SP, p):
    """ Total mass [kg]: gas from the state, payload and envelopes from p [E-DYN-1]. """
    return p.m_struct + m_ZP + m_SP


def displaced_volume(m_ZP, rho_a, p, cfg):
    """ Volume of air displaced by the balloon [m^3] [E-GAS-1, E-GAS-2].

        Isothermal gas at ambient pressure has density rho_a M_gas / M_atm,
        so its volume follows from the ambient density alone; the ZP sealed
        cap limits it to V_infl. """
    if cfg.thermal_model != "isothermal":
        raise NotImplementedError(cfg.thermal_model)
    return np.minimum(m_ZP * (p.M_atm / p.M_gas) / rho_a, p.V_infl)


def vertical_acceleration(w, w_a, rho_a, V, m_tot, p, cfg):
    """ dw/dt [m s^-2] from the momentum equation [E-DYN-1, E-DRG-1]:

        (m_tot + C_m rho_a V) dw/dt = (rho_a V - m_tot) g - 0.5 rho_a C_D A_top |w_r| w_r,
        w_r = w - w_a. """
    if cfg.drag_model != "constant" or cfg.area_model != "fixed":
        raise NotImplementedError((cfg.drag_model, cfg.area_model))
    w_r = w - w_a
    buoyancy = (rho_a * V - m_tot) * p.g
    drag = 0.5 * rho_a * p.C_D_top * p.A_top * np.abs(w_r) * w_r
    return (buoyancy - drag) / (m_tot + p.C_m * rho_a * V)


def rhs(t, y, p, cfg, env):
    """ Time derivative of the state array y (shape (N_STATE, n)).
        Variables the configuration does not evolve have zero tendency. """
    z = state.get(y, "z")
    w = state.get(y, "w")
    m_ZP = state.get(y, "m_ZP")
    a = env.sample(t, z)
    V = displaced_volume(m_ZP, a["rho"], p, cfg)
    m_tot = total_mass(m_ZP, state.get(y, "m_SP"), p)
    dydt = np.zeros_like(y)
    dydt[state.IDX["z"]] = w
    dydt[state.IDX["w"]] = vertical_acceleration(w, a["w"], a["rho"], V, m_tot, p, cfg)
    return dydt
