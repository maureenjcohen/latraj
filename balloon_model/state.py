""" Layout of the balloon state vector.

    The state is a float64 array of shape (N_STATE, n_particles). The layout
    holds every prognostic variable of the full model, so adding physics never
    changes it; variables a configuration does not evolve keep a zero tendency.
"""
import numpy as np

# Order of the state variables. Units: z [m], w, u, v [m/s], masses [kg],
# temperatures [K].
#   z, w        altitude and vertical velocity of the balloon
#   u, v        horizontal balloon velocity (evolved only with horizontal drag-lag)
#   m_ZP, m_SP  lifting-gas mass in the zero-pressure and superpressure balloons
#   T_ZP, T_SP  lifting-gas temperatures
#   T_e1..T_e4  envelope node temperatures (one node uses T_e1 only)
#   u_conv      convective gust variable [dimensionless], a unit-variance
#               Ornstein-Uhlenbeck process; the gust added to the vertical
#               wind is sigma * taper(z) * u_conv. Its statistics are those of
#               the wind seen by a balloon drifting with the flow, so it is
#               carried with the balloon rather than sampled from the air.
STATE_VARS = ("z", "w", "u", "v", "m_ZP", "m_SP", "T_ZP", "T_SP",
              "T_e1", "T_e2", "T_e3", "T_e4", "u_conv")
IDX = {name: i for i, name in enumerate(STATE_VARS)}
N_STATE = len(STATE_VARS)


def pack(n=None, **values):
    """ Build a state array from named values (scalars or length-n arrays).
        Variables not given are set to zero. n defaults to the length of the
        longest array given, or 1. """
    sizes = [np.size(v) for v in values.values() if np.ndim(v) > 0]
    if n is None:
        n = max(sizes) if sizes else 1
    unknown = set(values) - set(IDX)
    if unknown:
        raise KeyError(f"unknown state variables: {sorted(unknown)}")
    y = np.zeros((N_STATE, n), dtype=np.float64)
    for name, v in values.items():
        y[IDX[name]] = v
    return y


def get(y, name):
    """ Row of the state array for one variable (a view, not a copy). """
    return y[IDX[name]]
