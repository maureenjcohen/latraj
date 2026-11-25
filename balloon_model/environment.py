""" Atmosphere seen by the balloon.

    Every environment provides sample(t, z, lat=None, lon=None), returning a
    dict of float64 arrays with the shape of z:

        rho [kg m^-3], T [K], P [Pa], u, v, w [m s^-1]

    Later physics adds radiative keys (SW_dn, SW_up, SW_dir, LW_dn, LW_up,
    SZA) to the same dict. Column environments ignore lat and lon. The
    physics only ever calls sample(), so environments are interchangeable.
"""
import numpy as np
from .params import R_UNIVERSAL


def _full(value, z):
    return np.broadcast_to(np.asarray(value, dtype=np.float64), np.shape(z)).copy()


class IsothermalAtmosphere:
    """ Hydrostatic isothermal layer [E-DYN-2]: P = P_ref exp(-(z - z_ref)/H),
        H = R T / (M_a g), rho = P M_a / (R T). Uniform winds. """

    def __init__(self, T, P_ref, z_ref, M_atm, g, u=0.0, v=0.0, w=0.0):
        self.T, self.P_ref, self.z_ref = T, P_ref, z_ref
        self.M_atm, self.g = M_atm, g
        self.u, self.v, self.w = u, v, w
        self.H = R_UNIVERSAL * T / (M_atm * g)   # Scale height [m]

    def density_at(self, z):
        return self.P_ref * np.exp(-(np.asarray(z) - self.z_ref) / self.H) * self.M_atm / (R_UNIVERSAL * self.T)

    def altitude_of_density(self, rho):
        """ Altitude at which the density equals rho (inverse of density_at). """
        rho_ref = self.P_ref * self.M_atm / (R_UNIVERSAL * self.T)
        return self.z_ref + self.H * np.log(rho_ref / np.asarray(rho))

    def sample(self, t, z, lat=None, lon=None):
        z = np.asarray(z, dtype=np.float64)
        P = self.P_ref * np.exp(-(z - self.z_ref) / self.H)
        return {"rho": P * self.M_atm / (R_UNIVERSAL * self.T), "T": _full(self.T, z), "P": P,
                "u": _full(self.u, z), "v": _full(self.v, z), "w": _full(self.w, z)}


class UniformAtmosphere:
    """ Constant density, temperature, pressure and winds. Not hydrostatic:
        a test atmosphere in which a balloon of fixed volume has constant
        buoyancy at every altitude. """

    def __init__(self, rho, T, P, u=0.0, v=0.0, w=0.0):
        self.values = {"rho": rho, "T": T, "P": P, "u": u, "v": v, "w": w}

    def sample(self, t, z, lat=None, lon=None):
        return {k: _full(v, z) for k, v in self.values.items()}


class ProfileAtmosphere:
    """ Time-constant vertical profile, linearly interpolated in z (the same
        interpolation Parcels applies between model levels). Sampling outside
        the profile raises instead of extrapolating. """

    KEYS = ("rho", "T", "P", "u", "v", "w")

    def __init__(self, z_levels, **profiles):
        self.z = np.asarray(z_levels, dtype=np.float64)
        if np.any(np.diff(self.z) <= 0):
            raise ValueError("z_levels must increase strictly")
        missing = set(self.KEYS) - set(profiles)
        if missing:
            raise ValueError(f"profiles missing: {sorted(missing)}")
        self.profiles = {k: np.asarray(profiles[k], dtype=np.float64) for k in self.KEYS}

    def sample(self, t, z, lat=None, lon=None):
        z = np.asarray(z, dtype=np.float64)
        if np.any((z < self.z[0]) | (z > self.z[-1])):
            raise ValueError(f"altitude outside the profile [{self.z[0]}, {self.z[-1]}] m")
        return {k: np.interp(z, self.z, prof) for k, prof in self.profiles.items()}

    @classmethod
    def from_dataset(cls, ds, lon, lat, time_index=0):
        """ Column at (lon, lat) from a dataset written by create_nc_venus:
            variables U, V, W, RHO, TEMP, PRES on (time, height, lat, lon) and
            the altitude coordinate depth [m] on height. The column is linearly
            interpolated in lon and lat, as Parcels does on that grid. """
        col = ds.isel(time=time_index).interp(lon=lon, lat=lat)
        z = np.asarray(col["depth"].values, dtype=np.float64)
        names = {"rho": "RHO", "T": "TEMP", "P": "PRES", "u": "U", "v": "V", "w": "W"}
        profiles = {k: np.asarray(col[v].values, dtype=np.float64) for k, v in names.items()}
        order = np.argsort(z)
        return cls(z[order], **{k: p[order] for k, p in profiles.items()})
