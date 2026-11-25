""" Model parameters and configuration switches.

    BalloonParams holds every scalar parameter of the craft, its materials,
    the gases and the empirical coefficients. Any field may be a scalar or a
    per-particle array, so one call can carry an ensemble of designs.
    Parameters with no source value are None; physics that needs one of them
    raises instead of guessing a value.

    Temperature-dependent gas properties are functions, not parameters.
"""
from dataclasses import dataclass
import numpy as np

R_UNIVERSAL = 8.314462618   # Molar gas constant [J mol^-1 K^-1]
SIGMA_SB = 5.670374e-8      # Stefan-Boltzmann constant [W m^-2 K^-4]
M_HE = 4.0026e-3            # Molar mass of helium [kg mol^-1]


@dataclass
class BalloonParams:
    # --- A. Craft specification ---
    # m_payload makes the total mass 62 kg, the value used by the existing
    # balloon_vertical runs; the gondola and control module of the Izraelevitz
    # et al. (2024, J. Aircraft, preprint) prototype alone are 26.5 kg.
    m_payload: object = 34.95       # Gondola, payload and control module [kg]
    m_env_ZP: object = 21.3         # ZP envelope [kg]; Izraelevitz et al. (2024) give 21.3 kg for both envelopes and fittings together
    m_env_SP: object = 0.0          # SP envelope [kg]; included in m_env_ZP until the masses are known separately
    m_gas_ZP: object = 5.75         # Helium in the ZP balloon [kg]
    m_gas_SP: object = 0.0          # Helium in the SP balloon [kg]; no source value
    V_infl: object = 72.6           # ZP fully inflated volume [m^3]
    D_ZP_top: object = 5.0          # ZP upper sphere diameter [m]
    D_ZP_bot: object = 2.5          # ZP lower sphere diameter [m]
    D_SP: object = 2.5              # SP sphere diameter [m]
    A_top: object = 19.6            # Vertical drag reference area [m^2] (pi * 5^2 / 4)
    A_side: object = None           # Sideways drag reference area [m^2]; from the shape model
    dP_SP_max: object = 14.0e3      # SP structural pressure limit [Pa] (proof test)

    # --- B. Envelope material ---
    alpha_sol_ZP: object = 0.08     # ZP solar absorptance
    eps_IR_ZP: object = 0.52        # ZP infrared emissivity
    alpha_IR_ZP: object = None      # ZP infrared absorptance; None means equal to eps_IR_ZP (grey body)
    tau_sol: object = 0.0           # Solar transmittance (0 for a metallised film)
    r_sol: object = 0.0             # Solar reflectance of a transparent film
    tau_IR: object = 0.0            # Infrared transmittance
    r_IR: object = 0.0              # Infrared reflectance of a transparent film
    alpha_sol_SP: object = None     # SP solar absorptance; no source value
    eps_IR_SP: object = None        # SP infrared emissivity; no source value
    c_env: object = 1200.0          # Envelope specific heat [J kg^-1 K^-1]; a sweep default, swept over 1000-2300, not a measured value

    # --- C. Lifting gas ---
    M_gas: object = M_HE            # Molar mass of the lifting gas [kg mol^-1]

    # --- D. Ambient atmosphere ---
    M_atm: object = 43.45e-3        # Mean molar mass of the Venus atmosphere [kg mol^-1]
    g: object = 8.87                # Gravitational acceleration [m s^-2]

    # --- F. Empirical coefficients ---
    C_D_top: object = 0.8           # Vertical drag coefficient
    C_D_side: object = 1.0          # Sideways drag coefficient
    C_m: object = 0.2               # Virtual-mass coefficient
    C4: object = 1.5                # External convection multiplier
    C6: object = 3.0                # Internal convection multiplier
    C_d_vent: object = 0.6          # Vent discharge coefficient (sharp-edged orifice)

    # --- G. Gas transfer ---
    d_vent: object = 0.007          # SP-to-ZP vent orifice diameter [m]
    Vdot_pump: object = None        # Pump volumetric rate [m^3 s^-1]; no source value
    pump_adiabatic: bool = True     # Adiabatic heating of pumped gas

    @property
    def m_struct(self):
        """ Mass that is not lifting gas [kg]: payload and envelopes. The gas
            masses can change during a run, so the state carries them and
            m_gas_ZP, m_gas_SP here are initial values. """
        return self.m_payload + self.m_env_ZP + self.m_env_SP

    @property
    def R_gas(self):
        """ Specific gas constant of the lifting gas [J kg^-1 K^-1]. """
        return R_UNIVERSAL / self.M_gas

    @property
    def R_atm(self):
        """ Specific gas constant of the atmosphere [J kg^-1 K^-1]. """
        return R_UNIVERSAL / self.M_atm

    @property
    def cv_gas(self):
        """ Specific heat at constant volume of a monatomic gas [J kg^-1 K^-1]. """
        return 1.5 * self.R_gas

    @property
    def cp_gas(self):
        """ Specific heat at constant pressure of a monatomic gas [J kg^-1 K^-1]. """
        return 2.5 * self.R_gas

    @classmethod
    def balloon_vertical(cls, m_total, m_gas_ZP, **overrides):
        """ Parameters equivalent to the fieldset constants of the existing
            balloon_vertical kernel. That kernel takes the total mass as one
            number and the lift ratio as 10.86, so the payload absorbs the
            difference and M_atm is set to give the same ratio unless given. """
        if "m_payload" in overrides:
            raise TypeError("m_payload is derived from m_total; pass m_total instead")
        overrides.setdefault("M_atm", 10.86 * overrides.get("M_gas", M_HE))
        p = cls(m_gas_ZP=m_gas_ZP, **overrides)
        p.m_payload = m_total - (p.m_env_ZP + p.m_env_SP + p.m_gas_ZP + p.m_gas_SP)
        return p


# Allowed values of each switch; the first is the default. Values listed in
# PLANNED are part of the design but not implemented yet.
SWITCHES = {
    "thermal_model": ("isothermal",),
    "drag_model": ("constant",),
    "area_model": ("fixed",),
    "balloon_type": ("zero_pressure",),
    "zp_overpressure": ("sealed",),
    "horizontal": ("advect",),
}
PLANNED = {
    "thermal_model": ("one_node", "six_node"),
    "drag_model": ("Re",),
    "area_model": ("sphere", "shape"),
    "balloon_type": ("superpressure", "balloon_in_balloon"),
    "zp_overpressure": ("vent",),
    "horizontal": ("drag_lag",),
}


@dataclass
class ModelConfig:
    """ One switch per physics option; defaults give the minimal model, which
        reproduces the existing balloon_vertical kernel. """
    thermal_model: str = "isothermal"   # gas at ambient temperature and pressure
    drag_model: str = "constant"        # constant C_D_top
    area_model: str = "fixed"           # A_top fixed
    balloon_type: str = "zero_pressure" # single ZP balloon
    zp_overpressure: str = "sealed"     # gas mass fixed when the ZP is full
    horizontal: str = "advect"          # balloon moves with the horizontal wind

    def __post_init__(self):
        for name, allowed in SWITCHES.items():
            value = getattr(self, name)
            if value in PLANNED[name]:
                raise NotImplementedError(f"{name} = {value!r} is not implemented yet")
            if value not in allowed:
                raise ValueError(f"{name} = {value!r}; allowed: {allowed + PLANNED[name]}")
