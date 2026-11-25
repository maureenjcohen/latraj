import numpy as np
import pytest
import xarray as xr
from balloon_model.environment import IsothermalAtmosphere, UniformAtmosphere, ProfileAtmosphere
from balloon_model.params import R_UNIVERSAL

KEYS = {"rho", "T", "P", "u", "v", "w"}


def venus_isothermal(**winds):
    return IsothermalAtmosphere(T=300.0, P_ref=5.0e4, z_ref=55000.0, M_atm=43.45e-3, g=8.87, **winds)


def test_isothermal_is_hydrostatic_ideal_gas():
    """ dP/dz = -rho g (checked by a centred difference, truncation error
        ~ (dz/H)^2 / 6 ~ 4e-7 relative for dz = 10 m) and P = rho R T / M. """
    env = venus_isothermal()
    z = np.array([40000.0, 55000.0, 65000.0])
    a = env.sample(0.0, z)
    assert set(a) == KEYS
    dPdz = (env.sample(0.0, z + 10.0)["P"] - env.sample(0.0, z - 10.0)["P"]) / 20.0
    np.testing.assert_allclose(dPdz, -a["rho"] * env.g, rtol=1e-6)
    np.testing.assert_allclose(a["P"], a["rho"] * R_UNIVERSAL * a["T"] / env.M_atm, rtol=1e-14)


def test_isothermal_density_inverse():
    env = venus_isothermal()
    z = np.array([45000.0, 55000.0, 62000.0])
    np.testing.assert_allclose(env.altitude_of_density(env.density_at(z)), z, rtol=1e-14)


def test_sample_shapes_and_winds():
    env = venus_isothermal(u=-60.0, w=0.3)
    a = env.sample(0.0, np.full(4, 55000.0))
    for k in KEYS:
        assert a[k].shape == (4,)
    np.testing.assert_array_equal(a["u"], -60.0)
    np.testing.assert_array_equal(a["w"], 0.3)


def test_uniform_atmosphere():
    a = UniformAtmosphere(rho=0.9, T=300.0, P=5e4, w=0.1).sample(0.0, np.array([0.0, 1e5]))
    np.testing.assert_array_equal(a["rho"], 0.9)
    np.testing.assert_array_equal(a["w"], 0.1)


def profile_env():
    z = np.array([50000.0, 52000.0, 54000.0])
    return ProfileAtmosphere(z, rho=[1.2, 1.0, 0.9], T=[330.0, 320.0, 310.0], P=[1e5, 9e4, 8e4],
                             u=[0.0, 1.0, 2.0], v=[0.0, 0.0, 0.0], w=[0.1, 0.2, 0.3])


def test_profile_exact_at_levels_and_linear_between():
    env = profile_env()
    a = env.sample(0.0, np.array([50000.0, 52000.0, 54000.0, 51000.0, 53500.0]))
    np.testing.assert_allclose(a["rho"], [1.2, 1.0, 0.9, 1.1, 0.925], rtol=1e-14)
    np.testing.assert_allclose(a["w"], [0.1, 0.2, 0.3, 0.15, 0.275], rtol=1e-14)


@pytest.mark.parametrize("z", [49999.0, 54001.0])
def test_profile_raises_outside(z):
    with pytest.raises(ValueError, match="outside"):
        profile_env().sample(0.0, np.array([z]))


def test_profile_rejects_bad_input():
    with pytest.raises(ValueError):
        ProfileAtmosphere([1.0, 0.0], rho=[1, 1], T=[1, 1], P=[1, 1], u=[0, 0], v=[0, 0], w=[0, 0])
    with pytest.raises(ValueError, match="missing"):
        ProfileAtmosphere([0.0, 1.0], rho=[1, 1])


def test_profile_from_dataset():
    """ A column read from a create_nc_venus-layout dataset is bilinear in
        lon and lat and ordered by altitude. Fields are linear in lon, lat and
        z, so the interpolated column is known exactly. """
    lon = np.array([0.0, 10.0, 20.0])
    lat = np.array([-5.0, 0.0, 5.0])
    z = np.array([48000.0, 50000.0, 52000.0])
    Z, LAT, LON = np.meshgrid(z, lat, lon, indexing="ij")
    def field(a, b, c):
        return (["time", "height", "lat", "lon"], (a + b * LON + c * LAT + 1e-5 * Z)[None].astype(np.float32))
    ds = xr.Dataset(
        {"RHO": field(1.0, 0.001, 0.002), "TEMP": field(300.0, 0.1, 0.0), "PRES": field(9e4, 1.0, 0.0),
         "U": field(-50.0, 0.0, 0.0), "V": field(0.0, 0.0, 0.1), "W": field(0.0, 0.0, 0.0)},
        coords={"lon": lon, "lat": lat, "depth": ("height", z), "time": [0]})
    env = ProfileAtmosphere.from_dataset(ds, lon=5.0, lat=2.5)
    np.testing.assert_array_equal(env.z, z)
    expected_rho = 1.0 + 0.001 * 5.0 + 0.002 * 2.5 + 1e-5 * z
    np.testing.assert_allclose(env.profiles["rho"], expected_rho, rtol=1e-6)
    np.testing.assert_allclose(env.profiles["v"], 0.25 + 1e-5 * z, rtol=1e-6)
