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


@pytest.mark.parametrize("z", [49999.0, 54001.0, np.nan])
def test_profile_raises_outside(z):
    with pytest.raises(ValueError, match="outside"):
        profile_env().sample(0.0, np.array([z]))


def test_profile_rejects_bad_input():
    with pytest.raises(ValueError):
        ProfileAtmosphere([1.0, 0.0], rho=[1, 1], T=[1, 1], P=[1, 1], u=[0, 0], v=[0, 0], w=[0, 0])
    with pytest.raises(ValueError, match="missing"):
        ProfileAtmosphere([0.0, 1.0], rho=[1, 1])
    good = dict(rho=[1, 1], T=[1, 1], P=[1, 1], u=[0, 0], v=[0, 0], w=[0, 0])
    for z in ([], [0.0], [0.0, np.nan]):
        with pytest.raises(ValueError, match="z_levels"):
            ProfileAtmosphere(z, **good)
    with pytest.raises(ValueError, match="shape"):
        ProfileAtmosphere([0.0, 1.0], **{**good, "T": [1, 1, 1]})
    with pytest.raises(ValueError, match="non-finite"):
        ProfileAtmosphere([0.0, 1.0], **{**good, "rho": [1, np.nan]})


def test_profile_z_range():
    assert profile_env().z_range == (50000.0, 54000.0)


def column_dataset(z):
    """ Dataset in the create_nc_venus layout, two time levels, with fields
        linear in lon, lat and z (and +1 in every field at the second time),
        so any interpolated column is known exactly. """
    lon = np.array([0.0, 10.0, 20.0])
    lat = np.array([-5.0, 0.0, 5.0])
    Z, LAT, LON = np.meshgrid(z, lat, lon, indexing="ij")
    def field(a, b, c):
        f = a + b * LON + c * LAT + 1e-5 * Z
        return (["time", "height", "lat", "lon"], np.stack([f, f + 1.0]).astype(np.float32))
    return xr.Dataset(
        {"RHO": field(1.0, 0.001, 0.002), "TEMP": field(300.0, 0.1, 0.0), "PRES": field(9e4, 1.0, 0.0),
         "U": field(-50.0, 0.0, 0.0), "V": field(0.0, 0.0, 0.1), "W": field(0.0, 0.0, 0.0)},
        coords={"lon": lon, "lat": lat, "depth": ("height", z), "time": [0, 1]})


@pytest.mark.parametrize("time_index, offset", [(0, 0.0), (1, 1.0)])
def test_profile_from_dataset(time_index, offset):
    """ A column read from a create_nc_venus-layout dataset is bilinear in
        lon and lat, taken at the requested time, and ordered by altitude
        even when the file stores the levels top down. """
    z = np.array([52000.0, 48000.0, 50000.0])
    env = ProfileAtmosphere.from_dataset(column_dataset(z), lon=5.0, lat=2.5, time_index=time_index)
    z_sorted = np.sort(z)
    np.testing.assert_array_equal(env.z, z_sorted)
    expected_rho = 1.0 + 0.001 * 5.0 + 0.002 * 2.5 + 1e-5 * z_sorted + offset
    np.testing.assert_allclose(env.profiles["rho"], expected_rho, rtol=1e-6)
    np.testing.assert_allclose(env.profiles["v"], 0.25 + 1e-5 * z_sorted + offset, rtol=1e-6)


@pytest.mark.parametrize("lon, lat", [(270.0, 0.0), (5.0, 30.0)])
def test_profile_from_dataset_rejects_points_outside(lon, lat):
    with pytest.raises(ValueError, match="outside the dataset"):
        ProfileAtmosphere.from_dataset(column_dataset(np.array([48000.0, 50000.0])), lon=lon, lat=lat)
