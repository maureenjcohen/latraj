""" Tests for the create_nc_venus preprocessing script.

    Regression context: the height axis was once written in km instead of the
    metres Parcels expects (the *1e3 conversion was silently lost), which made
    every particle seed depth out of bounds and surfaced only as cryptic
    out-of-bounds/zero-density errors deep inside kernel execution. These tests
    pin the axis units at the source.
"""
import sys
import types

import numpy as np
import netCDF4 as nc
import pytest

# create_nc_venus reads config at import time, but config.py is gitignored and
# machine-specific. Use the real one when present; otherwise install a minimal
# stub so the module can be imported on any checkout (e.g. CI).
try:
    import config  # noqa: F401
except ModuleNotFoundError:
    stub = types.ModuleType('config')
    stub.INPUT_FILE = '/nonexistent/input.nc'
    stub.EXPERIMENT_NAME = 'UnitTest'
    stub.OUTPUT_DIR = '/nonexistent/output'
    stub.T_SELECT = (0, None)
    stub.H_SELECT = (0, None)
    sys.modules['config'] = stub

import create_nc_venus


@pytest.mark.parametrize("n_levels", [50, 78])
def test_level_heights_are_metres(n_levels):
    """ Each VPCM level-height set must be in metres, ascending, and complete.
        A km-scale axis (top ~97/146 instead of ~97000/146000) fails the range
        check. """
    h = np.asarray(create_nc_venus.LEVEL_HEIGHTS[n_levels])
    assert h.ndim == 1 and len(h) == n_levels
    assert np.all(np.diff(h) > 0), "height axis must be strictly ascending"
    assert h[0] < 100.0, "lowest model level should be within ~100 m of the surface"
    assert 5.0e4 < h[-1] < 1.0e6, (
        f"top model level is {h[-1]:g}; expected ~1e5 m. "
        "A value below ~1000 means the axis is in km, not metres."
    )


def test_heights_for_selects_by_level_count():
    """ heights_for mirrors venuslab/venusdata.py: pick the height array by the
        input file's vertical level count, and refuse to guess otherwise. """
    assert create_nc_venus.heights_for(50) is create_nc_venus.heights50
    assert create_nc_venus.heights_for(78) is create_nc_venus.heights78
    with pytest.raises(ValueError, match="60 vertical levels"):
        create_nc_venus.heights_for(60)


def test_make_file_roundtrip(tmp_path):
    """ Write a small synthetic run through make_file and verify the file that
        Parcels will actually read: metre-scale ascending Height axis with the
        right metadata, and scalar fields (RHO) on the same dims as the winds. """
    hghts = create_nc_venus.heights78[::20]  # spread subset, still metres
    lats = np.linspace(-60.0, 60.0, 3)
    lons = np.linspace(-180.0, 180.0, 5)
    n_t, time_len = 2, 3600
    shape = (n_t, len(hghts), len(lats), len(lons))

    u = np.full(shape, 1.0, dtype=np.float32)
    v = np.full(shape, -1.0, dtype=np.float32)
    w = np.full(shape, 0.1, dtype=np.float32)
    rho = np.full(shape, 0.85, dtype=np.float32)

    path = tmp_path / 'roundtrip.nc'
    ncout = nc.Dataset(path, 'w', format='NETCDF4')
    create_nc_venus.make_file(ncout, u, v, w, hghts, lats, lons,
                              n_t, time_len,
                              scalars={'RHO': (rho, 'kg/m3')})
    ncout.close()

    ds = nc.Dataset(path)
    try:
        height = ds['Height']
        h = height[:]
        assert height.units == 'm'
        assert height.positive == 'up'
        assert np.all(np.diff(h) > 0)
        assert h[-1] > 1.0e4, (
            f"top of written Height axis is {h[-1]:g}; a metre-scale planetary "
            "atmosphere file should extend beyond 10 km"
        )
        np.testing.assert_allclose(h, hghts, rtol=1e-6)

        for name in ('U', 'V', 'W', 'RHO'):
            assert name in ds.variables, f"{name} missing from output file"
            assert ds[name].dimensions == ('time', 'height', 'lat', 'lon')
        assert ds['RHO'].units == 'kg/m3'
        np.testing.assert_allclose(ds['RHO'][:], rho, rtol=1e-6)

        assert len(ds['Time'][:]) == n_t
    finally:
        ds.close()
