""" Tests for the shared Parcels v4 FieldSet helpers in fieldset_tools.py """
import numpy as np
import xarray as xr
import pytest
import parcels._sgrid as sgrid
from parcels import FieldSet, SphericalMesh
from fieldset_tools import add_zonal_halo, add_cell_areas, agrid_sgrid_metadata


def lon_dataset(lon):
    """ Minimal dataset with a longitude axis carrying the SGRID axis attribute """
    lon = np.asarray(lon, dtype=np.float32)
    return xr.Dataset(
        {"U": (["lon"], np.arange(lon.size, dtype=np.float32))},
        coords={"lon": ("lon", lon, {"axis": "X", "units": "degrees_east"})},
    )


def test_zonal_halo_extends_and_wraps():
    """ The halo adds n_cells on each side with longitudes shifted by 360 and the
        data copied from the opposite edge; the axis stays strictly increasing. """
    lon = np.arange(-180.0, 180.0, 3.75)
    ds = add_zonal_halo(lon_dataset(lon), 3)
    out = ds.lon.values
    assert out.size == lon.size + 6
    assert np.all(np.diff(out) > 0)
    np.testing.assert_allclose(out[:3], lon[-3:] - 360.0)
    np.testing.assert_allclose(out[-3:], lon[:3] + 360.0)
    np.testing.assert_array_equal(ds.U.values[:3], np.arange(lon.size)[-3:])
    np.testing.assert_array_equal(ds.U.values[-3:], np.arange(3))


def test_zonal_halo_keeps_axis_attribute():
    """ SGRID validation needs axis='X' on lon, which xr.concat would otherwise drop """
    ds = add_zonal_halo(lon_dataset(np.arange(0.0, 360.0, 3.75)), 3)
    assert ds.lon.attrs["axis"] == "X"


@pytest.mark.parametrize("lon", [np.linspace(-180.0, 180.0, 97), np.linspace(0.0, 360.0, 97)])
def test_zonal_halo_rejects_repeated_endpoint(lon):
    """ An axis that includes both cyclic endpoints would duplicate the seam """
    with pytest.raises(ValueError, match="repeats the cyclic endpoint"):
        add_zonal_halo(lon_dataset(lon), 3)


def test_cell_areas_match_analytic_and_preserve_time_interval():
    """ cell_areas equals R^2 dlat dlon cos(lat) on a uniform grid, is the same at
        every altitude, and combining it does not change the FieldSet time interval. """
    radius = 6.0518e6
    lon = np.arange(0.0, 360.0, 3.75, dtype=np.float32)
    lat = np.arange(-60.0, 61.0, 2.0, dtype=np.float32)
    alt = np.array([0.0, 50000.0, 90000.0], dtype=np.float32)
    shape = (2, alt.size, lat.size, lon.size)
    times = np.datetime64("1987-03-30T00:00:00", "ns") + np.arange(2) * np.timedelta64(1, "h")
    ds = xr.Dataset(
        {name: (["time", "height", "lat", "lon"], np.zeros(shape, dtype=np.float32)) for name in "UVW"},
        coords={
            "lon": ("lon", lon, {"axis": "X", "units": "degrees_east"}),
            "lat": ("lat", lat, {"axis": "Y", "units": "degrees_north"}),
            "depth": ("height", alt, {"axis": "Z", "units": "m", "positive": "up"}),
            "time": ("time", times, {"axis": "T"}),
        },
    )
    ds = sgrid._attach_sgrid_metadata(ds, agrid_sgrid_metadata())
    fieldset = FieldSet.from_sgrid_conventions(ds, mesh=SphericalMesh(radius=radius))
    interval = fieldset.time_interval
    fieldset = add_cell_areas(fieldset, radius)
    assert fieldset.time_interval == interval

    y = np.array([0.0, 30.0, -50.0])
    x = np.array([90.0, 180.0, 270.0])
    expected = (np.deg2rad(2.0) * radius) * (np.deg2rad(3.75) * radius * np.cos(np.deg2rad(y)))
    for z in alt:
        area = fieldset.cell_areas.eval(np.full(3, 1800.0), np.full(3, z), y, x)
        np.testing.assert_allclose(area, expected, rtol=1e-5)
