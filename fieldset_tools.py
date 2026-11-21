""" FieldSet helpers shared by the Parcels v4 run scripts """
# %%
import numpy as np
import xarray as xr
import parcels._sgrid as sgrid
from parcels import FieldSet, SphericalMesh

# %%
def agrid_sgrid_metadata():
    """ SGRID metadata for a rectilinear A-grid on (lon, lat, height) node
        dimensions, with lon/lat as node coordinates. All variables sit on the
        nodes, so the face-dimension entries only satisfy the SGRID schema. """
    return sgrid.SGrid2DMetadata(
        cf_role="grid_topology",
        topology_dimension=2,
        node_dimensions=("lon", "lat"),
        node_coordinates=("lon", "lat"),
        face_dimensions=(
            sgrid.FaceNodePadding("lon_c", "lon", sgrid.Padding.LOW),
            sgrid.FaceNodePadding("lat_c", "lat", sgrid.Padding.LOW),
            ),
        vertical_dimensions=(sgrid.FaceNodePadding("height_c", "height", sgrid.Padding.LOW),),
        )

# %%
def add_zonal_halo(ds, n_cells):
    """ Add halo in longitude so that advection can wrap around.
        Parcels4 currently doesn't wrap.
        Hopefully temporary and v4 won't need this in its final release. """
    lon = ds.lon.values
    west = ds.isel(lon=slice(-n_cells, None)).assign_coords(lon=lon[-n_cells:] - 360.0)
    east = ds.isel(lon=slice(0, n_cells)).assign_coords(lon=lon[:n_cells] + 360.0)
    out = xr.concat([west, ds, east], dim="lon")
    out["lon"].attrs = ds["lon"].attrs # assign_coords drops the axis attribute that SGRID validation requires
    return out

# %%
def add_cell_areas(fieldset, radius):
    """ Return the fieldset combined with a cell_areas field (needed for smagdiff),
        computed from the lat/lon grid of U and the planetary radius [m].
        Parcels v4 FieldSets are combined with +, so the caller must use the
        returned FieldSet. """
    lat = np.asarray(fieldset.U.grid.lat)
    lon = np.asarray(fieldset.U.grid.lon)
    xlon, ylat = np.meshgrid(lon, lat)
    dlat = np.deg2rad(np.abs(np.gradient(ylat, axis=0)))
    dlon = np.deg2rad(np.abs(np.gradient(xlon, axis=1)))
    dy = dlat*radius
    dx = dlon*radius*np.abs(np.cos(np.deg2rad(ylat)))
    area = (dy*dx).astype("float32")

    ca = xr.Dataset(
        {"cell_areas": (["time", "height","lat","lon"], area[None, None])},
        coords={
            "lon": ("lon", lon.astype("float32"), {"axis": "X"}),
            "lat": ("lat", lat.astype("float32"), {"axis": "Y"}),
            "depth": ("height", np.array([0.0], dtype="float32"), {"axis": "Z"}),
            "time": ("time", np.array([fieldset.time_interval.left]), {"axis":"T"})
        },
    )
    ca = sgrid._attach_sgrid_metadata(ca, agrid_sgrid_metadata())
    ca_fs = FieldSet.from_sgrid_conventions(ca, mesh=SphericalMesh(radius=radius))
    return fieldset + ca_fs
