import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def basic_heatmap(ds):
    """ Plots a heatmap of latitude vs longitude for all trajectories.
    Includes a colorbar showing the particle count in each bin."""
    fig, ax = plt.subplots(1,1)
    ds = ds.compute()
    
    # 1. Extract data for each trajectory
    traj = ds.sel(trajectory=traj_id)
    lon_raw = ds.lon.values.ravel()
    lat_raw = ds.lat.values.ravel()
    z_raw = ds.z.values.ravel()/1000
    stuck = ds.stuck.values.ravel()

    # 2. Apply boolean masks
    mask = (stuck == 0.0) & np.isfinite(lon_raw) & np.isfinite(lat_raw)
    lon_raw = lon_raw[mask]
    lat_raw = lat_raw[mask]
    z_raw = z_raw[mask]

    # 3. Calculate heatmap histogram
    heatmap_matrix, x_edge, y_edge = np.histogram2d(lon_raw, lat_raw, bins=[360, 180], range=[[-180, 180], -[90, 90]], density=False)
        
    # 4. Set heatmap edges
    xmin, xmax = x_edge[0], x_edge[-1]
    ymin, ymax = y_edge[0], y_edge[-1]

    # 5. Plot heatmap
    plt.imshow(
        heatmap_matrix.T,
        cmap = 'viridis',
        aspect = 'auto',
        interpolation = 'none',
        origin = 'lower',
        extent = [xmin, xmax, ymin, ymax]
    )
        
    # 6. Format heatmap layout and show plot
    plt.colorbar(label='Particle count')
    plt.xlabel('Longitude / deg')
    plt.ylabel('Latitude / deg')
    plt.title('Potential trajectories through Venus cloud decks')
    plt.show()

    
def double_heatmap(ds):
    """Plots two heatmaps for all trajectories, one showing latitude vs
    longitude and the other showing altitude vs longitude. Both heatmaps have
    colorbars showing the particle count per bin, though only one is titled."""
    fig, ax = plt.subplots(1,2, figsize=(11, 6), constrained_layout=True)
    ds = ds.compute()
    
    # 1. Retrieve data for each trajectory
    lon_raw = ds.lon.values.ravel()
    lat_raw = ds.lat.values.ravel()
    z_raw = ds.z.values.ravel()/1000
    stuck = ds.stuck.values.ravel()

    # 2. Apply boolean masks
    mask = (stuck == 0.0) & np.isfinite(lon_raw) & np.isfinite(lat_raw) & np.isfinite(z_raw)
    lon_raw = lon_raw[mask]
    lat_raw = lat_raw[mask]
    z_raw = z_raw[mask]

    # 3. Calculate heatmap histograms
    heatmap_latmatrix, xlat_edge, ylat_edge = np.histogram2d(
        lon_raw, 
        lat_raw, 
        bins=[360, 180], 
        range=[[-180, 180], [-90, 90]], # Hash out as needed to "zoom in"
        density=False
    )

    heatmap_zmatrix, xz_edge, yz_edge = np.histogram2d(
        lon_raw,
        z_raw,
        bins=[360, 150],
        range=[[-180, 180], [0, 80]], # Hash out as needed to "zoom in"
        density=False
    )
        
    # 3. Set heatmap edges
    xmin_lat, xmax_lat = xlat_edge[0], xlat_edge[-1]
    ymin_lat, ymax_lat = ylat_edge[0], ylat_edge[-1]

    xmin_z, xmax_z = xz_edge[0], xz_edge[-1]
    ymin_z, ymax_z = yz_edge[0], yz_edge[-1]

    # 4. Plot heatmaps
    plt.subplot(1, 2, 1)
    plt.imshow(
        heatmap_latmatrix.T,
        cmap = 'viridis',
        aspect = "auto",
        interpolation = 'none',
        origin = 'lower',
        extent = [xmin_lat, xmax_lat, ymin_lat, ymax_lat]
    )
    plt.xlabel('Longitude / deg')
    plt.ylabel('Latitude / deg')
    plt.colorbar(label='Particle count', orientation='horizontal')

    plt.subplot(1, 2, 2)
    plt.imshow(
        heatmap_zmatrix.T,
        cmap = 'viridis',
        aspect="auto",
        interpolation = 'none',
        origin = 'lower',
        extent = [xmin_z, xmax_z, ymin_z, ymax_z]
    )
    plt.xlabel('Longitude / deg')
    plt.ylabel('Altitude / km')
    plt.colorbar(label='Particle count', orientation='horizontal')

    # 5. Add heatmap title and show plots
    fig.suptitle('Potential trajectories through the Venus cloud decks')
    plt.show()