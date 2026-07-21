import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def basic_heatmap(ds):
    fig, ax = plt.subplots(1,1)
    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000
    heatmaps, x_edges, y_edges = [], [], []
    
    # 1. Extract data
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id).compute()
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        
        heatmap, x_edge, y_edge = np.histogram2d(
            lon_raw, 
            lat_raw, 
            bins=[360, 180], 
           # range=[[-180, 180], [-90, 90]], # Hash out as needed to "zoom in"
            density=False
        )
        heatmaps.append(heatmap)
        x_edges.append(x_edge)
        y_edges.append(y_edge)
        xmin = min(x.min() for x in x_edges)
        ymin = min(y.min() for y in y_edges)
        xmax = max(x.max() for x in x_edges)
        ymax = max(y.max() for y in y_edges)
        
    # 2. Create heatmap
    heatmap_matrix = np.sum(heatmaps, axis=0)
    plt.imshow(
        heatmap_matrix.T,
        cmap = 'viridis',
        aspect = 'auto',
        interpolation = 'none',
        origin = 'lower',
        extent = [xmin, xmax, ymin, ymax]
    )
        
    # 3. Format heatmap layout
    plt.colorbar(label='Particle counts')
    plt.xlabel('Longitude / deg')
    plt.ylabel('Latitude / deg')
    plt.title('Potential trajectories through Venus cloud decks')
    plt.show()



def double_heatmap(ds):
    fig, ax = plt.subplots(1,2, figsize=(11, None), constrained_layout=True)
    t0 = ds.time.values.min()
    heatmaps_lat, xlat_edges, ylat_edges = [], [], []
    heatmaps_z, xz_edges, yz_edges = [], [], []
    
    # 1. Extract data
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id).compute()
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        
        heatmap_lat, xlat_edge, ylat_edge = np.histogram2d(
            lon_raw, 
            lat_raw, 
            bins=[360, 180], 
            range=[[-180, 180], [-90, 90]], # Hash out as needed to "zoom in"
            density=False
        )

        heatmap_z, xz_edge, yz_edge = np.histogram2d(
            lon_raw,
            z_raw,
            bins=[360, 150],
            range=[[-180, 180], [30, 70]],
            density=False
        )
        
        heatmaps_lat.append(heatmap_lat)
        xlat_edges.append(xlat_edge)
        ylat_edges.append(ylat_edge)
        xmin_lat = min(x.min() for x in xlat_edges)
        ymin_lat = min(y.min() for y in ylat_edges)
        xmax_lat = max(x.max() for x in xlat_edges)
        ymax_lat = max(y.max() for y in ylat_edges)

        heatmaps_z.append(heatmap_z)
        xz_edges.append(xz_edge)
        yz_edges.append(yz_edge)
        xmin_z = min(x.min() for x in xz_edges)
        ymin_z = min(y.min() for y in yz_edges)
        xmax_z = max(x.max() for x in xz_edges)
        ymax_z = max(y.max() for y in yz_edges)
        
        
    #2. Create heatmaps
    heatmap_latmatrix = np.sum(heatmaps_lat, axis=0)
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

    heatmap_zmatrix = np.sum(heatmaps_z, axis=0)
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

    # 3. Format graph layouts
    plt.colorbar(label='Particle counts')
    fig.suptitle('Potential trajectories through the Venus cloud decks')
    plt.show()