import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def basic_heatmap(ds):
    fig, ax = plt.subplots(1,1)
    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000
    heatmaps = []
    
    # 1. Extract data
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id).compute()
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        
        heatmap, x_edges, y_edges = np.histogram2d(
            lon_raw, 
            lat_raw, 
            bins=[360, 180], 
            density=False
        )
        heatmaps.append(heatmap)
        
    #2. Create heatmap
    heatmap_matrix = np.concatenate([heatmap for heatmap in heatmaps])
    plt.imshow(
        heatmap_matrix.T,
        cmap = 'viridis',
        aspect = "auto",
        interpolation = 'none',
        origin = 'lower',
        extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
    )
        
    # 2. Format graph layout
    plt.colorbar(label='Particle counts')
    plt.xlabel('Longitude / deg')
    plt.ylabel('Latitude / deg')
    plt.title('Potential trajectories through Venus cloud decks')
    plt.show()



def double_heatmap(ds):
    ## To be finished later! not currently in use
    fig, ax = plt.subplots(1,2)
    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000
    
    # 1. Extract data
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id).compute()
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        
        heatmap_lat, x_edges, y_edges = np.histogram2d(
            lon_raw, 
            lat_raw, 
            bins=200, 
            density=False
        )

        heatmap_z, x_edges, y_edges = np.histogram2d(
            lon_raw, 
            z_raw, 
            bins=[360, 70], 
            density=False
        )

        start_lon = traj.lon.values[0]
        start_lat = traj.lat.values[0]
        start_z   = traj.z.values[0]/1000
        
        #2. Create heatmap
        plt.subplot(1, 2, 1)
        plt.imshow(
            heatmap_lat.T,
            cmap = 'viridis',
            interpolation = 'none',
            extent = [-180, 180, -90, 90],
            aspect = "equal"
        )
        plt.xlabel('Longitude / deg')
        plt.ylabel('Latitude / deg')
        plt.title('Potential trajectories through Venus cloud decks')

        plt.subplot(1, 2, 2)
        plt.imshow(
            heatmap_z.T,
            cmap = 'viridis',
            interpolation = 'none',
            extent = [-180, 180, 0, 70],
            aspect = "equal"
        )
        plt.xlabel('Longitude / deg')
        plt.ylabel('Altitude / km')
        plt.title('Potential trajectories through Venus cloud decks')

        
    # 2. Format graph layout
    plt.colorbar(label='Particle counts')
    plt.show()