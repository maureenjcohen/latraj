import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def basic_heatmap(ds):
    fig, ax = plt.subplots(1,1)
    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000
    
    # 1. Extract data
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id).compute()
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        
        lon_diffs = np.abs(np.diff(lon_raw))
        jump_indices = np.where(lon_diffs > 180)[0] + 1
        elapsed_days_raw = (traj.time.values - t0) / np.timedelta64(1, 'D')
        
        lon_clean = np.insert(lon_raw, jump_indices, np.nan)
        lat_clean = np.insert(lat_raw, jump_indices, np.nan)
        z_clean = np.insert(z_raw, jump_indices, np.nan)
        days_clean = np.insert(elapsed_days_raw, jump_indices, np.nan)
        heatmap, x_edges, y_edges = np.histogram2d(lon_raw, lat_raw, bins = [360, 180])

        start_lon = traj.lon.values[0]
        start_lat = traj.lat.values[0]
        start_z   = traj.z.values[0]/1000
        
        #2. Create heatmap
        plt.imshow(
            heatmap.T,
            cmap = 'viridis',
            colorizer = days_clean
        )
        
    # 2. Format graph layout
    plt.colorbar(label='Time (Days)')
    plt.xlabel('Longitude / deg')
    plt.ylabel('Latitude / deg')
    plt.title('Potential trajectories through Venus cloud decks')
    plt.show()