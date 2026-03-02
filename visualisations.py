# %%
import xarray as xr
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Assuming your dataset is loaded as 'ds'
# ds = xr.open_zarr('path_to_your_file.zarr')


# %%
def traj3d(ds):
    t0 = ds.time.values.min()
    fig = go.Figure()

    # Loop through each trajectory ID
    for i, traj_id in enumerate(ds.trajectory.values):
        
        # 1. Isolate and compute the data for this specific trajectory
        traj = ds.sel(trajectory=traj_id).compute()
        elapsed_days = (traj.time.values - t0) / np.timedelta64(1, 'D')

        # 2. Add the 3D line to the plot
        fig.add_trace(go.Scatter3d(
            x=traj.lon.values,
            y=traj.lat.values,
            z=traj.z.values*1e-3,
            mode='lines',
            hoverinfo='skip', # Disable hover functionality for the static plot
            
            line=dict(
                width=6,
                color=elapsed_days, 
                colorscale='Plasma',               
                showscale=True if i == 0 else False,
                colorbar=dict(
                    title="Time<br>(days)", 
                    thickness=15, 
                    len=0.6, 
                    x=0.7,
                    tickfont=dict(size=12) # Ensure labels are readable in print
                ) if i == 0 else None
            )
        ))

    # 3. Define the static camera viewing angle
    # You may need to tweak the 'eye' coordinates to get the perfect perspective
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.4, y=1.4, z=0.4) 
    )

    # 4. Format the 3D environment for a print document
    fig.update_layout(
        title=dict(
            text="Potential trajectories through the Venus cloud decks",
            x=0.5, 
            y=0.7,
            font=dict(size=20, family="Arial") # Use standard document fonts
        ),
        scene=dict(
            xaxis_title="Longitude / deg",
            yaxis_title="Latitude / deg",
            zaxis_title="Altitude / km",
            
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=0.5),
            camera=camera,
            
            # Pure white backgrounds are best for document integration
            bgcolor='white', 
            xaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
            yaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
            zaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
        ),
        margin=dict(l=0, r=0, b=0, t=60), 
        paper_bgcolor='white',
        plot_bgcolor='white',
        showlegend=False
    )

    # 5. Export as a high-resolution static image
    # Note: This requires the 'kaleido' package installed in your Python environment
    #fig.write_image("trajectory_grant_figure.png", width=1200, height=800, scale=3)

    # You can still call fig.show() in your notebook just to preview the camera angle
    fig.show()

# %%
def doubletraj(ds):
    fig = make_subplots(
    rows=1, cols=2,
    specs=[[{'type': 'scene'}, {'type': 'scene'}]],
    horizontal_spacing=0.0 # Brings the two plots slightly closer together
    )

    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000


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

        start_lon = traj.lon.values[0]
        start_lat = traj.lat.values[0]
        start_z   = traj.z.values[0]/1000

        # 3. Add the trace to the specific subplot column (i + 1)
        fig.add_trace(go.Scatter3d(
            x=lon_clean,
            y=lat_clean,
            z=z_clean,
            mode='lines',
            name=f'Trajectory {traj_id-7}',
            hoverinfo='skip',
            
            line=dict(
                width=5, # Slightly thinner lines look better in smaller subplots
                color=days_clean, 
                colorscale='plasma_r',
                
                # We only want one colorbar, so we attach it to the second plot 
                # (which sits on the right side of the figure)
                showscale=True if i == 1 else False,
                colorbar=dict(
                    title="Time<br>(Days)", 
                    thickness=15, 
                    len=0.6, 
                    x=0.46, # Pushed just outside the rightmost plot
                    tickfont=dict(size=12) 
                ) if i == 1 else None
            )
        ), row=1, col=i+1) # <-- Crucial: Tells Plotly which subplot to put this in

        fig.add_trace(go.Scatter3d(
        x=[start_lon],
        y=[start_lat],
        z=[start_z],
        mode='markers',
        marker=dict(
            size=6,         # Adjust size as needed for visibility
            color='red', 
            symbol='circle'
        ),
        showlegend=False,   # Keep the legend clean
        hoverinfo='skip'    # Disable hover since it's a static document
    ), row=1, col=i+1)      # <-- Crucial: Ensures the dot goes to the correct subplot

    # 4. Define the shared camera angle
    shared_camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.3, y=-1.3, z=0.5) 
    )

    # 5. Create a shared scene layout dictionary to avoid repeating code
    shared_scene_layout = dict(
        xaxis_title="Longitude / deg",
        yaxis_title="Latitude / deg",
        zaxis_title="Altitude / km",
        aspectmode='manual',
        aspectratio=dict(x=1, y=1, z=0.5),
        camera=shared_camera, # Apply the exact same camera angle to both
        bgcolor='white', 
        xaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
        yaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
        zaxis=dict(backgroundcolor="white", gridcolor="lightgrey", range=[z_min, z_max]),
    )

    # 6. Apply the formatting to both scenes (scene1 and scene2)
    fig.update_layout(
        title=dict(
        text="Potential trajectories through the Venus cloud decks",
        font=dict(size=24, family="Arial"),
        x=0.5,  # Centers the text horizontally across the whole figure
        y=0.7  # Pushes it up to the very top edge
    ),
        showlegend=False,
        scene=dict(
        **shared_scene_layout, 
        domain=dict(x=[0.0, 0.54], y=[0.0, 1.0])
    ),   # Applies to the left subplot
        scene2=dict(
        **shared_scene_layout, 
        domain=dict(x=[0.42, 1.0], y=[0.0, 1.0])
    ), # Applies to the right subplot
        
        margin=dict(l=0, r=0, b=0, t=100),
        paper_bgcolor='white',
        plot_bgcolor='white',
        
        # Optional: Format the subplot titles
        font=dict(family="Arial")
    )

    # 7. Export the wide figure
    # Use a wider aspect ratio (e.g., 1400x700) to accommodate side-by-side plots nicely
    #fig.write_image("trajectory_subplots.png", width=1400, height=700, scale=3)

    fig.show()
# %%
