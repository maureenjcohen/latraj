# %%
import xarray as xr
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from parcels import read_particlefile
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import seaborn as sns

# %%
def load_particlefile(path):
    df = read_particlefile(path).to_pandas()
    df = df.rename(columns={"x":"lon", "y":"lat", "t":"time",
                            "particle_id": "trajectory"})
    df = df.sort_values(["trajectory","time"])
    df["obs"] = df.groupby("trajectory").cumcount()
    return df.set_index(["trajectory","obs"]).to_xarray()

# %%
def traj3d(ds, traj_id):
    t0 = ds.time.values.min()
    z_min = ds.z.values.min()/1000
    z_max = ds.z.values.max()/1000
    fig = go.Figure()
        
# 1. Isolate and compute the data for this specific trajectory
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

    # 2. Add the 3D line to the plot
    fig.add_trace(go.Scatter3d(
        x=lon_clean,
        y=lat_clean,
        z=z_clean,
        mode='lines',
        name=f'Trajectory',
        hoverinfo='skip', # Disable hover functionality for the static plot
        
        line=dict(
            width=5,
            color=days_clean, 
            colorscale='plasma_r',               
            showscale=True,
            colorbar=dict(
                title="Time<br>(days)", 
                thickness=15, 
                len=0.6, 
                x=0.66,
                tickfont=dict(size=12) # Ensure labels are readable in print
            )
        )
    ))

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
    ))

    # 3. Define the static camera viewing angle
    # You may need to tweak the 'eye' coordinates to get the perfect perspective
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.3, y=-1.3, z=0.5) 
    )

    # 4. Format the 3D environment for a print document
    fig.update_layout(
        title=dict(
            text="Trajectory",
            x=0.5, 
            y=0.65,
            font=dict(size=24, family="Arial") # Use standard document fonts
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
            zaxis=dict(backgroundcolor="white", gridcolor="lightgrey", range=[z_min, z_max]),
        ),
        margin=dict(l=0, r=0, b=0, t=100), 
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
def compare_traj(ds1, ds2, traj_id):
    t0 = ds1.time.values.min()
    z_min = min(ds1.z.values.min()/1000, ds2.z.values.min()/1000)
    z_max = max(ds1.z.values.max()/1000,ds2.z.values.max()/1000)
    fig = go.Figure()
        
    label = ["Control", "Assim"]
    #style = ["lines", "markers"]
    symbols = ["circle", "x"]
# 1. Isolate and compute the data for this specific trajectory
    for i, ds in enumerate([ds1, ds2]):
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

        # 2. Add the 3D line to the plot
        if symbols[i] == "circle":
            fig.add_trace(go.Scatter3d(
                x=lon_clean,
                y=lat_clean,
                z=z_clean,
                mode="markers",
                name=f'{label[i]}',
                hoverinfo='skip', # Disable hover functionality for the static plot
                showlegend=False,

                marker=dict(
                    symbol=symbols[i],
                    size=4,        
                    color=days_clean,        
                    colorscale='plasma_r',   # choose a colorscale
                    opacity=0.8,
                    line=dict(width=2, color="Black"),
                    colorbar=dict(
                         title="Time<br>(days)", 
                         thickness=15, 
                         len=0.6, 
                         x=0.92,
                         tickfont=dict(size=12) # Ensure labels are readable in print
                     )
                )
            ))
        elif symbols[i] == "x":
                fig.add_trace(go.Scatter3d(
                x=lon_clean,
                y=lat_clean,
                z=z_clean,
                mode="markers",
                name=f'{label[i]}',
                hoverinfo='skip', # Disable hover functionality for the static plot
                showlegend=False,
                
                marker=dict(
                    symbol=symbols[i],
                    size=2,        
                    color=days_clean,        
                    colorscale='plasma_r',   # choose a colorscale
                    opacity=0.8,
                    line=dict(width=2, color="Black")
                )
            ))
                        
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
        ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="markers",
        marker=dict(symbol="circle", size=8, color="black"),
        name="Control",
        showlegend=True,
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None],
        mode="markers",
        marker=dict(symbol="x", size=8, color="black"),
        name="Assim",
        showlegend=True,
    ))
    # 3. Define the static camera viewing angle
    # You may need to tweak the 'eye' coordinates to get the perfect perspective
    zoom=1.3
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.3*zoom, y=-1.3*zoom, z=1.5*zoom) 
    )

    # 4. Format the 3D environment for a print document
    fig.update_layout(
        title=dict(
            text=f"Start coords: {float(start_lon):.2f} lon, {float(start_lat):.2f} lat, {float(start_z):.2f} km",
            x=0.5, 
            y=0.85,
            font=dict(size=12, family="Arial") # Use standard document fonts
        ),
        scene=dict(
            xaxis_title="Longitude / deg",
            yaxis_title="Latitude / deg",
            zaxis_title="Altitude / km",
            
            aspectmode='manual',
            aspectratio=dict(x=2, y=1, z=0.5),
            camera=camera,
            
            # Pure white backgrounds are best for document integration
            bgcolor='white', 
            xaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
            yaxis=dict(backgroundcolor="white", gridcolor="lightgrey"),
            zaxis=dict(backgroundcolor="white", gridcolor="lightgrey", range=[z_min, z_max]),
        ),
        margin=dict(l=0, r=0, b=0, t=10), 
        paper_bgcolor='white',
        plot_bgcolor='white',
        xaxis=dict(visible=False),
        yaxis=dict(visible=False)
    )

    # 5. Export as a high-resolution static image
    # Note: This requires the 'kaleido' package installed in your Python environment
    #fig.write_image(f"compare_traj_{traj_id}.png", width=1200, height=800, scale=3)

    # You can still call fig.show() in your notebook just to preview the camera angle
    fig.show()

# %%
def _edges_from_centres(c):
    """Cell edges as midpoints between (ascending) centres, extrapolated at the ends."""
    c = np.asarray(c, dtype=float)
    mids = 0.5 * (c[:-1] + c[1:])
    first = c[0] - (mids[0] - c[0])
    last  = c[-1] + (c[-1] - mids[-1])
    return np.concatenate([[first], mids, [last]])

# %%
def count_heatmap(ds1, ds2, z_range, lat_range=(-90,-30),
                  grid_path='/exomars/projects/mc5526/lagrangian_trajectory/Mars_inputs/control/control.nc', 
                  lon_step=16, lat_step=4, cmap="cividis", scale="log",
                  vmin=None,
                  save=False, savepath='/exomars/projects/mc5526/lagrangian_trajectory/scratch_plots/'):
    z_low, z_high = z_range
    labels = ["Control", "Assim"]

    # 1. Native grid centres -> ascending edges (histogram2d needs ascending bins)
    grid = xr.open_dataset(grid_path)
    lon_c = grid.lon.values                       # already ascending
    lat_c = grid.lat.values
    lat_order = np.argsort(lat_c)                        # descending -> ascending
    lat_asc = lat_c[lat_order]

    lon_edges = _edges_from_centres(lon_c)
    lat_edges = _edges_from_centres(lat_asc)

    # 2. Bin every observation into a native gridbox
    grids = []
    for ds in (ds1, ds2):
        lon = ds.lon.values.ravel()
        lat = ds.lat.values.ravel()
        z   = ds.z.values.ravel() / 1000.0              # m -> km
        mask = ((z >= z_low) & (z <= z_high)
                & np.isfinite(lon) & np.isfinite(lat))
        counts, _, _ = np.histogram2d(lon[mask], lat[mask],
                                    bins=[lon_edges, lat_edges])
        grids.append(counts.T[::-1])                     # rows=lat, north on top

    lat_disp = lat_asc[::-1]
    if lat_range is not None:
        lo, hi = min(lat_range), max(lat_range)
        keep = (lat_disp >= lo) & (lat_disp <= hi)   # boolean over rows
        grids    = [g[keep] for g in grids]
        lat_disp = lat_disp[keep]

    vmax = max(g.max() for g in grids)

    # 3. Thinned tick labels (north-on-top ordering for lat)
    xlabels = [f"{v:g}" if i % lon_step == 0 else "" for i, v in enumerate(lon_c)]
    ylabels = [f"{v:g}" if i % lat_step == 0 else "" for i, v in enumerate(lat_disp)]

    if scale == "log":
        norm = LogNorm(vmin=vmin if vmin is not None else 1, vmax=vmax)   # 0-count cells masked (see set_bad)
    elif scale == "linear":
        norm = Normalize(vmin=vmin if vmin is not None else 0, vmax=vmax)  # 0-count cells at bottom of ramp
    else:
        raise ValueError(f"scale must be 'log' or 'linear', got {scale!r}")
    lat_edges_disp = _edges_from_centres(lat_disp)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(cmap_obj(0.0))

    # 4. Draw the two heatmaps side by side
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, grid_counts, label in zip(axes, grids, labels):
        im = ax.pcolormesh(lon_edges, lat_edges_disp, grid_counts,
                            cmap=cmap_obj, norm=norm)
        ax.set_xticks(lon_c[::lon_step])
        ax.set_xticklabels([f"{v:.0f}" for v in lon_c[::lon_step]])
        ax.set_yticks(lat_disp[::lat_step])
        ax.set_yticklabels([f"{v:.0f}" for v in lat_disp[::lat_step]])
        fig.colorbar(im, ax=ax, label="Particle count")
        ax.set_title(f"{label}: {z_low:g}–{z_high:g} km")
        ax.set_xlabel("Longitude / deg")
        ax.set_ylabel("Latitude / deg")

    plt.tight_layout()
    if save:
        plt.savefig(savepath + f'heatmap_{z_range[0]}_to_{z_range[1]}.png', format='png', bbox_inches='tight')
    plt.show()
    return
# %%
