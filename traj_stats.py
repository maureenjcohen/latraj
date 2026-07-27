import xarray as xr
import numpy as np

# maybe-temporary config settings for upper and lower boundaries
lat_lowbound = 0
lat_highbound = 40
z_lowbound = 48
z_highbound = 55

def sum_ranges(ds):
    ds = ds.compute()
    lat_last, z_last = [], []
    lat_inrange, lat_outofrange = 0, 0
    z_inrange, z_outofrange = 0, 0
    
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000

        lat_last.append(lat_raw[-1])
        z_last.append(z_raw[-1])

    for value in lat_last:
        if value > lat_highbound or value < lat_lowbound:
            lat_outofrange += 1
        elif value < lat_highbound and value > lat_lowbound:
            lat_inrange += 1

    for value in z_last:
        if value > z_highbound or value < z_lowbound:
            z_outofrange += 1
        elif value < z_highbound and value > z_lowbound: 
            z_inrange += 1

    total = len(lat_last)
    lat_inpercent = (lat_inrange/total)*100
    lat_outpercent = (lat_outofrange/total)*100

    z_inpercent = (z_inrange/total)*100
    z_outpercent = (z_outofrange/total)*100
    
    print(f"{lat_inpercent}% of particles stayed in the latitude range.")
    print(f"{lat_outpercent}% of particles ended up outside the latitude range.")
    print(f"{z_inpercent}% of particles stayed in the altitude range.")
    print(f"{z_outpercent}% of particles ended up outside the altitude range.")
        

def means(ds):
    ds = ds.compute()
    lat_last, z_last = [], []
    lat_means, z_means = [], []
    lat_stdvs, z_stdvs = [], []
    x, y = 0, 0

    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000

        lat_mean = np.mean(lat_raw)
        lat_stdv = np.std(lat_raw)
        z_mean = np.mean(z_raw)
        z_stdv = np.std(z_raw)

        lat_last.append(lat_raw[-1])
        z_last.append(z_raw[-1])

        lat_means.append(lat_mean)
        lat_stdvs.append(lat_stdv)
        z_means.append(z_mean)
        z_stdvs.append(z_stdv)

    mean_lat_last = np.mean(lat_last)
    stdv_lat_last = np.std(lat_last)
    mean_z_last  = np.mean(z_last)
    stdv_z_last = np.std(z_last)

    print(f"Mean final position: Latitude {mean_lat_last} degrees, altitude {mean_z_last} km.")
    print(f"Standard deviation on final position: Latitude {stdv_lat_last} degrees, altitude {stdv_z_last} km.")
    
    print("\nMean values:")
    for lat_mean, z_mean in zip(lat_means, z_means):
        x += 1
        print(f"Trajectory: {x}, Mean Latitude: {lat_mean} degrees, Mean Altitude: {z_mean} km.")
        
    print("\nStandard Deviations:")
    for lat_stdv, z_stdv in zip(lat_stdvs, z_stdvs):
        y += 1
        print(f"Trajectory: {y}, Latitude: {lat_stdv} degrees, Altitude: {z_stdv} km.")

        
        