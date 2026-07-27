import xarray as xr
import numpy as np
import datetime as dt
import config

latitude_boundaries = {
    'South Poles': [-90, -60],
    'South Midlatitudes': [-60, -40],
    'South Equator': [-40, 0],
    'North Equator': [0, 40],
    'North Midlatitudes': [40, 60],
    'North Poles': [60, 90]
}

altitude_boundaries = {
    'Deep Atmosphere': [0, 48000],
    'Convective Clouds': [48000, 55000],
    'Upper Clouds': [55000, 70000]
}

def get_boundaries():
    lat_lowbound, lat_highbound, z_lowbound, z_highbound = 0, 0, 0, 0
    for region, coords in latitude_boundaries.items():
        if config.PARTICLE_LAT[0] >= coords[0] and config.PARTICLE_LAT[0] <= coords[1]:
            lat_lowbound, lat_highbound = coords[0], coords[1]
    for region, coords in altitude_boundaries.items():
        if config.PARTICLE_DEPTH[0] >= coords[0] and config.PARTICLE_DEPTH[0] <= coords[1]:
            z_lowbound, z_highbound = coords[0], coords[1]
    return lat_lowbound, lat_highbound, z_lowbound, z_highbound
    raise TypeError('Does not accept ints or floats')
    raise ValueError('Input value too high or too low')

def sum_ranges(ds):
    ds = ds.compute()
    lat_last, z_last = [], []
    lat_inrange, lat_outofrange = 0, 0
    z_inrange, z_outofrange = 0, 0
    lat_lowbound, lat_highbound, z_lowbound, z_highbound = get_boundaries()
    print(f"Latitude range: ({lat_lowbound}, {lat_highbound}) degrees, altitude range: ({z_lowbound/1000}, {z_highbound/1000}) km.")
    
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values
        stuck = traj.stuck.values

        mask = stuck == 0.0
        lon_raw = lon_raw[mask]
        lat_raw = lat_raw[mask]
        z_raw = z_raw[mask]

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


def ejection_count(ds):
    ds = ds.compute()
    stuck_values = []
    time_diffs = []
    stuck_sum = 0
    x = 0
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        time = traj.time.values
        stuck = traj.stuck.values
        stuck_values.append(stuck[0])
        
        mask = stuck == 1
        time = time[mask]
        
        if not np.size(time) == 0:
            timedelta = time[-1] - time[0]
            timedelta = int(timedelta) / (24*3600000000000)
            time_diffs.append(timedelta)
        else:
            timedelta = 0
            time_diffs.append(timedelta)
        
        if 1 in stuck:
            stuck_sum += 1
            
    total = len(stuck_values)
    stuck_percent = (stuck_sum/total)*100
    print(f"{stuck_percent}% of particles were ejected at the pole.")
    for timedelta in time_diffs:
        x+=1
        print(f"Trajectory {x}: Particle trajectory spent {60-timedelta} days in range, {timedelta} days out of range.")


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
        stuck = traj.stuck.values

        mask = stuck == 0.0
        lon_raw = lon_raw[mask]
        lat_raw = lat_raw[mask]
        z_raw = z_raw[mask]
        
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

            
        