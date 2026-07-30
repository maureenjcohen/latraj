import xarray as xr
import numpy as np
import datetime as dt
import config
from prettytable import PrettyTable

latitude_boundaries = {
    'South Poles': [-90, -59],
    'South Midlatitudes': [-60, -39],
    'South Equator': [-40, -1],
    'North Equator': [0, 39],
    'North Midlatitudes': [40, 59],
    'North Poles': [60, 90]
}

altitude_boundaries = {
    'Deep Atmosphere': [0, 48000],
    'Convective Clouds': [48001, 55000],
    'Upper Clouds': [55001, 70000]
}


def get_boundaries():
    """ Retrieves the upper and lower boundaries of the region 
    the particle was initially placed in by extracting the initial
    position from config.py and comparing it with dictionaries.
    Takes the first element of a list, assumes that all elements are identical.
    """
    lat_lowbound, lat_highbound, z_lowbound, z_highbound = None, None, None, None
    if isinstance(config.PARTICLE_LAT, (int, float)):
        config.PARTICLE_LAT = [config.PARTICLE_LAT]
    if isinstance(config.PARTICLE_DEPTH, (int, float)):
        config.PARTICLE_DEPTH = [config.PARTICLE_DEPTH]
    for region, coords in latitude_boundaries.items():
        if config.PARTICLE_LAT[0] >= coords[0] and config.PARTICLE_LAT[0] <= coords[1]:
            lat_lowbound, lat_highbound = coords[0], coords[1]
            if lat_lowbound is None or lat_highbound is None:
                raise ValueError('PARTICLE_LAT is outside all latitude regions')
    for region, coords in altitude_boundaries.items():
        if config.PARTICLE_DEPTH[0] >= coords[0] and config.PARTICLE_DEPTH[0] <= coords[1]:
            z_lowbound, z_highbound = coords[0], coords[1]
            if z_lowbound is None or z_highbound is None:
                raise ValueError('PARTICLE_DEPTH is outside all altitude regions')
    return lat_lowbound, lat_highbound, z_lowbound, z_highbound

    
def final_ranges(ds):
    """ Calculates the percentage of trajectories which end inside or outside
    of their latitude and altitude regions. """
    ds = ds.compute()
    lat_lowbound, lat_highbound, z_lowbound, z_highbound = get_boundaries()
    lat_last, z_last = [], []
    lat_inrange, lat_outofrange = 0, 0
    z_inrange, z_outofrange = 0, 0
    
    # 1. Retrieve data for each trajectory 
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lat_raw = traj.lat.values
        z_raw = traj.z.values
        stuck = traj.stuck.values

    # 2. Apply boolean masks to filter out trajectories which got stuck:
        mask = (stuck == 0.0) & np.isfinite(lat_raw) & np.isfinite(z_raw)
        lat_raw = lat_raw[mask]
        z_raw = z_raw[mask]
        if lat_raw.size == 0:
            continue

        lat_last.append(lat_raw[-1])
        z_last.append(z_raw[-1])
    
    #3. Calculate sums of trajectories ending in-range or out-of-range
    for value1, value2 in zip(lat_last, z_last):
        if value1 >= lat_highbound or value1 <= lat_lowbound:
            lat_outofrange += 1
        elif value1 <= lat_highbound and value1 >= lat_lowbound:
            lat_inrange += 1

        if value2 >= z_highbound or value2 <= z_lowbound:
            z_outofrange += 1
        elif value2 <= z_highbound and value2 >= z_lowbound: 
            z_inrange += 1

    # 4. Calculate percentage of trajectories which ended in-range or out-of-range
    total = len(ds.trajectory.values)
    lat_inpercent = (lat_inrange/total)*100
    lat_outpercent = (lat_outofrange/total)*100
    z_inpercent = (z_inrange/total)*100
    z_outpercent = (z_outofrange/total)*100

    # 5. Print results
    print(f"Latitude range: ({lat_lowbound}, {lat_highbound}) degrees. Altitude range: ({z_lowbound/1000}, {z_highbound/1000}) km.")    
    print(f"{lat_inpercent}% of trajectories ended within the latitude range.")
    print(f"{lat_outpercent}% of trajectories ended outside the latitude range.")
    print(f"{z_inpercent}% of trajectories ended within the altitude range.")
    print(f"{z_outpercent}% of trajectories ended outside the altitude range.")


def range_count(ds):
    """ Calculates the percentage of trajectories which leave their latitude
    or altitude regions, and for how long. """
    ds = ds.compute()
    lat_lowbound, lat_highbound, z_lowbound, z_highbound = get_boundaries()
    lat_time_diffs, z_time_diffs = [], []
    lat_inrange, lat_outofrange = 0, 0
    z_inrange, z_outofrange = 0, 0    
    x = 0

    # 1. Retrieve data for each trajectory:
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        time = traj.time.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values
        stuck = traj.stuck.values

    # 2. Apply boolean masks to isolate values where the trajectory goes out-of-range:
        mask = (stuck == 0.0) & np.isfinite(lat_raw) & np.isfinite(z_raw)
        time = time[mask]
        dt_days = (time[1] - time[0]) / np.timedelta64(1, 'D')
        lat_mask = np.where((lat_lowbound >= lat_raw[mask]) | (lat_highbound <= lat_raw[mask]), True, False)
        z_mask = np.where((z_lowbound >= z_raw[mask]) | (z_highbound <= z_raw[mask]), True, False)
        lat_time, z_time = time[lat_mask], time[z_mask]

    # 3. Calculate the timedelta for trajectories which went out-of-range:     
        if not np.size(lat_time) == 0:
            lat_timedelta = np.size(lat_time)*dt_days
            lat_outofrange += 1
        else:
            lat_timedelta = 0
            lat_inrange += 1 

        if not np.size(z_time) == 0:
            z_timedelta = np.size(z_time)*dt_days
            z_outofrange += 1
        else:
            z_timedelta = 0
            z_inrange += 1 
        
        lat_time_diffs.append(lat_timedelta)
        z_time_diffs.append(z_timedelta)

    # 5. Calculate the percentage of trajectories which went out-of-range or stayed in-range:
    total = len(ds.trajectory.values)
    lat_goodpercent = (lat_inrange/total)*100
    lat_badpercent = (lat_outofrange/total)*100
    z_goodpercent = (z_inrange/total)*100
    z_badpercent = (z_outofrange/total)*100

    # 6. Print percentage results:
    print(f"Latitude range: ({lat_lowbound}, {lat_highbound}) degrees. Altitude range: ({z_lowbound/1000}, {z_highbound/1000}) km.")
    print(f"{lat_goodpercent}% of particles stayed within the latitude region. {lat_badpercent}% of particles left the latitude region.")
    print(f"{z_goodpercent}% of particles stayed within the altitude region. {z_badpercent}% of particles left the altitude region.")

    # 7. Create tables showing the amount of time a trajectory spent in-range or out-of-range:
    lat_table = PrettyTable(["Trajectory no.", "Days spent within latitude region", "Days spent out of latitude region"])
    z_table = PrettyTable(["Trajectory no.", "Days spent within altitude region", "Days spent out of altitude region"])
    
    for lat_timedelta, z_timedelta in zip(lat_time_diffs, z_time_diffs):
        x+=1
        lat_table.add_row([x, config.RUNTIME_DAYS-lat_timedelta, lat_timedelta])
        z_table.add_row([x, config.RUNTIME_DAYS-z_timedelta, z_timedelta])
        
    print(lat_table)
    print(z_table)


def ejection_count(ds):
    """ Calculates the percentage of trajectories which get ejected
    near the poles, and for how long. """
    ds = ds.compute()
    time_diffs = []
    stuck_sum = 0
    
    # 1. Retrieve data from each trajectory:
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        time = traj.time.values
        stuck = traj.stuck.values

    # 2. Apply boolean masks to isolate stuck values:
        mask = (stuck == 1)
        time = time[mask]

    # 3. Calculate timedelta for trajectories which got ejected:
        if not np.size(time) == 0:
            timedelta = (time[-1] - time[0]) / np.timedelta64(1, 'D')
            stuck_sum += 1
        else:
            timedelta = 0
        time_diffs.append((traj_id, timedelta))

    # 4. Calculate the percentage of trajectories which were ejected:
    stuck_percent = (stuck_sum/len(ds.trajectory.values))*100
    print(f"{stuck_percent}% of particles were ejected at the pole.")

    # 5. Create table showing the amount of time spent ejected:
    stuck_table = PrettyTable(["Trajectory no.", "Days spent within boundaries", "Days spent ejected"])
    for traj_id, timedelta in time_diffs:
        stuck_table.add_row([traj_id+1, config.RUNTIME_DAYS-timedelta, timedelta])
    print(stuck_table)


def means(ds):
    """ Calculates the mean and standard deviation for latitude and altitude
    for each trajectory, and the mean and standard deviation on their final
    values of latitude and altitude. """
    ds = ds.compute()
    lat_last, z_last = [], []
    lat_means, z_means = [], []
    lat_stdvs, z_stdvs = [], []
    x, y = 0, 0

    # 1. Retrieve data for each trajectory:
    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000
        stuck = traj.stuck.values

    # 2. Apply boolean masks to filter out trajectories which got stuck:
        mask = (stuck == 0.0) & np.isfinite(lat_raw) & np.isfinite(z_raw)
        lat_raw = lat_raw[mask]
        z_raw = z_raw[mask]
        if lat_raw.size == 0:
            continue

    # 3. Calculate means & standard deviations for each trajectory:
        lat_mean = np.mean(lat_raw)
        lat_stdv = np.std(lat_raw)
        z_mean = np.mean(z_raw)
        z_stdv = np.std(z_raw)

        lat_means.append(lat_mean)
        lat_stdvs.append(lat_stdv)
        z_means.append(z_mean)
        z_stdvs.append(z_stdv)

        lat_last.append(lat_raw[-1])
        z_last.append(z_raw[-1])

    # 4. Calculate means & standard deviations for the final positions:
    mean_lat_last = np.mean(lat_last)
    stdv_lat_last = np.std(lat_last)
    mean_z_last  = np.mean(z_last)
    stdv_z_last = np.std(z_last)

    print(f"Mean final position: Latitude {mean_lat_last} degrees, altitude {mean_z_last} km.")
    print(f"Standard deviation on final position: Latitude {stdv_lat_last} degrees, altitude {stdv_z_last} km.")

    # 5. Create table showing the means and standard deviations:
    means_table = PrettyTable(["Trajectory no.", "Mean Latitude (deg)", "Mean Altitude (km)"])
    stdv_table = PrettyTable(["Trajectory no.", "Latitude Standard Deviation (deg)", "Altitude Standard Deviation (km)"])
    for lat_mean, z_mean in zip(lat_means, z_means):
        x += 1
        means_table.add_row([x, lat_mean, z_mean])
        
    for lat_stdv, z_stdv in zip(lat_stdvs, z_stdvs):
        y += 1
        stdv_table.add_row([y, lat_stdv, z_stdv])

    print(means_table)
    print(stdv_table)