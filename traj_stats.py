import xarray as xr
import numpy as np

# maybe-temporary config settings for upper and lower range boundaries
lat_lowbound = 40
lat_highbound = 60
z_lowbound = 0
z_highbound = 48

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
    """
    for value1, value2 in lat_last, z_last:
        if value1 < 40 and value1 > -40 and value2 < 48:
            in_range += 1
        elif (value1 > 40 or value1 < -40) and value2 > 48:
            out_range += 1
    """
    total = len(lat_last)
    lat_inpercent = (lat_inrange/total)*100
    lat_outpercent = (lat_outofrange/total)*100

    z_inpercent = (z_inrange/total)*100
    z_outpercent = (z_outofrange/total)*100
    
    print(f"{lat_inpercent}% of particles stayed in the latitude range.")
    print(f"{lat_outpercent}% of particles ended up outside the latitude range.")
    print(f"{z_inpercent}% of particles stayed in the altitude range.")
    print(f"{z_outpercent}% of particles ended up outside the altitude range.")
    #print(f'{in_range} particles stayed completely in range.')
    #print(f'{out_range} particles went completely outside of the range.')
        



def mean(ds):
    ds = ds.compute()

    for i, traj_id in enumerate(ds.trajectory.values):
        traj = ds.sel(trajectory=traj_id)
        lon_raw = traj.lon.values
        lat_raw = traj.lat.values
        z_raw = traj.z.values/1000

        lat_mean = np.mean(lat_raw)
        lat_stdv = np.std(lat_raw)
        z_mean = np.mean(z_raw)
        z_stdv = np.std(z_raw)

        
        