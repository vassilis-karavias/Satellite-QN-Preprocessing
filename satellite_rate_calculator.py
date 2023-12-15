import pandas as pd
import os
import numpy as np

####################################### Local Data ####################################################################
# site_height = {"London": 11, "Madrid": 660, "Edinburgh": 47, "Athens": 20, "Barcelona": 12, "Paris": 35, "Nantes": 20,
#                "Lyon": 237, "Lisbon": 2, "Bern": 540, "Freiburg": 280, "Milan": 120, "Florence": 50, "Naples": 100,
#                "Frankfurt": 110, "Dortmund":86 , "Berlin": 35, "Ioannina": 480}

# time_zone = {"London": 0, "Madrid": 1, "Edinburgh": 0, "Athens": 2, "Barcelona": 1, "Paris": 1, "Nantes": 1,
#                "Lyon": 1, "Lisbon": 1, "Bern": 1, "Freiburg": 1, "Milan": 1, "Florence": 1, "Naples": 1,
#                "Frankfurt": 1, "Dortmund":1 , "Berlin": 1, "Ioannina": 2}
##################################### Global Data #####################################################################
site_height = {"London": 11, "Frankfurt": 112, "Graz": 490, "Johanessburg": 1753, "Sao Paulo": 760, "Tokyo": 40, "Auckland": 200,
                    "New Delhi": 216, "Mumbai": 14, "Bangalore": 920,
                    "Xinglong": 890, "Nanshan": 2028, "Perth": 0, "Brisbane": 32, "Melborne": 31, "Albany": 45, "Denver": 1609, "Nashville": 169}

time_zone = {"London": 0, "Frankfurt": 1, "Graz": 1, "Johanessburg": 2, "Sao Paulo": -3, "Tokyo": 9, "Auckland": 13,
                    "New Delhi": 5.5, "Mumbai": 5.5, "Bangalore": 5.5,
                    "Xinglong": 8, "Nanshan": 6, "Perth": 8, "Brisbane": 10, "Melborne": 11, "Albany": -4, "Denver": -6, "Nashville": -4}

def get_elevation_rates(file_paths):
    heights = ["0.0", "0.05", "0.1", "0.2", "0.3", "0.5"]
    elevation_heights = {}
    for h in heights:
        elevation_df = pd.read_csv(file_paths + h + "km.csv")
        elevation_heights[float(h)] = elevation_df
    return elevation_heights


def get_time_elevation_angle(file_path, satellite_site, satellite_RAAN):
    time_angle = pd.read_csv(file_path + "_" + satellite_site + "_" + str(satellite_RAAN) + ".csv")
    return time_angle


def simple_rates_calculator(elevation_heights, file_path, satellite_site, satellite_RAAN):
    time_angle_df = get_time_elevation_angle(file_path, satellite_site, satellite_RAAN)
    elevation_df = elevation_heights[min(list(elevation_heights.keys()), key = lambda x: abs(x * 1000-site_height[satellite_site]))]
    time_prev = None
    key_size_tot = 0.0
    for ind in time_angle_df.index:
        time = time_angle_df["currentTime"][ind]
        if time_prev == None:
            time_prev = time
            continue
        time_step = time - time_prev
        if time_step > 1000:
            time_prev = None
            continue
        elevation = time_angle_df["elevationAngle"][ind]
        elevation_current = elevation_df.iloc[(elevation_df["Elevation (deg)"] - elevation).abs().argsort()[:1]]
        key_size_tot += time_step * elevation_current["Rate (b)"].values[0]
        time_prev = time
    key_rate = key_size_tot / (24 * 60 * 60)
    return key_rate


def simple_clash_method_rates(elevation_heights, file_path, satellite_sites, satellite_RAAN):
    time_angle_dfs = {}
    times = []
    for site in satellite_sites:
        time_angle_df = get_time_elevation_angle(file_path, site, satellite_RAAN)
        time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
        times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
        times = list(set(times) | set(times_current))
        time_angle_dfs[site] = time_angle_df_non_zero_capacity
    times = sorted(times)
    number_at_times = {}
    for time in times:
        num = 0
        for site in time_angle_dfs.keys():
            if time_angle_dfs[site]["currentTime"].isin([time]).any():
                num += 1
        number_at_times[time] = num
    key_rates_sites = {}
    for site in time_angle_dfs.keys():
        elevation_df = elevation_heights[
            min(list(elevation_heights.keys()), key=lambda x: abs(x * 1000 - site_height[site]))]
        time_prev = None
        key_size_tot = 0.0
        for ind in time_angle_dfs[site].index:
            time = time_angle_dfs[site]["currentTime"][ind]
            number_at_time = number_at_times[time]
            if time_prev == None:
                time_prev = time
                continue
            time_step = time - time_prev
            if time_step > 1000:
                time_prev = None
                continue
            elevation = time_angle_dfs[site]["elevationAngle"][ind]
            elevation_current = elevation_df.iloc[(elevation_df["Elevation (deg)"] - elevation).abs().argsort()[:1]]
            key_size_tot += time_step * elevation_current["Rate (b)"].values[0]/ number_at_time
            time_prev = time
        key_rate = key_size_tot / (24 * 60 * 60)
        key_rates_sites[site] = key_rate
    return key_rates_sites



def simple_clash_method_rates_segragation(loss_data, file_path, satellite_sites, satellite_RAAN):
    time_angle_dfs = {}
    times = []
    for site in satellite_sites:
        time_angle_df = get_time_elevation_angle(file_path, site, satellite_RAAN)
        time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 8]
        times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
        times = list(set(times) | set(times_current))
        time_angle_dfs[site] = time_angle_df_non_zero_capacity
    times = sorted(times)
    number_at_times = {}
    for time in times:
        num = 0
        for site in time_angle_dfs.keys():
            if time_angle_dfs[site]["currentTime"].isin([time]).any():
                num += 1
        number_at_times[time] = num
    segragation = []
    for i in range(len(loss_data)):
        t = loss_data[i][0]
        elev_angle = loss_data[i][1]
        if elev_angle * 180 / np.pi >= 8:
            segragation.append((i,1/number_at_times[t]))
    return segragation


def simple_clash_method_rates_2(elevation_heights, file_path, satellite_sites, satellite_RAAN):
    time_angle_dfs = {}
    times = []
    for site in satellite_sites:
        time_angle_df = get_time_elevation_angle(file_path, site, satellite_RAAN)
        time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
        times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
        times = list(set(times) | set(times_current))
        time_angle_dfs[site] = time_angle_df_non_zero_capacity
    times = sorted(times)
    # create a timeset
    number_at_times = {}
    time_list_prev = []
    prev_time = None
    prev_num = 0
    for time in times:
        num = 0
        for site in time_angle_dfs.keys():
            if time_angle_dfs[site]["currentTime"].isin([time]).any():
                num += 1
        if prev_time == None:
            prev_time = time
            prev_num = num
        if time - prev_time > 100:
            for t in time_list_prev:
                number_at_times[t] = prev_num
            time_list_prev = [time]
            prev_time = time
            prev_num = num
        else:
            time_list_prev.append(time)
            prev_time = time
            if num > prev_num:
                prev_num = num
    for t in time_list_prev:
        number_at_times[t] = prev_num
    key_rates_sites = {}
    for site in time_angle_dfs.keys():
        elevation_df = elevation_heights[
            min(list(elevation_heights.keys()), key=lambda x: abs(x * 1000 - site_height[site]))]
        time_prev = None
        key_size_tot = 0.0
        for ind in time_angle_dfs[site].index:
            time = time_angle_dfs[site]["currentTime"][ind]
            number_at_time = number_at_times[time]
            if time_prev == None:
                time_prev = time
                continue
            time_step = time - time_prev
            if time_step > 1000:
                time_prev = None
                continue
            elevation = time_angle_dfs[site]["elevationAngle"][ind]
            elevation_current = elevation_df.iloc[(elevation_df["Elevation (deg)"] - elevation).abs().argsort()[:1]]
            key_size_tot += time_step * elevation_current["Rate (b)"].values[0]/ number_at_time
            time_prev = time
        key_rate = key_size_tot / (24 * 60 * 60)
        key_rates_sites[site] = key_rate
    return key_rates_sites


def simple_clash_method_rates_2_segragation(loss_data, file_path, satellite_sites, satellite_RAAN, ground_station, n =1):
    time_angle_dfs = {}
    times = []
    for site in satellite_sites:
        time_angle_df = get_time_elevation_angle(file_path, site, satellite_RAAN)
        time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
        times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
        times = list(set(times) | set(times_current))
        time_angle_dfs[site] = time_angle_df_non_zero_capacity
    times = sorted(times)
    # create a timeset
    number_at_times = {}
    time_list_prev = []
    prev_time = None
    prev_num = 0
    for time in times:
        num = 0
        for site in time_angle_dfs.keys():
            if time_angle_dfs[site]["currentTime"].isin([time]).any():
                num += 1
        if prev_time == None:
            prev_time = time
            prev_num = num
        if time - prev_time > 100:
            for t in time_list_prev:
                number_at_times[t] = prev_num
            time_list_prev = [time]
            prev_time = time
            prev_num = num
        else:
            time_list_prev.append(time)
            prev_time = time
            if num > prev_num:
                prev_num = num
    for t in time_list_prev:
        number_at_times[t] = prev_num
    segragation = []
    prev_t = None
    prev_i = None
    time_zone_dark_low = 1800 + time_zone[ground_station]  * 3600
    time_zone_dark_high = 45000 + time_zone[ground_station] * 3600

    for i in range(len(loss_data)):
        t = loss_data[i][0]
        elev_angle = loss_data[i][1]
        if elev_angle * 180 / np.pi >= 0.2:
            if prev_t == None:
                prev_t = t
                prev_i = i
                continue
            if abs(t - prev_t) > 2:
                if time_zone_dark_low < 0:
                    if 0<prev_t<time_zone_dark_high or 86400 + time_zone_dark_low < prev_t< 86400:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
                    else:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))
                elif time_zone_dark_high > 86400:
                    if 0<prev_t<time_zone_dark_high- 86400 or time_zone_dark_low < prev_t< 86400:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
                    else:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))
                else:
                    if prev_t <45000 + time_zone_dark_high and prev_t > 1800 + time_zone_dark_low:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
                    else:
                        segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))
                prev_t = None
                prev_i = None
                continue
            prev_t = t
            prev_i = i
    if prev_t != None:
        if time_zone_dark_low < 0:
            if 0 < prev_t < time_zone_dark_high or 86400 + time_zone_dark_low < prev_t < 86400:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
            else:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))
        elif time_zone_dark_high > 86400:
            if 0 < prev_t < time_zone_dark_high - 86400 or time_zone_dark_low < prev_t < 86400:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
            else:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))
        else:
            if prev_t < time_zone_dark_high and prev_t > time_zone_dark_low:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "dark"))
            else:
                segragation.append((prev_i, min(n/number_at_times[prev_t], 1), "light"))

    return segragation



def no_clash_method_rates_segragation(loss_data, file_path, satellite_sites, satellite_RAAN, ground_station):
    time_angle_dfs = {}
    times = []
    for site in satellite_sites:
        time_angle_df = get_time_elevation_angle(file_path, site, satellite_RAAN)
        time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
        times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
        times = list(set(times) | set(times_current))
        time_angle_dfs[site] = time_angle_df_non_zero_capacity
    times = sorted(times)
    # create a timeset
    segragation = []
    prev_t = None
    prev_i = None
    time_zone_dark_low = 1800 + time_zone[ground_station] * 3600
    time_zone_dark_high = 45000 + time_zone[ground_station] * 3600
    for i in range(len(loss_data)):
        t = loss_data[i][0]
        elev_angle = loss_data[i][1]
        if elev_angle * 180 / np.pi >= 0.2:
            if prev_t == None:
                prev_t = t
                prev_i = i
                continue
            if abs(t - prev_t) > 2:
                if time_zone_dark_low < 0:
                    if 0<prev_t<time_zone_dark_high or 86400 + time_zone_dark_low < prev_t< 86400:
                        segragation.append((prev_i, 1, "dark"))
                    else:
                        segragation.append((prev_i, 1, "light"))
                elif time_zone_dark_high > 86400:
                    if 0<prev_t<time_zone_dark_high- 86400 or time_zone_dark_low < prev_t< 86400:
                        segragation.append((prev_i, 1, "dark"))
                    else:
                        segragation.append((prev_i, 1, "light"))
                else:
                    if prev_t <time_zone_dark_high and prev_t > time_zone_dark_low:
                        segragation.append((prev_i, 1, "dark"))
                    else:
                        segragation.append((prev_i, 1, "light"))
                prev_t = None
                prev_i = None
                continue
            prev_t = t
            prev_i = i
    if prev_t != None:
        if time_zone_dark_low < 0:
            if 0 < prev_t < time_zone_dark_high or 86400 + time_zone_dark_low < prev_t < 86400:
                segragation.append((prev_i, 1, "dark"))
            else:
                segragation.append((prev_i, 1, "light"))
        elif time_zone_dark_high > 86400:
            if 0 < prev_t < time_zone_dark_high - 86400 or time_zone_dark_low < prev_t < 86400:
                segragation.append((prev_i, 1, "dark"))
            else:
                segragation.append((prev_i, 1, "light"))
        else:
            if prev_t < time_zone_dark_high and prev_t > time_zone_dark_low:
                segragation.append((prev_i, 1, "dark"))
            else:
                segragation.append((prev_i, 1, "light"))
    return segragation



if __name__ == "__main__":
    elevation_heights = get_elevation_rates("elevation_rates_h_")
    RAAN = [80, 90, 100, 110 , 120, 130]
    for R in RAAN:
        # key_rate = simple_rates_calculator(elevation_heights, file_path = "satelliteElevAngle", satellite_site= "Paris", satellite_RAAN= R)
        # print("Rate for RAAN " + str(R) + ":" + str(key_rate))
        key_rate_sites = simple_clash_method_rates(elevation_heights, file_path="satelliteElevAngle",
                                                     satellite_sites=list(site_height.keys()), satellite_RAAN=R)
        for site in key_rate_sites.keys():
            print("Rate for RAAN " + str(R) + " for site " + site + ": " + str(key_rate_sites[site]))
        key_rate_sites = simple_clash_method_rates_2(elevation_heights, file_path = "satelliteElevAngle", satellite_sites = list(site_height.keys()), satellite_RAAN = R)
        for site in key_rate_sites.keys():
            print("Rate for RAAN " + str(R) +" for site " + site + ": " + str(key_rate_sites[site]))