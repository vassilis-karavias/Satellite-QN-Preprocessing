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



def get_time_elevation_angle(file_path, satellite_site, satellite_RAAN):
    time_angle = pd.read_csv(file_path + "_" + satellite_site + "_" + str(satellite_RAAN) + ".csv")
    return time_angle




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

