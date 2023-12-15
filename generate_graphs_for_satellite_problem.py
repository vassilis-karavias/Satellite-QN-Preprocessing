import numpy as np
import pandas as pd
from graph_tool.all import *
import satellites_graphs.satellite_rate_calculator as satellite_rate_calculator
from satellites_graphs.SatQuMA.channel.time_dependent_loss import (get_losses, get_losses_input_file)
from satellites_graphs.SatQuMA.key.get_key import SKL_sys_loop
import copy

########################################################## Local Data #################################################################

# time_zone = {"London": 0, "Madrid": 1, "Edinburgh": 0, "Athens": 2, "Barcelona": 1, "Paris": 1, "Nantes": 1,
#                "Lyon": 1, "Lisbon": 1, "Bern": 1, "Freiburg": 1, "Milan": 1, "Florence": 1, "Naples": 1,
#                "Frankfurt": 1, "Dortmund":1 , "Berlin": 1, "Ioannina": 2}
               
########################################################## Global Data ################################################################

time_zone = {"London": 0, "Frankfurt": 1, "Graz": 1, "Johanessburg": 2, "Sao Paulo": -3, "Tokyo": 9, "Auckland": 13,
                    "New Delhi": 5.5, "Mumbai": 5.5, "Bangalore": 5.5,
                    "Xinglong": 8, "Nanshan": 6, "Perth": 8, "Brisbane": 10, "Melborne": 11, "Albany": -4, "Denver": -6, "Nashville": -4}

class Network_Structure():

    def __init__(self, ground_core_sites, ground_station_sites, Tij):
        """
        Initialise the Network Structure problem. This needs a set of ground core sites,
        a set of ground station sites and a location-to-location traffic matrix. A list of ground station
        sites possible for each ground core network can either be specified or calculated based on strategies.

        :param ground_core_sites: A dictionary of ground core sites with {site: [lat, long, [ground station
        sites that the graph connects the site to]]} : |L| length
        :param ground_station_sites: A dictionary of the ground station site {site: [lat, long]}:
        :param Tij: Trasmission requirement Matrix between ground core sites: dict with {(site1,site2) : Tij}
        """
        # inputs of L, S, Tij. P is given by the satellites derived in calculation and Rij calculated in
        # satellite_rate_calculator. delta_ij calculated in this, as is K
        self.ground_core_sites = ground_core_sites
        self.ground_station_sites = ground_station_sites
        for site1, site2 in Tij.keys():
            if site1 not in self.ground_core_sites.keys() or site2 not in self.ground_core_sites.keys():
                print("Tij contains sites not in core sites")
                raise ValueError
        self.Tij = Tij


    def great_circle_distance_calculation(self, lat_1, lat_2, lon_1, lon_2):
        # Use Haversine formula to calculate the distance between the sites
        hav_lat = np.power(np.sin((lat_2 - lat_1)/2),2)
        hav_lon = np.power(np.sin((lon_2 - lon_1) / 2), 2)
        hav_dist = hav_lat + np.cos(lat_1) * np.cos(lat_2) * hav_lon
        if abs(hav_dist) <= 1:
            return 2 * 6371000 * np.arcsin(np.sqrt(hav_dist))
        else:
            return 2 * 6371000 * np.arcsin(1)

    def calculate_physical_distance_between_core_sites_and_stations(self):
        distance_dict = {}
        for core_site in self.ground_core_sites.keys():
            for ground_station in self.ground_station_sites.keys():
                lat_1 = self.ground_core_sites[core_site][0]
                lon_1 = self.ground_core_sites[core_site][1]
                lat_2 = self.ground_station_sites[ground_station][0]
                lon_2 = self.ground_station_sites[ground_station][1]
                dist = self.great_circle_distance_calculation(lat_1 =lat_1, lat_2 = lat_2, lon_1 = lon_1, lon_2 = lon_2)
                distance_dict[core_site, ground_station] = dist
        self.distance_dict = distance_dict


    def calculate_cost_core_network_graph(self, cost_per_km):
        try:
            distances_array = []
            for core, ground in self.distance_dict.keys():
                if ground in self.available_sites[core]:
                    distances_array.append(self.distance_dict[core, ground])
            distance_core_sites_array = []
            disance_core_sites_dict = {}
            for core in self.ground_core_sites.keys():
                ground = self.available_sites[core][0]
                for core_2 in self.ground_core_sites.keys():
                    if core_2 != core and self.available_sites[core_2][0] == ground:
                        lat_1 = self.ground_core_sites[core][0]
                        lon_1 = self.ground_core_sites[core][1]
                        lat_2 = self.ground_core_sites[core_2][0]
                        lon_2 = self.ground_core_sites[core_2][1]
                        dist = self.great_circle_distance_calculation(lat_1=lat_1, lat_2=lat_2, lon_1=lon_1,
                                                                      lon_2=lon_2)
                        distance_core_sites_array.append(dist)
                        disance_core_sites_dict[core, core_2] = dist
            return (np.sum(distances_array) + 0.5 * np.sum(distance_core_sites_array)) * cost_per_km/1000
        except:
            self.calculate_physical_distance_between_core_sites_and_stations()
            distances_array = []
            for core, ground in self.distance_dict.keys():
                if ground in self.available_sites[core]:
                    distances_array.append(self.distance_dict[core, ground])
            for core in self.ground_core_sites.keys():
                ground = self.available_sites[core][0]
                for core_2 in self.ground_core_sites.keys():
                    if core_2 != core and self.available_sites[core_2][0] == ground:
                        lat_1 = self.ground_core_sites[core][0]
                        lon_1 = self.ground_core_sites[core][1]
                        lat_2 = self.ground_core_sites[core_2][0]
                        lon_2 = self.ground_core_sites[core_2][1]
                        dist = self.great_circle_distance_calculation(lat_1=lat_1, lat_2=lat_2, lon_1=lon_1,
                                                                      lon_2=lon_2)
                        distance_core_sites_array.append(dist)
            return (np.sum(distances_array) + 0.5 * np.sum(distance_core_sites_array)) * cost_per_km / 1000

    def return_average_distance(self):
        try:
            distances_array = []
            for core, ground in self.distance_dict.keys():
                if ground in self.available_sites[core]:
                    distances_array.append(self.distance_dict[core,ground])
            return np.mean(distances_array), np.std(distances_array)
        except:
            self.calculate_physical_distance_between_core_sites_and_stations()
            distances_array = []
            for core, ground in self.distance_dict.keys():
                if ground in self.available_sites[core]:
                    distances_array.append(self.distance_dict[core,ground])
            return np.mean(distances_array), np.std(distances_array)



    def select_available_sites_based_on_distance(self, distance):
        # distance if the distance allowed
        available_sites = {}
        for core_site in self.ground_core_sites.keys():
            available_sites[core_site] = []
        if not hasattr(self, "distance_dict"):
            self.calculate_physical_distance_between_core_sites_and_stations()
        for core_site, ground_station in self.distance_dict.keys():
            ground_station_for_core_sites = self.ground_core_sites[core_site][2]
            if ground_station not in ground_station_for_core_sites:
                continue
            else:
                # ground station can be reached by node check if within distance
                if self.distance_dict[core_site, ground_station] <= distance * 1000:
                    available_sites[core_site].append(ground_station)
        self.available_sites = available_sites
        return available_sites


    def select_available_sites_based_on_given_sites(self, available_sites):
        """
        Sets the available_sites variable based on a given set of available sites.
        Checks there is a physical connection to all available sites before setting
        variable. Raises ValueError if the available sites does not have direct connection
        :param available_sites: dictionary of available ground stations for each core site {core_site: [ground_stations]}
        """
        for core_site in available_sites.keys():
            for ground_station in available_sites[core_site]:
                if ground_station not in self.ground_core_sites[core_site][2]:
                    print("Ground station provided is not physically connected to core site")
                    raise ValueError
        self.available_sites = available_sites

    def select_single_available_site(self):
        available_sites = {}
        min_distance = {}
        for core_site in self.ground_core_sites.keys():
            available_sites[core_site] = None
            min_distance[core_site] = np.infty
        if not hasattr(self, "distance_dict"):
            self.calculate_physical_distance_between_core_sites_and_stations()
        min_ground_station = None
        for core_site, ground_station in self.distance_dict.keys():
            ground_station_for_core_sites = self.ground_core_sites[core_site][2]
            if ground_station not in ground_station_for_core_sites:
                continue
            elif core_site == "Beijing":
                if self.distance_dict[core_site, ground_station] <= min_distance[core_site]:
                    available_sites[core_site]= [ground_station]
                    min_distance[core_site] = self.distance_dict[core_site, ground_station]
            else:
                # ground station can be reached by node check if within distance
                if self.distance_dict[core_site, ground_station] <= min_distance[core_site]:
                    available_sites[core_site]= [ground_station]
                    min_distance[core_site] = self.distance_dict[core_site, ground_station]
        self.available_sites = available_sites
        return available_sites


    def select_available_sites_based_on_connectivity(self):
        # selects sites based on if there is a connection on the physical network or not.
        # The limit of select_available_sites_based_on_distance if distance = infty
        available_sites = {}
        for core_site in self.ground_core_sites.keys():
            available_sites[core_site] = []
        if not hasattr(self, "distance_dict"):
            self.calculate_physical_distance_between_core_sites_and_stations()
        for core_site, ground_station in self.distance_dict.keys():
            ground_station_for_core_sites = self.ground_core_sites[core_site][2]
            if ground_station not in ground_station_for_core_sites:
                continue
            else:
                available_sites[core_site].append(ground_station)
        self.available_sites = available_sites
        return available_sites


    def get_commodity_set(self):
        commodity_set = []
        for ground_station_1 in self.ground_station_sites.keys():
            for ground_station_2 in self.ground_station_sites.keys():
                if ground_station_1 != ground_station_2:
                    if (ground_station_2, ground_station_1) not in commodity_set:
                        commodity_set.append(ground_station_1, ground_station_2)
        self.commodity_set = commodity_set
        return commodity_set


    def calculate_satellite_rate_based_on_expected_transmission_requirements(self,elevation_heights, file_path, satellite_sites, satellite_RAAN):
        #### get effective transmission_requirements
        number_core_sites_per_gs_site = {}
        for site in self.ground_station_sites.keys():
            number_sites = 0
            for core_site in self.ground_core_sites.keys():
                if site in self.available_sites[core_site]:
                    number_sites += 1
            number_core_sites_per_gs_site[site] = number_sites
        Tij_eff = {}
        for site in self.ground_station_sites.keys():
            Tij_eff_site = 0.0
            for core_site in self.ground_core_sites.keys():
                if site in self.available_sites[core_site]:
                    Tij_eff_site_core_site = 0.0
                    for core_site_2 in self.ground_core_sites.keys():
                        if core_site_2 != core_site:
                            if site not in self.available_sites[core_site_2]:
                                Tij_eff_site_core_site += self.Tij[core_site, core_site_2] + self.Tij[core_site_2, core_site]
                            else:
                                Tij_eff_site_core_site += (self.Tij[core_site, core_site_2] + self.Tij[core_site_2, core_site]) * (1-1/(self.available_sites[core_site_2]))
                    Tij_eff_site += Tij_eff_site_core_site / (2 * len(self.available_sites[core_site]))
            Tij_eff[site] = Tij_eff_site

        time_angle_dfs = {}
        times = []
        for site in satellite_sites:
            time_angle_df = satellite_rate_calculator.get_time_elevation_angle(file_path, site, satellite_RAAN)
            time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
            times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
            times = list(set(times) | set(times_current))
            time_angle_dfs[site] = time_angle_df_non_zero_capacity
        times = sorted(times)
        # create a timeset
        sites_at_times = {}
        time_list_prev = []
        prev_time = None
        for time in times:
            if prev_time == None:
                prev_time = time
                sites = set()
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        sites = sites.union(set(site))
            if time - prev_time > 10:
                for t in time_list_prev:
                    sites_at_times[t] = sites
                time_list_prev = [time]
                prev_time = time
                sites = set()
            else:
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        sites = sites.union(set(site))
                time_list_prev.append(time)
                prev_time = time
        for t in time_list_prev:
            sites_at_times[t] = sites
        key_rates_sites = {}
        for site in time_angle_dfs.keys():
            elevation_df = elevation_heights[
                min(list(elevation_heights.keys()), key=lambda x: abs(x * 1000 - satellite_rate_calculator.site_height[site]))]
            time_prev = None
            key_size_tot = 0.0
            for ind in time_angle_dfs[site].index:
                time = time_angle_dfs[site]["currentTime"][ind]

                if time_prev == None:
                    time_prev = time
                    continue
                time_step = time - time_prev
                if time_step > 1000:
                    time_prev = None
                    continue
                eff_rate_time = Tij_eff[site]
                eff_mod_rate_time = 0.0
                for site_2 in sites_at_times[time]:
                    eff_mod_rate_time += Tij_eff[site_2]
                eff_rate_time = eff_rate_time / eff_mod_rate_time
                elevation = time_angle_dfs[site]["elevationAngle"][ind]
                elevation_current = elevation_df.iloc[(elevation_df["Elevation (deg)"] - elevation).abs().argsort()[:1]]
                key_size_tot += time_step * elevation_current["Rate (b)"].values[0] * eff_rate_time
                time_prev = time
            key_rate = key_size_tot / (24 * 60 * 60)
            key_rates_sites[site] = key_rate
        return key_rates_sites


    def calculate_satellite_rate_based_on_expected_transmission_requirements_segragation(self, file_path, satellite_sites, satellite_RAAN):
        #### get effective transmission_requirements
        number_core_sites_per_gs_site = {}
        for site in self.ground_station_sites.keys():
            number_sites = 0
            for core_site in self.ground_core_sites.keys():
                if site in self.available_sites[core_site]:
                    number_sites += 1
            number_core_sites_per_gs_site[site] = number_sites
        Tij_eff = {}
        for site in self.ground_station_sites.keys():
            Tij_eff_site = 0.0
            # Tij eff = sum_site1 (sum_site2 Tij[(site1,site2)] * p(site1 & site2 not in the same GS| site 1 is in current GS) * p(site1 is in current GS))
            # if site 1 not in current GS list then p(site1 is in current GS) = 0 and can be excluded from sum
            # if site 1 is in current GS list, assuming equal prob of partitioning p(site1 is in current GS) = 1/|available sites for GS|
            # if site 1 is in current GS, p(site2 also not in GS) is given by : 1 if site2 not in set of available
            # sites for GS or 1-1/(number of GS per core site[site2]) if it is
            for core_site in self.ground_core_sites.keys():
                if site in self.available_sites[core_site]:
                    Tij_eff_site_core_site = 0.0
                    for core_site_2 in self.ground_core_sites.keys():
                        if core_site_2 != core_site:
                            if site not in self.available_sites[core_site_2]:
                                Tij_eff_site_core_site += self.Tij[core_site, core_site_2] + self.Tij[core_site_2, core_site]
                            else:
                                Tij_eff_site_core_site += (self.Tij[core_site, core_site_2] + self.Tij[core_site_2, core_site]) * (1-1/(len(self.available_sites[core_site_2])))
                    Tij_eff_site += Tij_eff_site_core_site / (len(self.available_sites[core_site]))
            Tij_eff[site] = Tij_eff_site

        time_angle_dfs = {}
        times = []
        for site in satellite_sites:
            time_angle_df = satellite_rate_calculator.get_time_elevation_angle(file_path, site, satellite_RAAN)
            time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
            times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
            times = list(set(times) | set(times_current))
            time_angle_dfs[site] = time_angle_df_non_zero_capacity
        times = sorted(times)
        # create a timeset
        sites_at_times = {}
        time_list_prev = []
        prev_time = None
        for time in times:
            if prev_time == None:
                prev_time = time
                sites = []
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        if site not in sites:
                            sites.append(site)
            if time - prev_time > 10:
                for t in time_list_prev:
                    sites_at_times[t] = sites
                time_list_prev = [time]
                prev_time = time
                sites = []
            else:
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        if site not in sites:
                            sites.append(site)
                time_list_prev.append(time)
                prev_time = time
        for t in time_list_prev:
            sites_at_times[t] = sites
        segragations_sites = {}
        for site in time_angle_dfs.keys():
            time_prev = None
            segragation = []
            time_zone_dark_low = 1800 + time_zone[site] * 3600
            time_zone_dark_high = 45000 + time_zone[site]* 3600
            for ind in time_angle_dfs[site].index:
                time = time_angle_dfs[site]["currentTime"][ind]

                if time_prev == None:
                    time_prev = time
                    continue
                time_step = time - time_prev
                if time_step > 2:
                    eff_rate_time = Tij_eff[site]
                    eff_mod_rate_time = 0.0
                    for site_2 in sites_at_times[time_prev]:
                        eff_mod_rate_time += Tij_eff[site_2]
                    if eff_mod_rate_time > 0.01:
                        eff_rate_time = eff_rate_time / eff_mod_rate_time
                    else:
                        eff_rate_time = 0.05
                    if time_zone_dark_low < 0:
                        if 0 < time_prev < time_zone_dark_high or 86400 + time_zone_dark_low < time_prev < 86400:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))
                    elif time_zone_dark_high > 86400:
                        if 0 < time_prev < time_zone_dark_high - 86400 or time_zone_dark_low < time_prev < 86400:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))
                    else:
                        if time_prev < time_zone_dark_high and time_prev > time_zone_dark_low:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))

                    time_prev = None
                    continue
                time_prev = time
            if time_prev != None:
                eff_rate_time = Tij_eff[site]
                eff_mod_rate_time = 0.0
                for site_2 in sites_at_times[time_prev]:
                    eff_mod_rate_time += Tij_eff[site_2]
                if eff_mod_rate_time > 0.01:
                    eff_rate_time = eff_rate_time / eff_mod_rate_time
                else:
                    eff_rate_time = 0.01
                if time_zone_dark_low < 0:
                    if 0 < time_prev < time_zone_dark_high or 86400 + time_zone_dark_low < time_prev < 86400:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))
                elif time_zone_dark_high > 86400:
                    if 0 < time_prev < time_zone_dark_high - 86400 or time_zone_dark_low < time_prev < 86400:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))
                else:
                    if time_prev < time_zone_dark_high and time_prev > time_zone_dark_low:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))

            segragations_sites[site] = segragation
        return segragations_sites




    def calculate_satellite_rate_based_on_key_rate_segragation(self, file_path, satellite_sites, satellite_RAAN,
                                                               count, ci, ni, x, x0i, xb, theta_max, dt_range,
                                                               main_params, opt_params,
                                                               args_fixed, bounds, cons, options, header, opt_head,
                                                               fulldata, optdata,
                                                               multidata, sys_params, sysLoss, satellite_trajectory_file, f_atm):
        #### get effective transmission_requirements
        number_core_sites_per_gs_site = {}
        for site in self.ground_station_sites.keys():
            number_sites = 0
            for core_site in self.ground_core_sites.keys():
                if site in self.available_sites[core_site]:
                    number_sites += 1
            number_core_sites_per_gs_site[site] = number_sites
        time_angle_dfs = {}
        times = []
        for site in satellite_sites:
            time_angle_df = satellite_rate_calculator.get_time_elevation_angle(file_path, site, satellite_RAAN)
            time_angle_df_non_zero_capacity = time_angle_df.loc[time_angle_df["elevationAngle"] > 0.2]
            times_current = time_angle_df_non_zero_capacity["currentTime"].unique()
            times = list(set(times) | set(times_current))
            time_angle_dfs[site] = time_angle_df_non_zero_capacity
        times = sorted(times)
        # Evaluate Key Rates for each pass for each satellite station pair
        FSeff_sat = {}
        for site in self.ground_station_sites.keys():
            loss_data = get_losses_input_file(satellite_trajectory_file + f"_{site}_{satellite_RAAN}.csv", main_params['loss'], f_atm,
                                              main_params['out']['tPrint'],
                                              main_params['out']['out_path'])
            FSeff = loss_data[:, 2]
            FSeff_sat[site] = FSeff
        K_ij_P = {}
        prev_time = None
        prev_ind = None
        current_pass_times = None
        for site in self.ground_station_sites.keys():
            prev_time = None
            prev_ind = None
            time_zone_dark_low = 1800 + time_zone[site] * 3600
            time_zone_dark_high = 45000 + time_zone[site]* 3600
            for ind in time_angle_dfs[site].index:
                time = time_angle_dfs[site]["currentTime"][ind]
                if prev_time == None:
                    prev_time = time
                    current_pass_times = time
                if time - prev_time > 100:
                    if prev_ind == None:
                        if time_zone_dark_low < 0:
                            if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                                segragation= [(ind, 1, "dark")]
                            else:
                                segragation= [(ind, 1, "light")]
                        elif time_zone_dark_high > 86400:
                            if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                                segragation= [(ind, 1, "dark")]
                            else:
                                segragation= [(ind, 1, "light")]
                        else:
                            if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                                segragation= [(ind, 1, "dark")]
                            else:
                                segragation= [(ind, 1, "light")]

                    else:
                        if time_zone_dark_low < 0:
                            if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "dark")]
                            else:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "light")]
                        elif time_zone_dark_high > 86400:
                            if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "dark")]
                            else:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "light")]
                        else:
                            if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "dark")]
                            else:
                                segragation= [(prev_ind, 0, "dark"),(ind, 1, "light")]
                    prev_ind = ind
                    current_pass_times = (current_pass_times, prev_time)
                    args_fixed[3] = FSeff_sat[site]
                    i = 0
                    max_rate = 0.0
                    while i <2 and max_rate < 0.001:
                        x = [np.random.uniform(low= xb[0][0], high = xb[0][1]),np.random.uniform(low= xb[1][0], high = xb[1][1]),
                             np.random.uniform(low= xb[2][0], high = xb[2][1]), np.random.uniform(low= xb[3][0], high = xb[3][1]),
                             np.random.uniform(low= xb[4][0], high = xb[4][1])]
                        main_params_temp = copy.deepcopy(main_params)
                        # main_params_temp['opt']['tOptimise'] = False


                        ##### args_fixed need to be obtained for each satellite individually....
                        ##### Currently this just uses the first satellite in the list....
                        multidata, count, max_rate = SKL_sys_loop(count, ci, ni, x, x0i, xb, theta_max, dt_range, main_params_temp, opt_params,
                                     args_fixed, bounds, cons, options, header, opt_head, fulldata, optdata,
                                     multidata, sys_params, sysLoss, segragation)
                        i += 1
                    K_ij_P[site,current_pass_times] = max_rate
                    prev_time = None
                else:
                    prev_time = time
            if prev_time != None:
                if prev_ind == None:
                    if time_zone_dark_low < 0:
                        if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                            segragation = [(ind, 1, "dark")]
                        else:
                            segragation = [(ind, 1, "light")]
                    elif time_zone_dark_high > 86400:
                        if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                            segragation = [(ind, 1, "dark")]
                        else:
                            segragation = [(ind, 1, "light")]
                    else:
                        if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                            segragation = [(ind, 1, "dark")]
                        else:
                            segragation = [(ind, 1, "light")]
                else:
                    if time_zone_dark_low < 0:
                        if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "dark")]
                        else:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "light")]
                    elif time_zone_dark_high > 86400:
                        if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "dark")]
                        else:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "light")]
                    else:
                        if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "dark")]
                        else:
                            segragation = [(prev_ind, 0, "dark"), (ind, 1, "light")]
                prev_ind = ind
                current_pass_times = (current_pass_times, prev_time)
                args_fixed[3] = FSeff_sat[site]
                i = 0
                max_rate = 0.0
                while i < 2 and max_rate < 0.001:
                    x = [np.random.uniform(low=xb[0][0], high=xb[0][1]), np.random.uniform(low=xb[1][0], high=xb[1][1]),
                         np.random.uniform(low=xb[2][0], high=xb[2][1]), np.random.uniform(low=xb[3][0], high=xb[3][1]),
                         np.random.uniform(low=xb[4][0], high=xb[4][1])]
                    main_params_temp = copy.deepcopy(main_params)
                    # main_params_temp['opt']['tOptimise'] = False
                    multidata, count, max_rate = SKL_sys_loop(count, ci, ni, x, x0i, xb, theta_max, dt_range,
                                                              main_params_temp, opt_params,
                                                              args_fixed, bounds, cons, options, header, opt_head,
                                                              fulldata, optdata,
                                                              multidata, sys_params, sysLoss, segragation)
                    i += 1
                K_ij_P[site, current_pass_times] = max_rate
                prev_time = None
        # create a timeset
        sites_at_times = {}
        time_list_prev = []
        prev_time = None
        current_pass_times = []
        for time in times:
            if prev_time == None:
                prev_time = time
                sites = []
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        if site not in sites:
                            sites.append(site)
            if time - prev_time > 2:
                for t in time_list_prev:
                    sites_at_times[t] = sites
                time_list_prev = [time]
                prev_time = time
                sites = []
            else:
                for site in time_angle_dfs.keys():
                    if time_angle_dfs[site]["currentTime"].isin([time]).any():
                        if site not in sites:
                            sites.append(site)
                time_list_prev.append(time)
                prev_time = time
        for t in time_list_prev:
            sites_at_times[t] = sites
        segragations_sites = {}
        for site in time_angle_dfs.keys():
            time_prev = None
            segragation = []
            current_pass_times = []
            time_zone_dark_low = 1800 + time_zone[site] * 3600
            time_zone_dark_high = 45000 + time_zone[site]* 3600
            for ind in time_angle_dfs[site].index:
                time = time_angle_dfs[site]["currentTime"][ind]

                if time_prev == None:
                    time_prev = time
                    current_pass_times.append(time)
                    continue
                time_step = time - time_prev
                if time_step > 2:
                    eff_mod_rate_time = 0.0
                    for site_current, time_scale_current in K_ij_P:
                        if site_current == site:
                            if any(time_scale_current[0] <= t <= time_scale_current[1] for t in current_pass_times):
                                eff_rate_time = K_ij_P[site_current, time_scale_current]
                                eff_mod_rate_time += K_ij_P[site_current, time_scale_current]
                        elif site_current in sites_at_times[time_prev]:
                            if any(time_scale_current[0] <= t <= time_scale_current[1] for t in current_pass_times):
                                eff_mod_rate_time += K_ij_P[site_current, time_scale_current]
                    if eff_mod_rate_time > 0.01:
                        eff_rate_time = eff_rate_time / eff_mod_rate_time
                    else:
                        eff_rate_time = 0.05
                    if time_zone_dark_low < 0:
                        if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))
                    elif time_zone_dark_high > 86400:
                        if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))
                    else:
                        if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                            segragation.append((ind, eff_rate_time, "dark"))
                        else:
                            segragation.append((ind, eff_rate_time, "light"))


                    current_pass_times = []
                    time_prev = None
                    continue
                else:
                    current_pass_times.append(time)
                time_prev = time
            if time_prev != None:
                eff_mod_rate_time = 0.0
                for site_current, time_scale_current in K_ij_P:
                    if site_current == site:
                        if any(time_scale_current[0] <= t <= time_scale_current[1] for t in current_pass_times):
                            eff_rate_time = K_ij_P[site_current, time_scale_current]
                            eff_mod_rate_time += K_ij_P[site_current, time_scale_current]
                    elif site_current in sites_at_times[time_prev]:
                        if any(time_scale_current[0] <= t <= time_scale_current[1] for t in current_pass_times):
                            eff_mod_rate_time += K_ij_P[site_current, time_scale_current]
                if eff_mod_rate_time > 0.01:
                    eff_rate_time = eff_rate_time / eff_mod_rate_time
                else:
                    eff_rate_time = 0.05
                if time_zone_dark_low < 0:
                    if 0 < prev_time < time_zone_dark_high or 86400 + time_zone_dark_low < prev_time < 86400:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))
                elif time_zone_dark_high > 86400:
                    if 0 < prev_time < time_zone_dark_high - 86400 or time_zone_dark_low < prev_time < 86400:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))
                else:
                    if prev_time < time_zone_dark_high and prev_time > time_zone_dark_low:
                        segragation.append((ind, eff_rate_time, "dark"))
                    else:
                        segragation.append((ind, eff_rate_time, "light"))


            segragations_sites[site] = segragation
        return segragations_sites