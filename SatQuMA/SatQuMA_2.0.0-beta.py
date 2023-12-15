# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:20:58 2021

@author: Duncan McArthur

Modified by:

@author: Vasileios Karavias
"""

from os.path import join

import numpy as np
import pandas as pd
import csv
import os
from input.read_input import get_input
from input.write_input import (print_header, print_input)
from channel.atmosphere.atmos_data import get_f_atm
from optimize.optimiser import (set_constraints, opt_arrays, out_heads)
from key.get_key import (SKL_main_loop, SKL_main_loop_2)
from output.outputs import (out_arrays, get_timings, format_time)
from generate_graphs_for_satellite_problem import Network_Structure
from SatQuMA.satellite_utils import (get_Tij, get_satellite_ground_stations, get_ground_core_stations)
from store_satellite_data_simple_graph import Store_Network_Problem

###############################################################################
site_height = {"London": 11, "Madrid": 66, "Edinburgh": 47, "Athens": 20, "Barcelona": 12, "Paris": 35, "Nantes": 20,
               "Lyon": 237, "Lisbon": 2, "Bern": 540, "Freiburg": 280, "Milan": 120, "Florence": 50, "Naples": 100,
               "Frankfurt": 110, "Dortmund":86 , "Berlin": 35, "Ioannina": 480}

site_height = {"London": 11, "Frankfurt": 112, "Graz": 490, "Johanessburg": 1753, "Sao Paulo": 760, "Tokyo": 40,
                    "Auckland": 200,
                    "New Dehli": 216, "Mumbai": 14, "Bangalore": 920,
                    "Xinglong": 890, "Nanshan": 2028, "Perth": 0, "Brisbane": 32, "Melborne": 31, "Albany": 45,
                    "Denver": 1609}
###############################################################################

# Set the inpufile path and filename(s);
# Path to input files, empty string means location of execution
inputpath    = 'inputfiles'

# File name for main input parameters
filename     = 'input.txt'

# File name for advanced/optional input parameters (optional)
filename_adv = 'input-adv.txt'

###############################################################################

# Print the SatQuMA header
# print_header()
#
# # Start clock and CPU timer
# tc00, tp00 = get_timings()
#
# # Retrieve input parameters
# protocol, main_params, adv_params = get_input(join(inputpath,filename),
#                                               join(inputpath,filename_adv))
#
# # Print input parameters
# if main_params['out']['tPrint']:
#     print('-'*80)
#     print('Security protocol:',protocol)
#     print('-'*80)
#     print_input(main_params,'Main')
#     print('-'*80)
#     print_input(adv_params,'Adv')
#     print('-'*80)
#
# ###############################################################################
#
# # Get atmospheric function (if required), else a dummy function is returned
# f_atm = get_f_atm(main_params['loss'])
#
# # Get bounds/constraints for the optimiser
# bounds, cons, options, x, xb = set_constraints(main_params['opt'],
#                                                main_params['fixed'],
#                                                adv_params['opt'],protocol)
#
# # Initialise counting/arrays
# ni, ci , x0i     = opt_arrays(main_params['iter'])
#
# # Generate output file headers based on input parameters
# header, opt_head = out_heads(protocol,main_params['opt']['tOptimise'],
#                              adv_params['opt']['method'])
# # Initialise output data arrays - arrays are sized based on the length of the
# # headers generated above
# fulldata, optdata, multidata = out_arrays(ni,header,opt_head,
#                                           main_params['opt']['tOptimise'],
#                                           main_params['out']['tdtOptData'])
#
#
# # ground_core_sites = {"London": [np.radians(51.5072),np.radians(-0.1276),["London"]], "Madrid" : [np.radians(40.4168), np.radians(-3.7038), ["Madrid"]],
# #                      "Barcelona": [np.radians(41.3874), np.radians(2.1686), ["Madrid"]]}
# #
# #
# # ground_station_sites = {"London": [np.radians(51.5072),np.radians(-0.1276)], "Madrid": [np.radians(40.4168), np.radians(-3.7038)]}
# #
# # Tij  = {("London", "Madrid"): 1350, ("London", "Barcelona"): 720, ("Madrid", "Barcelona"): 240, ("Madrid", "London"): 1350, ("Barcelona","London"): 720, ("Barcelona","Madrid"): 240}
#
#
# ground_core_sites = get_ground_core_stations()
#
# ground_station_sites = get_satellite_ground_stations()
#
# Tij = get_Tij()
#
# network = Network_Structure(ground_core_sites, ground_station_sites, Tij)
#
#
#
# # Get key!!!
# SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
#                   options, header, opt_head, fulldata, optdata, multidata, satellite_trajectory_file = "new_satelliteElevAngle_London_110.csv", allocation_method = "ClashFullBypass", file_path="new_satelliteElevAngle",
#                                                      satellite_sites=list(site_height.keys()), satellite_RAAN=110)
# SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
#                   options, header, opt_head, fulldata, optdata, multidata, satellite_trajectory_file = "new_satelliteElevAngle_London_110.csv", allocation_method = "ClashFullBypassTransmissionModulated", file_path="new_satelliteElevAngle",
#                                                      satellite_sites=list(site_height.keys()), satellite_RAAN=110, network_structure=network)
# "corrected_satellite_rates_random_allocation_n_2"
# SKL_main_loop(main_params,adv_params,x,x0i,xb,ci,ni,f_atm,bounds,cons,
#               options,header,opt_head,fulldata,optdata,multidata)
six_ground_stations = ["London", "Madrid", "Paris", "Berlin", "Athens", "Naples"]
seventeen_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
               "Lyon", "Lisbon", "Bern", "Freiburg", "Milan", "Florence", "Naples",
                "Dortmund" , "Berlin", "Ioannina"]
sixteen_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
               "Lyon", "Lisbon", "Bern", "Freiburg", "Florence", "Naples",
                "Dortmund" , "Berlin", "Ioannina"]
fifteen_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
                "Lisbon", "Bern", "Freiburg", "Florence", "Naples",
                "Dortmund" , "Berlin", "Ioannina"]
fourteen_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
                "Bern", "Freiburg", "Florence", "Naples",
                "Dortmund" , "Berlin", "Ioannina"]
thirteen_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
                "Bern", "Freiburg", "Florence", "Naples",
                 "Berlin", "Ioannina"]
twelve_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
                "Bern", "Florence", "Naples",
                 "Berlin", "Ioannina"]
eleven_ground_stations = ["London", "Madrid", "Edinburgh", "Athens", "Barcelona", "Paris", "Nantes",
                "Bern", "Florence", "Naples","Berlin"]
ten_ground_stations = ["London", "Madrid",  "Athens", "Barcelona", "Paris", "Nantes",
                "Bern", "Florence", "Naples","Berlin"]
nine_ground_stations = ["London", "Madrid",  "Athens","Paris", "Nantes",
                "Bern", "Florence", "Naples","Berlin"]

eight_ground_stations = ["London", "Madrid",  "Athens","Paris", "Nantes",
                "Bern", "Naples","Berlin"]
seven_ground_stations = ["London", "Madrid",  "Athens","Paris", "Nantes",
                "Bern", "Naples"]
five_ground_stations = ["London", "Madrid",  "Athens","Paris", "Naples"]
four_ground_stations = ["London", "Madrid",  "Athens","Paris"]
# 80,90,100, 110, 120,
three_ground_stations = ["London", "Madrid", "Athens"]
core_station_list =["London", "Birmingham", "Manchester", "Paris", "Marseille", "Frankfurt", "Dusseldorf", "Amsterdam", "Madrid", "Zurich", "Milan", "Stockholm", "Copenhagen", "Sofia", "Bucharest"]


def get_resulting_keys(file_input, file_output, ground_station_list):
    transmission_keys = pd.read_csv(file_input)
    possible_ids = transmission_keys["ID"].unique()
    key_rates = []
    for id in possible_ids:
        transmission_keys = transmission_keys[transmission_keys["ID"] == id].drop(["ID"], axis=1)
        for index, row in transmission_keys.iterrows():
            if row["ground_core_station"] in ground_station_list:
                key_rates.append({"ID": id, "ground_core_station": row["ground_core_station"],
                                  "satellite_path": row["satellite_path"], "capacity": row["capacity"]})
    dictionary_fieldnames = ["ID", "ground_core_station", "satellite_path", "capacity"]
    if os.path.isfile(file_output + '.csv'):
        with open(file_output + '.csv', mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
            writer.writerows(key_rates)
    else:
        with open(file_output + '.csv', mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
            writer.writeheader()
            writer.writerows(key_rates)


###############################################################################
# UNCOMMENT THIS FOR SATELLITE RATE CALCULATIONS FOR SIMPLE CASE
# Method of Segragation: choice from ["None", "ClashFixedTime", "ClashFullBypass", "ClashFullBypassTransmissionModulated", "ClashFullBypassKeyModulated"]
# Method of network segragation: choice from ["SingleConnection","DistanceBased", "ConnectionBased" ]
# net = Store_Network_Problem(allocation_method = "ClashFullBypassTransmissionModulated", RAAN_list = [80, 90, 100, 110, 120, 130], core_station_list =core_station_list, ground_station_list = sixteen_ground_stations, network_segragation = "SingleConnection", distance = None, n = 1)
# net.store_rate_satellite_ground_station_file(file_path = "Data_Centres_Satellite_Rate_Transmission_Modulated_fewer_core_sites_16", graph_id = 0)
# net.store_simple_Tij_file(file_path = "Data_Centres_Tij_Transmission_Modulated_fewer_core_sites_16", graph_id = 0)
# # mean, std = net.network.return_average_distance()
# # print("Mean Distance: " + str(mean) + ". STD Distance: " + str(std))
# # Stop clock and CPU timers
# get_resulting_keys(file_input = "corrected_satellites_rates_transmission_modulated_fewer_core_sites_18.csv", file_output ="corrected_satellites_rates_transmission_modulated_fewer_core_sites_3", ground_station_list= three_ground_stations)
# get_resulting_keys(file_input = "corrected_satellites_rates_transmission_modulated_fewer_core_sites_18.csv", file_output ="corrected_satellites_rates_transmission_modulated_fewer_core_sites_15", ground_station_list= fifteen_ground_stations)




# UNCOMMENT THIS FOR SATELLITE RATE CALCULATIONS FOR FULL OPTIMISATION
# net = Store_Network_Problem(allocation_method = "ClashFullBypassTransmissionModulated", RAAN_list = [80, 90, 100, 110, 120,130], core_station_list =None, ground_station_list = None, network_segragation = "ConnectionBased", distance = 650, n = 1)
# # net.store_rate_satellite_ground_station_file(file_path = "corrected_satellites_rates_transmission_modulated_full_range_fewer_core_sites_3", graph_id = 0)
# net.store_Tij_file(file_path = "Data_Centres_Full_Problem_Tij_satellites_path_transmission_modulated_full_range_18", graph_id = 0)
# net.store_key_dict_file(file_path = "Data_Centres_key_dict_transmission_modulated_full_range_18", graph_id = 0)
# net.store_core_site_keys(file_path = "Data_Centres_core_sites_transmission_modulated_full_range_18", graph_id = 0)

# tc11, tp11 = get_timings()
#
# tc   = tc11-tc00      # Calculation duration from clock
# tp   = tp11-tp00      # Calculation duration from CPU
# print('\nFinal clock timer:',format_time(tc))
# print('Final CPU timer:  ',format_time(tp),'\n')
# print('All done!')






#################################### DISTANCE CALCULATION


# ground_core_sites = {"London": [np.radians(51.5072),np.radians(-0.1276),["London"]], "Madrid" : [np.radians(40.4168), np.radians(-3.7038), ["Madrid"]],
#                      "Barcelona": [np.radians(41.3874), np.radians(2.1686), ["Madrid"]]}
#
#
# ground_station_sites = {"London": [np.radians(51.5072),np.radians(-0.1276)], "Madrid": [np.radians(40.4168), np.radians(-3.7038)]}
#
# Tij  = {("London", "Madrid"): 1350, ("London", "Barcelona"): 720, ("Madrid", "Barcelona"): 240, ("Madrid", "London"): 1350, ("Barcelona","London"): 720, ("Barcelona","Madrid"): 240}

# core_station_list = None
# ground_station_list = three_ground_stations
# #
# ground_core_sites = get_ground_core_stations(core_station_list =core_station_list, ground_station_list = ground_station_list)
#
# ground_station_sites = get_satellite_ground_stations(ground_station_list = ground_station_list)
#
# Tij = get_Tij(core_station_list = core_station_list)
#
# network = Network_Structure(ground_core_sites, ground_station_sites, Tij)
#
#
# network.select_single_available_site()
# cost = network.calculate_cost_core_network_graph(cost_per_km=0.02)
# print(str(cost))
# cost = network.calculate_cost_core_network_graph(cost_per_km=0.001)
# print(str(cost))
# cost = network.calculate_cost_core_network_graph(cost_per_km=0.0001)
# print(str(cost))



#### Global Network Problem ######

nine_gs = ["Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "Mumbai","Xinglong", "Melborne", "Denver"]
ten_gs =  ["Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "Mumbai","Xinglong", "Melborne", "Denver"]
eleven_gs  = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "Mumbai","Xinglong", "Melborne", "Denver"]
twelve_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Xinglong", "Melborne", "Denver"]
thirteen_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Bangalore", "Xinglong", "Melborne", "Denver"]
fourteen_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Bangalore", "Xinglong", "Melborne", "Perth","Denver"]
fifteen_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Bangalore", "Xinglong", "Melborne", "Perth", "Brisbane","Denver"]
sixteen_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Bangalore", "Xinglong", "Melborne", "Perth","Brisbane", "Albany","Denver"]
seventeen_gs = ["London","Frankfurt", "Graz", "Johanessburg", "Sao Paulo", "Tokyo", "Auckland", "New Delhi", "Mumbai","Bangalore", "Xinglong",  "Melborne", "Perth","Brisbane", "Albany","Nashville","Denver"]
# net = Store_Network_Problem(allocation_method = "ClashFullBypassTransmissionModulated", RAAN_list = [80,90,100,110,120,130], core_station_list =None, ground_station_list = None, network_segragation = "SingleConnection", distance = None, n = 2)
# net.store_rate_satellite_ground_station_file(file_path = "Data_Centres_Global_30_Year_Future_Satellite_Rate_Transmission_Modulated", graph_id = 0)
# net.store_simple_Tij_file(file_path = "Data_Centres_Global_30_Year_Future_Tij_Transmission_Modulated_18", graph_id = 0)


net = Store_Network_Problem(allocation_method = "ClashFullBypassTransmissionModulated", RAAN_list = [80, 90, 100, 110, 120,130], core_station_list =None, ground_station_list = nine_gs, network_segragation = "ConnectionBased", distance = None, n = 1)
# # net.store_rate_satellite_ground_station_file(file_path = "corrected_satellites_rates_transmission_modulated_full_range_fewer_core_sites_3", graph_id = 0)
net.store_Tij_file(file_path = "Data_Centres_Global_10_Years_Full_Problem_Tij_satellites_path_transmission_modulated_full_range", graph_id = 0)
net.store_key_dict_file(file_path = "Data_Centres_Global_10_Years_key_dict_transmission_modulated_full_range", graph_id = 0)
net.store_core_site_keys(file_path = "Data_Centres_Global_10_Years_core_sites_transmission_modulated_full_range", graph_id = 0)


core_station_list = None
ground_station_list = seventeen_gs
#
ground_core_sites = get_ground_core_stations(core_station_list =core_station_list, ground_station_list = ground_station_list)

ground_station_sites = get_satellite_ground_stations(ground_station_list = ground_station_list)

Tij = get_Tij(core_station_list = core_station_list)

network = Network_Structure(ground_core_sites, ground_station_sites, Tij)


network.select_single_available_site()
cost = network.calculate_cost_core_network_graph(cost_per_km=0.02)
print(str(cost))
# cost = network.calculate_cost_core_network_graph(cost_per_km=0.001)
# print(str(cost))
cost = network.calculate_cost_core_network_graph(cost_per_km=0.0001)
print(str(cost))
