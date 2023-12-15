from satellites_graphs.generate_graphs_for_satellite_problem import Network_Structure
from satellites_graphs.SatQuMA.satellite_utils import (get_ground_core_stations, get_satellite_ground_stations,get_Tij)
from satellites_graphs.SatQuMA.input.read_input import get_input
from satellites_graphs.SatQuMA.input.write_input import (print_header, print_input)
from satellites_graphs.SatQuMA.channel.atmosphere.atmos_data import get_f_atm
from satellites_graphs.SatQuMA.optimize.optimiser import (set_constraints, opt_arrays, out_heads)
from satellites_graphs.SatQuMA.key.get_key import (SKL_main_loop, SKL_main_loop_2)
from satellites_graphs.SatQuMA.output.outputs import (out_arrays, get_timings, format_time)
import numpy as np
import csv
import os
from os.path import join


class Store_Network_Problem:


    def __init__(self, allocation_method, RAAN_list, core_station_list =None, ground_station_list = None, network_segragation = "ConnectionBased", distance = None, n = 1, is_local = False, file_path = "global_satelliteElevAngle"):
        if is_local:
            self.site_height = {"London": 11, "Madrid": 66, "Edinburgh": 47, "Athens": 20, "Barcelona": 12, "Paris": 35,
                            "Nantes": 20,
                            "Lyon": 237, "Lisbon": 2, "Bern": 540, "Freiburg": 280, "Milan": 120, "Florence": 50,
                            "Naples": 100,
                            "Frankfurt": 110, "Dortmund": 86, "Berlin": 35, "Ioannina": 480}
        else:
            self.site_height = {"London": 11, "Frankfurt": 490, "Graz": 112, "Johanessburg": 1753, "Sao Paulo": 760, "Tokyo": 40,
                       "Auckland": 200,
                       "New Delhi": 216, "Mumbai": 14, "Bangalore": 920,
                       "Xinglong": 890, "Nanshan": 2028, "Perth": 0, "Brisbane": 32, "Melborne": 31, "Albany": 45,
                       "Denver": 1609, "Nashville": 169}
        ###############################################################################

        # Set the inpufile path and filename(s);
        # Path to input files, empty string means location of execution
        inputpath = "inputfiles"

 

        # File name for main input parameters
        filename = 'input.txt'

        # File name for advanced/optional input parameters (optional)
        filename_adv = 'input-adv.txt'

        ###############################################################################

        # Print the SatQuMA header
        print_header()

        # Start clock and CPU timer
        tc00, tp00 = get_timings()

        # Retrieve input parameters
        protocol, main_params, adv_params = get_input(join(inputpath, filename),
                                                      join(inputpath, filename_adv))

        # Print input parameters
        if main_params['out']['tPrint']:
            print('-' * 80)
            print('Security protocol:', protocol)
            print('-' * 80)
            print_input(main_params, 'Main')
            print('-' * 80)
            print_input(adv_params, 'Adv')
            print('-' * 80)

        ###############################################################################

        # Get atmospheric function (if required), else a dummy function is returned
        f_atm = get_f_atm(main_params['loss'])

        # Get bounds/constraints for the optimiser
        bounds, cons, options, x, xb = set_constraints(main_params['opt'],
                                                       main_params['fixed'],
                                                       adv_params['opt'], protocol)

        # Initialise counting/arrays
        ni, ci, x0i = opt_arrays(main_params['iter'])

        # Generate output file headers based on input parameters
        header, opt_head = out_heads(protocol, main_params['opt']['tOptimise'],
                                     adv_params['opt']['method'])
        # Initialise output data arrays - arrays are sized based on the length of the
        # headers generated above
        fulldata, optdata, multidata = out_arrays(ni, header, opt_head,
                                                  main_params['opt']['tOptimise'],
                                                  main_params['out']['tdtOptData'])

        # ground_core_sites = {"London": [np.radians(51.5072),np.radians(-0.1276),["London"]], "Madrid" : [np.radians(40.4168), np.radians(-3.7038), ["Madrid"]],
        #                      "Barcelona": [np.radians(41.3874), np.radians(2.1686), ["Madrid"]]}
        #
        #
        # ground_station_sites = {"London": [np.radians(51.5072),np.radians(-0.1276)], "Madrid": [np.radians(40.4168), np.radians(-3.7038)]}
        #
        # Tij  = {("London", "Madrid"): 1350, ("London", "Barcelona"): 720, ("Madrid", "Barcelona"): 240, ("Madrid", "London"): 1350, ("Barcelona","London"): 720, ("Barcelona","Madrid"): 240}

        ground_core_sites = get_ground_core_stations(core_station_list =core_station_list, ground_station_list = ground_station_list)

        ground_station_sites = get_satellite_ground_stations(ground_station_list = ground_station_list)

        Tij = get_Tij(core_station_list = core_station_list)

        network = Network_Structure(ground_core_sites, ground_station_sites, Tij)
        max_rates = {}
        for RAAN in  RAAN_list:
            max_rates[RAAN] = {}
            segragations = None
            for ground_station in ground_station_sites.keys():
        # # Get key!!!
                max_rate = 0.0
                i = 0
                while i < 3 and max_rate < 0.0001:
                    x = [np.random.uniform(low=xb[0][0], high=xb[0][1]), np.random.uniform(low=xb[1][0], high=xb[1][1]),
                         np.random.uniform(low=xb[2][0], high=xb[2][1]), np.random.uniform(low=xb[3][0], high=xb[3][1]),
                         np.random.uniform(low=xb[4][0], high=xb[4][1])]
                    print(f"Current RAAN {RAAN}, Current ground station {ground_station}")
                    if segragations == None:
                        if ground_station_list != None:
        #
                            max_rate = SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
                                            options, header, opt_head, fulldata, optdata, multidata, ground_station = ground_station,
                                            satellite_trajectory_file=file_path + f"_{ground_station}_{RAAN}.csv",
                                            allocation_method=allocation_method, file_path=file_path,
                                            satellite_sites=ground_station_list, satellite_RAAN=RAAN, network_structure=network, method_for_network_segragation = network_segragation, distance = distance, n = n)
                            if allocation_method == "ClashFullBypassKeyModulated" or allocation_method == "ClashFullBypassTransmissionModulated":
                                max_rate, segragations = max_rate
        #
                        else:
                            max_rate = SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
                                            options, header, opt_head, fulldata, optdata, multidata, ground_station = ground_station,
                                             satellite_trajectory_file=file_path + f"_{ground_station}_{RAAN}.csv",
                                            allocation_method=allocation_method, file_path=file_path,
                                            satellite_sites=list(self.site_height.keys()), satellite_RAAN=RAAN, network_structure=network, method_for_network_segragation = network_segragation, distance = distance, n = n)
                            if allocation_method == "ClashFullBypassKeyModulated" or allocation_method == "ClashFullBypassTransmissionModulated":
                                max_rate, segragations = max_rate
                    else:
                        if ground_station_list != None:
        #
                            max_rate = SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
                                                        options, header, opt_head, fulldata, optdata, multidata, ground_station = ground_station,
                                                       satellite_trajectory_file=file_path + f"_{ground_station}_{RAAN}.csv",
                                                       allocation_method=allocation_method, file_path=file_path,
                                                       satellite_sites=ground_station_list, satellite_RAAN=RAAN,
                                                       network_structure=network,
                                                       method_for_network_segragation=network_segragation, distance=distance, segragations = segragations, n = n )
                            if allocation_method == "ClashFullBypassKeyModulated" or allocation_method == "ClashFullBypassTransmissionModulated":
                                max_rate, segragations = max_rate
                        else:
                            max_rate = SKL_main_loop_2(main_params, adv_params, x, x0i, xb, ci, ni, f_atm, bounds, cons,
                                                       options, header, opt_head, fulldata, optdata, multidata, ground_station = ground_station,
                                                       satellite_trajectory_file=file_path + f"_{ground_station}_{RAAN}.csv",
                                                       allocation_method=allocation_method, file_path=file_path",
                                                       satellite_sites=list(self.site_height.keys()), satellite_RAAN=RAAN,
                                                       network_structure=network,
                                                       method_for_network_segragation=network_segragation, distance=distance, segragations = segragations, n = n)
                            if allocation_method == "ClashFullBypassKeyModulated" or allocation_method == "ClashFullBypassTransmissionModulated":
                                max_rate, segragations = max_rate
                    i += 1
                max_rates[RAAN][ground_station] = max_rate

        network.calculate_physical_distance_between_core_sites_and_stations()
        if network_segragation == "ConnectionBased":
            network.select_available_sites_based_on_connectivity()
        elif network_segragation == "DistanceBased":
            network.select_available_sites_based_on_distance(distance=distance)
        elif network_segragation == "SingleConnection":
            network.select_single_available_site()

        self.max_rates = max_rates
        self.Tij = Tij
        self.network = network
        self.ground_core_sites = ground_core_sites
        self.ground_station_sites = ground_station_sites


    def store_rate_satellite_ground_station_file(self, file_path, graph_id):
        rates_dictionary = []
        for RAAN in self.max_rates.keys():
            for station in self.max_rates[RAAN].keys():
                rates_dictionary.append({"ID": graph_id, "ground_core_station": station, "satellite_path": RAAN,
                                         "capacity": self.max_rates[RAAN][station]})

        dictionary_fieldnames = ["ID", "ground_core_station", "satellite_path","capacity"]
        if os.path.isfile(file_path + '.csv'):
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(rates_dictionary)
        else:
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(rates_dictionary)

    def store_simple_Tij_file(self, file_path, graph_id):
        ## Need to have initialised class with network_segragation == "SingleConnection" for correct functionality
        available_site_per_core = self.network.available_sites
        Tij = {}
        ground_stations_to_core_sites = {}
        for ground_station in self.ground_station_sites:
            ground_stations_to_core_sites[ground_station] = []
        for core_site in available_site_per_core:
            ground_stations_to_core_sites[available_site_per_core[core_site][0]].append(core_site)
        for ground_station in self.ground_station_sites:
            for ground_station_2 in self.ground_station_sites:
                if ground_station != ground_station_2:
                    Tij[ground_station, ground_station_2] = 0.0
                    for core_site in ground_stations_to_core_sites[ground_station]:
                        for core_site_2 in ground_stations_to_core_sites[ground_station_2]:
                            Tij[ground_station, ground_station_2] += self.Tij[core_site, core_site_2]

        Tij_dictionary = []
        for ground_station, ground_station_2 in Tij.keys():
            Tij_dictionary.append({"ID": graph_id, "ground_core_station_1": ground_station, "ground_core_station_2": ground_station_2,
                                         "Tij": Tij[ground_station, ground_station_2]})

        dictionary_fieldnames = ["ID", "ground_core_station_1", "ground_core_station_2", "Tij"]
        if os.path.isfile(file_path + '.csv'):
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(Tij_dictionary)
        else:
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(Tij_dictionary)



    def store_Tij_file(self, file_path, graph_id):
        Tij_dictionary = []
        for core_station, core_station_2 in self.Tij.keys():
            Tij_dictionary.append({"ID": graph_id, "core_site_1": core_station, "core_site_2": core_station_2,
                                         "Tij": self.Tij[core_station, core_station_2]})

        dictionary_fieldnames = ["ID", "core_site_1", "core_site_2", "Tij"]
        if os.path.isfile(file_path + '.csv'):
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(Tij_dictionary)
        else:
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(Tij_dictionary)

    def store_key_dict_file(self, file_path, graph_id):
        Tij =  []
        for ground_station in self.ground_station_sites:
            for ground_station_2 in self.ground_station_sites:
                if ground_station != ground_station_2:
                    if (ground_station_2, ground_station) not in Tij and (ground_station, ground_station_2) not in Tij:
                        Tij.extend([(ground_station_2, ground_station)])
        Tij_dictionary = []
        for core_station, core_station_2 in Tij:
            Tij_dictionary.append({"ID": graph_id, "ground_core_station_1": core_station, "ground_core_station_2": core_station_2})

        dictionary_fieldnames = ["ID", "ground_core_station_1", "ground_core_station_2"]
        if os.path.isfile(file_path + '.csv'):
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(Tij_dictionary)
        else:
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(Tij_dictionary)


    def store_core_site_keys(self,file_path, graph_id):
        max_len = 0
        for key in self.network.available_sites.keys():
            if len(self.network.available_sites[key]) > max_len:
                max_len = len(self.network.available_sites[key])
        core_site_keys = []
        for key in self.network.available_sites.keys():
            core_site_keys.append({"ID": graph_id, "core_site": key})
            for i in range(max_len):
                if len(self.network.available_sites[key])> i:
                    core_site_keys[-1][f"ground_station_allowed_{i}"] = self.network.available_sites[key][i]
                else:
                    core_site_keys[-1][f"ground_station_allowed_{i}"] = None

        dictionary_fieldnames = ["ID", "core_site"]
        for i in range(max_len):
            dictionary_fieldnames.append(f"ground_station_allowed_{i}")
        if os.path.isfile(file_path + '.csv'):
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writerows(core_site_keys)
        else:
            with open(file_path + '.csv', mode='a') as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=dictionary_fieldnames)
                writer.writeheader()
                writer.writerows(core_site_keys)



