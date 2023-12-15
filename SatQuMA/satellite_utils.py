import numpy as np

##################################################################### Local Data #####################################

# sites_lat_lon = { "Athens": [37.9838, 23.7275], "Thessaloniki":[40.6401, 22.9444],"London": [51.5072, -0.1276],
#                   "Birmingham": [52.48652, -1.8904], "Manchester": [53.4808,-2.2426], "Leeds": [53.8008, -1.5491],
#                   "Newcastle": [54.9783, -1.6178], "Portsmouth": [50.8198, -1.0880], "Paris": [48.8566,2.3522],
#                   "Lyon": [45.7640,4.8357], "Lille": [50.6292, 3.0573], "Marseille": [43.2965,5.3698],
#                   "Nantes": [47.2184, -1.5538], "Frankfurt": [50.1109,8.6921], "Berlin": [52.5200, 13.4050],
#                   "Munich": [48.1351, 11.5820], "Dusseldorf": [51.2277, 6.7735], "Nurenburg": [49.4521,11.0767],
#                   "Leipzig": [51.3397,12.3731], "Stuttgart": [48.7758,9.1829], "Amsterdam": [52.3676,4.9041],
#                   "Madrid": [40.4168, -3.7038], "Barcelona": [41.3874,2.1686],"Malaga": [36.7213, -4.4213],
#                   "Valencia": [39.4699, -0.3763], "Porto": [41.1579,-8.6291],"Lisbon": [38.7223, -9.1393],
#                   "Zurich": [47.3769,8.5417], "Geneva": [46.2044, 6.1432], "Bern": [46.9480, 7.4474],
#                   "Lucerne": [47.0502, 8.3093], "Milan": [45.4642,9.19000], "Florence": [43.7696, 11.2558],
#                   "Venice": [45.4408, 12.3155], "Rome": [41.9028, 12.4964], "Brussels": [50.8476,4.3572],
#                   "Luxembourg": [49.8153, 6.1296], "Stockholm": [59.3293, 18.0686], "Copenhagen": [55.6761, 12.5683],
#                   "Oslo": [59.9139, 10.7522], "Sofia": [42.6977, 23.3219], "Bucharest": [44.4268, 26.1025]}
#
#
#
# # sites_lat_lon = {"London": [51.5072, -0.1276], "Madrid": [40.4168, -3.7038], "Edinburgh": [55.9533, -3.188267],
# #                  "Barcelona": [41.3874,2.1686], "Paris": [48.8566,2.3522], "Nantes": [47.2184, -1.5538],
# #                  "Lyon": [45.7640,4.8357], "Lisbon": [38.7223, -9.1393], "Bern": [46.9480, 7.4474],
# #                  "Freiburg": [47.9990,7.8421], "Milan": [45.4642,9.19000], "Florence": [43.7696, 11.2558],
# #                  "Naples": [40.8518, 14.2681], "Frankfurt": [50.1109,8.6921], "Dortmund": [51.5136,7.4653],
# #                  "Berlin": [52.5200, 13.4050], "Athens": [37.9838, 23.7275], "Ioannina": [39.6650,20.8537],
# #                  "Bristol": [51.4545,-2.5879], "Cambridge": [52.1951,0.1313], "Cardiff": [51.4837,-3.1681],
# #                  "Liverpool": [53.4084,-2.9916], "Birmingham": [52.48652, -1.8904], "Glasgow": [55.8617,-4.2583],
# #                  "Orleans": [47.9030, 1.9093], "Bordeaux": [44.8378, -0.5792], "Toulouse": [43.6047,1.4442],
# #                  "Marseille": [43.2965,5.3698], "Bilbao": [43.2630, -2.9350], "Seville": [37.3891, -5.9845],
# #                  "Malaga": [36.7213, -4.4213], "Valencia": [39.4699, -0.3763], "Zaragoza": [41.6488, -0.8891],
# #                  "Porto": [41.1579,-8.6291], "Vigo": [42.2406, -8.7270], "Munich": [48.1351, 11.5820],
# #                  "Monaco": [43.7384, 7.4246], "Bern": [46.9480, 7.4474], "Geneva": [46.2044, 6.1432],
# #                  "Brussels": [50.8476,4.3572], "Amsterdam": [52.3676,4.9041], "Zurich": [47.3769,8.5417],
# #                  "Stuttgart": [48.7758,9.1829], "Nuremberg": [49.4521,11.0767], "Hamburg": [53.5488,9.9872],
# #                  "Bremen": [53.0793, 8.8017], "Leipzig": [51.3397,12.3731], "Cologne": [50.9375, 6.9603],
# #                  "Genoa": [44.4056,8.9463], "Bologna": [44.4949, 11.3426], "Bari": [41.1171,16.8719],
# #                  "Nice": [43.7102,7.2620], "Ancona": [43.6158,13.5189], "Kalamata": [37.0366,22.1144],
# #                  "Patras": [38.2466,21.7346], "Volos": [39.3666, 22.9507], "Thessaloniki":[40.6401, 22.9444],
# #                  "Kavala": [40.9376, 24.4129]
# #                  }
#
#
# ground_station_potential_sites = {"London": [51.5072, -0.1276], "Madrid": [40.4168, -3.7038], "Edinburgh": [55.9533, -3.188267],
#                  "Barcelona": [41.3874,2.1686], "Paris": [48.8566,2.3522], "Nantes": [47.2184, -1.5538],
#                  "Lyon": [45.7640,4.8357], "Lisbon": [38.7223, -9.1393], "Bern": [46.9480, 7.4474],
#                  "Freiburg": [47.9990,7.8421], "Milan": [45.4642,9.19000], "Florence": [43.7696, 11.2558],
#                  "Naples": [40.8518, 14.2681], "Frankfurt": [50.1109,8.6921], "Dortmund": [51.5136,7.4653],
#                  "Berlin": [52.5200, 13.4050], "Athens": [37.9838, 23.7275], "Ioannina": [39.6650,20.8537]}
#
# mainland_gs_sites = ["Madrid", "Barcelona", "Paris", "Nantes","Lyon", "Lisbon", "Bern","Freiburg", "Milan", "Florence",
#                  "Naples", "Frankfurt", "Dortmund","Berlin"]
# uk_gs = ["London", "Edinburgh"]
# greece_gs = ["Athens", "Ioannina"]
# connections_allowed = {"Athens": greece_gs, "Thessaloniki":greece_gs,"London": uk_gs,
#                   "Birmingham": uk_gs, "Manchester": uk_gs, "Leeds": uk_gs,
#                   "Newcastle": uk_gs, "Portsmouth": uk_gs, "Paris": mainland_gs_sites,
#                   "Lyon": mainland_gs_sites, "Lille": mainland_gs_sites, "Marseille": mainland_gs_sites,
#                   "Nantes": mainland_gs_sites, "Frankfurt": mainland_gs_sites, "Berlin": mainland_gs_sites,
#                   "Munich": mainland_gs_sites, "Dusseldorf":mainland_gs_sites, "Nurenburg":mainland_gs_sites,
#                   "Leipzig": mainland_gs_sites, "Stuttgart": mainland_gs_sites, "Amsterdam": mainland_gs_sites,
#                   "Madrid": mainland_gs_sites, "Barcelona": mainland_gs_sites,"Malaga": mainland_gs_sites,
#                   "Valencia": mainland_gs_sites, "Porto": mainland_gs_sites,"Lisbon": mainland_gs_sites,
#                   "Zurich": mainland_gs_sites, "Geneva": mainland_gs_sites, "Bern": mainland_gs_sites,
#                   "Lucerne": mainland_gs_sites, "Milan": mainland_gs_sites, "Florence": mainland_gs_sites,
#                   "Venice": mainland_gs_sites, "Rome": mainland_gs_sites, "Brussels": mainland_gs_sites,
#                   "Luxembourg": mainland_gs_sites, "Stockholm": mainland_gs_sites, "Copenhagen": mainland_gs_sites,
#                   "Oslo": mainland_gs_sites, "Sofia": greece_gs, "Bucharest": greece_gs}
#
# datacentres_site = {"Athens": 16, "Thessaloniki":5,"London": 147,
#                   "Birmingham": 39, "Manchester": 36, "Leeds": 11,
#                   "Newcastle": 12, "Portsmouth": 18, "Paris": 59,
#                   "Lyon": 20, "Lille": 16, "Marseille": 29,
#                   "Nantes": 20, "Frankfurt": 72, "Berlin": 22,
#                   "Munich": 20, "Dusseldorf":39, "Nurenburg":23,
#                   "Leipzig": 18, "Stuttgart": 17, "Amsterdam": 122,
#                   "Madrid": 26, "Barcelona": 21,"Malaga": 9,
#                   "Valencia": 14, "Porto": 12,"Lisbon": 21,
#                   "Zurich": 36, "Geneva": 10, "Bern": 7,
#                   "Lucerne": 17, "Milan": 45, "Florence": 13,
#                   "Venice": 13, "Rome": 11, "Brussels": 40,
#                   "Luxembourg": 15, "Stockholm": 43, "Copenhagen": 37,
#                   "Oslo": 24, "Sofia": 27, "Bucharest": 25}



##################################################################### Global Data #####################################

ground_station_potential_sites = {"London": [51.5072, -0.1276], "Frankfurt": [50.1109, 8.6821], "Graz": [47.0707, 15.4395], "Johanessburg": [-26.2041, 28.0473],
                       "Sao Paulo": [-23.5558, -46.6496], "Tokyo": [35.653, 139.839], "Auckland": [-36.8509, 174.7645],
                       "New Delhi": [28.6139, 77.2090], "Mumbai": [19.0760, 72.8777], "Bangalore": [12.9716, 77.5946],
                       "Xinglong": [40.395833, 117.5775], "Nanshan": [43.475278, 87.176667], "Perth": [-31.9523, 115.8613],
                       "Brisbane": [-27.4705, 153.0260], "Melborne": [-37.8136, 144.9631], "Albany": [42.6526, -73.7562],
                       "Denver": [39.7392, -104.9903], "Nashville": [36.1627, -86.7816]}

datacentres_site = {"London": 151, "Manchester": 66, "Amsterdam": 162, "Stockholm": 127, "Zurich": 83, "Frankfurt": 173,
                    "Berlin": 66, "Paris": 111, "Madrid": 112, "Milan": 105, "Johanessburg": 20, "Sao Paulo": 75,
                    "Buenos Aires": 22, "Tokyo": 54, "Auckland": 35, "New Delhi": 30, "Mumbai": 52, "Bangalore": 50,
                   "Beijing":  29, "Shanghai": 33, "Hong Kong": 16, "Perth": 18, "Adelade": 18, "Melborne": 30, "Sydney": 36,
                    "Brisbane": 25, "LA": 331, "Seattle": 157, "Dallas": 285, "Chicago": 259, "Washington": 635, "Jacksonville": 187}


western_europe = ["London", "Manchester", "Amsterdam", "Stockholm", "Zurich", "Frankfurt", "Berlin", "Paris", "Madrid",
                  "Milan"]
africa = ["Johanessburg"]
latin_america = ["Sao Paulo", "Buenos Aires"]
asia_pacific = ["Tokyo", "Auckland", "New Delhi", "Mumbai", "Bangalore",
                   "Beijing", "Shanghai", "Hong Kong", "Perth", "Adelade", "Melborne", "Sydney",
                    "Brisbane"]
north_america = ["LA", "Seattle", "Dallas", "Chicago", "Washington", "Jacksonville"]

datacentres_site_20_years = {}
for key in western_europe:
    datacentres_site_20_years[key] = datacentres_site[key] * (1+0.95 * 4)
for key in africa:
    # technically the prediction in this region is 10 times growth in 5 years but this is likely a
    # result of the small initial number of datacentres so we cap the growth to 3 times for the first
    # 20 years and
    datacentres_site_20_years[key] = datacentres_site[key] * (1+2*4)
for key in latin_america:
    datacentres_site_20_years[key] = datacentres_site[key] * (1+1.23*4)
for key in asia_pacific:
    datacentres_site_20_years[key] = datacentres_site[key] * (1 + 1.43*4)
for key in north_america:
    datacentres_site_20_years[key] = datacentres_site[key] * (1 + 0.36*4)

datacentres_site_10_years = {}
for key in western_europe:
    datacentres_site_10_years[key] = datacentres_site[key] * (1+0.95*2)
for key in africa:
    # technically the prediction in this region is 10 times growth in 5 years but this is likely a
    # result of the small initial number of datacentres so we cap the growth to 3 times for the first
    # 20 years and
    datacentres_site_10_years[key] = datacentres_site[key] * (1+2*2)
for key in latin_america:
    datacentres_site_10_years[key] = datacentres_site[key] * (1+1.23*2)
for key in asia_pacific:
    datacentres_site_10_years[key] = datacentres_site[key] * (1 + 1.43*2)
for key in north_america:
    datacentres_site_10_years[key] = datacentres_site[key] * (1 + 0.36*2)

datacentres_site_30_years = {}
for key in western_europe:
    datacentres_site_30_years[key] = datacentres_site[key] * (1+0.95*6)
for key in africa:
    # technically the prediction in this region is 10 times growth in 5 years but this is likely a
    # result of the small initial number of datacentres so we cap the growth to 3 times for the first
    # 20 years and
    datacentres_site_30_years[key] = datacentres_site[key] * (1+2* 6)
for key in latin_america:
    datacentres_site_30_years[key] = datacentres_site[key] * (1+1.23* 6)
for key in asia_pacific:
    datacentres_site_30_years[key] = datacentres_site[key] * (1 + 1.43*6)
for key in north_america:
    datacentres_site_30_years[key] = datacentres_site[key] * (1 + 0.36*6)

europe_gs = ["London", "Frankfurt", "Graz"]
china_gs = ["Xinglong", "Nanshan"]
india_gs= ["New Delhi", "Mumbai", "Bangalore"]
austrailia_gs = ["Perth", "Brisbane", "Melborne"]
us_gs = ["Albany", "Denver", "Nashville"]
connections_allowed = {"London": europe_gs, "Manchester": europe_gs, "Amsterdam": europe_gs, "Stockholm": europe_gs, "Zurich": europe_gs, "Frankfurt": europe_gs,
                    "Berlin": europe_gs, "Paris": europe_gs, "Madrid": europe_gs, "Milan": europe_gs, "Johanessburg": ["Johanessburg"], "Sao Paulo": ["Sao Paulo"],
                    "Buenos Aires": ["Sao Paulo"], "Tokyo": ["Tokyo"], "Auckland": ["Auckland"], "New Delhi": india_gs, "Mumbai": india_gs, "Bangalore": india_gs,
                   "Beijing":  china_gs, "Shanghai": china_gs, "Hong Kong": china_gs, "Perth": austrailia_gs, "Adelade": austrailia_gs, "Melborne": austrailia_gs,
                    "Sydney": austrailia_gs,
                    "Brisbane": austrailia_gs, "LA": us_gs, "Seattle": us_gs, "Dallas": us_gs, "Chicago": us_gs, "Washington": us_gs, "Jacksonville": us_gs}

sites_lat_lon = {"London": [51.5072, -0.1276], "Manchester": [53.4808,-2.2426], "Amsterdam": [52.3676,4.9041],
                 "Stockholm": [59.3293, 18.0686], "Zurich": [47.3769,8.5417], "Frankfurt": [50.1109,8.6921],
                 "Berlin": [52.5200, 13.4050], "Paris": [48.8566,2.3522], "Madrid": [40.4168, -3.7038],
                  "Milan": [45.4642,9.19000], "Johanessburg": [-26.2041, 28.0473], "Sao Paulo": [-23.5558, -46.6496],
                    "Buenos Aires": [-34.6037, -58.3816], "Tokyo": [35.653, 139.839], "Auckland": [-36.8509, 174.7645],
                  "New Delhi": [28.6139, 77.2090], "Mumbai": [19.0760, 72.8777], "Bangalore": [12.9716, 77.5946],
                   "Beijing":  [39.9042, 116.4074], "Shanghai": [31.2304, 121.4737], "Hong Kong": [22.3193, 114.1694], "Perth": [-31.9523, 115.8613],
                 "Adelade": [-34.9285, 138.6007], "Melborne": [-37.8136, 144.9631], "Sydney": [-33.8688, 151.2093],
                    "Brisbane": [-27.4705, 153.0260], "LA": [34.055, -118.24], "Seattle": [47.6061, -122.3328], "Dallas": [32.7767, -96.7970],
                 "Chicago": [41.8781, 87.6298], "Washington": [38.9072, -77.0369], "Jacksonville": [30.3322, -81.6557]}

############################################################################################################################################################


def get_ground_core_stations(core_station_list =None, ground_station_list = None):
    if ground_station_list != None:
        connections_allowed_new_dict = {}
        for key in connections_allowed.keys():
            current_connections = connections_allowed[key]
            connections_allowed_new = []
            for gs in current_connections:
                if gs in ground_station_list:
                    connections_allowed_new.append(gs)
            connections_allowed_new_dict[key] = connections_allowed_new
    else:
        connections_allowed_new_dict = connections_allowed
    if core_station_list != None:
        ground_core_stations = {}
        for key in sites_lat_lon.keys():
            if key in core_station_list:
                if len(connections_allowed_new_dict[key]) >0:
                    ground_core_stations[key] = [np.radians(sites_lat_lon[key][0])
                                        , np.radians(sites_lat_lon[key][1]), connections_allowed_new_dict[key]]
    else:
        ground_core_stations = {}
        for key in sites_lat_lon.keys():
            if len(connections_allowed_new_dict[key]) > 0:
                ground_core_stations[key] = [np.radians(sites_lat_lon[key][0]), np.radians(sites_lat_lon[key][1]),
                                             connections_allowed_new_dict[key]]
    return ground_core_stations


def get_satellite_ground_stations(ground_station_list = None):
    if ground_station_list != None:
        satellite_ground_stations = {}
        for key in ground_station_potential_sites.keys():
            if key in ground_station_list:
                satellite_ground_stations[key] = [np.radians(ground_station_potential_sites[key][0]), np.radians(ground_station_potential_sites[key][1])]
        return satellite_ground_stations
    else:
        satellite_ground_stations = {}
        for key in ground_station_potential_sites.keys():
                satellite_ground_stations[key] = [np.radians(ground_station_potential_sites[key][0]),
                                                  np.radians(ground_station_potential_sites[key][1])]
        return satellite_ground_stations

def get_Tij(core_station_list = None):
    if core_station_list != None:
        Tij = {}
        for core_station in datacentres_site.keys():
            if core_station in core_station_list:
                for core_station_2 in datacentres_site.keys():
                    if core_station_2 in core_station_list and core_station_2 != core_station:
                        Tij[(core_station, core_station_2)] = datacentres_site[core_station] * datacentres_site[core_station_2] * 0.02
    else:
        Tij = {}
        for core_station in datacentres_site.keys():
            for core_station_2 in datacentres_site.keys():
                if core_station_2 != core_station:
                    Tij[(core_station, core_station_2)] = datacentres_site[core_station] * datacentres_site[core_station_2] * 0.02
    return Tij
