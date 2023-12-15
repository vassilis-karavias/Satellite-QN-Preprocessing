import csv
import pandas as pd
import os
import numpy as np
def import_files(time_angle_file_path, time_rates_file_path):
    time_angle = pd.read_csv(time_angle_file_path)
    time_rate = pd.read_csv(time_rates_file_path)
    time_angle = time_angle[["# Time (s)", "Elevation (rad)"]]
    time_rate= time_rate[["# dt (s)", "SKL (b)"]]
    time_angle_positive = time_angle[time_angle["# Time (s)"] >= 0]
    return time_angle_positive, time_rate


def extrapolate_from_surrounding_data(time_rate, time_current, time_step):
    grad = max(0,(time_rate.loc[time_rate["# dt (s)"] == time_current + time_step]["SKL (b)"].values[0] - time_rate.loc[time_rate["# dt (s)"] == time_current - time_step]["SKL (b)"].values[0]) / (2 * time_step))
    return time_rate.loc[time_rate["# dt (s)"] == time_current - time_step]["SKL (b)"].values[0] + grad * time_step




def convert_to_elevation_rate(time_angle_positive, time_rate):
    time_step = time_rate["# dt (s)"][1]
    elevation_rate_dict = {}
    for i in range(len(time_angle_positive)):
        time = time_angle_positive.loc[i,"# Time (s)"]
        if time % time_step == 0:
            elevation = time_angle_positive.loc[i, "Elevation (rad)"] * 180 / np.pi
            if i != 0 and i - 1 < len(time_angle_positive):
                if time_rate.loc[time_rate["# dt (s)"] ==  time + time_step]["SKL (b)"].values[0] < 0.0001:
                    rate_extrapolated = extrapolate_from_surrounding_data(time_rate, time_current= time + time_step, time_step=time_step)
                    rate = max(0, (rate_extrapolated-time_rate.loc[time_rate["# dt (s)"] == time]["SKL (b)"].values[0]) / (2 * time_step))
                elif time_rate.loc[time_rate["# dt (s)"] ==  time]["SKL (b)"].values[0] < 0.0001:
                    rate_extrapolated = extrapolate_from_surrounding_data(time_rate, time_current=time , time_step=time_step)
                    rate = max(0, (time_rate.loc[time_rate["# dt (s)"] == time + time_step]["SKL (b)"].values[0] -
                                   rate_extrapolated) / (2 * time_step))
                else:
                    rate = max(0, (time_rate.loc[time_rate["# dt (s)"] ==  time + time_step]["SKL (b)"].values[0] - time_rate.loc[time_rate["# dt (s)"] ==  time]["SKL (b)"].values[0]) / (2 * time_step))
                elevation_rate_dict[elevation] = rate
    return elevation_rate_dict


def store_data_in_csv(elevation_rate_dict, file_path):
    dictionary = []
    for key in elevation_rate_dict.keys():
        dictionary.append({"Elevation (deg)": key, "Rate (b)": elevation_rate_dict[key]})
    if os.path.isfile(file_path):
        fieldnames = ["Elevation (deg)", "Rate (b)"]
        with open(file_path, mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writerows(dictionary)
    else:
        fieldnames = ["Elevation (deg)", "Rate (b)"]
        with open(file_path, mode='a') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(dictionary)


if __name__ == "__main__":
    time_angle_file_path = "FS_loss_th_m_90.00_wl_785nm_h_566.89km_h1_0.5km_aT_0.15m_aR_0.5m_w0_0.15m.csv"
    time_rates_file_path = "out_th_m_90.00_Pec_3e-06_QBERI_0.01_h1_0.5km.csv"
    elevation_rates_file_path = "elevation_rates_h_0.5km.csv"

    time_angle_positive, time_rate = import_files(time_angle_file_path, time_rates_file_path)
    elevation_rate_dict = convert_to_elevation_rate(time_angle_positive, time_rate)
    store_data_in_csv(elevation_rate_dict, file_path= elevation_rates_file_path)


