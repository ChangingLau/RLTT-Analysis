#!/usr/bin/python3

# Goal: to calculate Resilient Modulus according data in csv files.

# import modules

import pandas as pd 
import numpy as np 
import datetime as dt 
import os
import datetime
import copy

class Calculation:
    # initiate some variables
    def __init__(self, root_path, file_names, gauge_length, samp_interval, ):
        # sampling interval
        self.interval = samp_interval
        # root path
        self.root_path = root_path
        # empty dataframe for original data
        self.ori_empty_df = pd.DataFrame(
            {
                "Sequence": [],
                "Acq Cycle": [],
                "Time": [],
                "Actuator Displacement": [],
                "Load": [],
                "Pressure": [],
                "Axial LVDT1": [],
                "Axial LVDT2": [],
            }
        )
        # empty dataframe for mr data
        self.mr_empty_dict = {
            "Sequence": [],
            "Acq Cycle": [],
            "Pressure": [],
            "Cyclic Axial Load": [],
            "Cyclic Axial Stress": [],
            "Axial Permanent Deform": [],
            "Axial Permanent Strain": [],
            "Axial Resilient Deform": [],
            "Axial Resilient Strain": [],
            "Axial Resilient Modulus": [],
        }
        # file names, a dictionary, key = folder_name, value = list of filenames
        self.file_names = file_names
        # folder number
        self.folder_no = len(self.file_names)
        # gauge length of each specimen
        self.gauge_length = pd.read_csv(gauge_length, sep=",", index_col=0)
        self.gauge_length["Height"].astype('float64')  # set column value as float
        # area of top pad, area = pi * r ^ 2
        self.pad_area = np.pi * 0.1 ** 2 / 4
        # sequence list
        self.sequence = list(range(1, 16))
        # timeseries of a sequence
        self.seq_timeseries = [
            time_i * self.interval for time_i in range(int(100 / self.interval))
        ]
        # length of timeseries of a sequence (100 s)
        self.seq_timeseries_len = len(self.seq_timeseries)
        # timeseries of a cycle
        self.cycle_timeseries = [
            time_j * self.interval for time_j in range(int(1 / self.interval))
        ]
        # length of timeseries of a cycle (100 s)
        self.cycle_timeseries_len = len(self.cycle_timeseries)
        # sequence boundaries
        self.seq_point = list(range(
            0, self.seq_timeseries_len * 15 + 1, self.seq_timeseries_len
        ))
        # start time
        self.start_time = datetime.datetime.now()
        # log dictionary
        self.log = {
            "Folder": [],
            "Combination": [],
            "CalDeform":[],
            "CalMr":[],
        }

    def set_time0(self, ):
        # initiate start_time
        self.start_time = datetime.datetime.now()

    def duration(self, time0, ):
        # determine duration
        return datetime.datetime.now() - time0

    def readcsv(self, ori, sep=",", ):
        # read csv and return it as dataframe 
        return pd.read_csv(ori, sep=sep, )  

    def writecsv(self, content, target, index=False, ):
        # write csv with dataframe
        csv_name = target.split(".")
        content.to_csv("{}.{}".format(csv_name[0], csv_name[-1]), index=index, )  

    def get_ori_data(self, ):
        print("\nCombining Data...")
        # iterate fold names in keys of self.file_names
        count_i = 1
        for folder_i in self.file_names:
            # initiate start time
            self.set_time0()
            # initiate an empty data frame for a new file folder
            temp0 = self.ori_empty_df.copy()
            # print hint
            print("Processing[{}/{}]>> {}".format(count_i, self.folder_no, folder_i, ))
            # log
            self.log["Folder"].append(folder_i)
            # iterate file names in value of dict self.file_names
            csv_ori_name = os.path.join(self.root_path, "{}_ori.csv".format(folder_i), )
            # combine csv files if "folder.csv" doesnt exist
            if not os.path.isfile(csv_ori_name):  # not os.path.isfile(csv_ori_name)
                for file_name_i in self.file_names[folder_i]:
                    try:
                        # try read data from csv file, and concate it with temp0
                        temp0 = pd.concat(
                            [temp0, self.readcsv(file_name_i).iloc[1:, :]], 
                            ignore_index=True,
                            ).astype('float64')
                    except:
                        # print hint
                        print("!!! Combining Data: {}, {}".format(folder_i, file_name_i, ))
                        FileNotFoundError
                # write files into csv files, and named with fold names
                self.writecsv(pd.DataFrame(temp0), csv_ori_name)
            # log duration
            self.log["Combination"].append(self.duration(self.start_time))
            count_i += 1
        # print hint
        print("File collection done!")

    def calculate_data(self, ):
        # print hint
        print("\nClaculating Deformation...")
        count_j = 1
        for folder_j in self.file_names:
            # initiate start time
            self.set_time0()
            # print hint
            print("Processing[{}/{}]>> {}".format(count_j, self.folder_no, folder_j, ))
            csv_ori_path = os.path.join(self.root_path, "{}_ori.csv".format(folder_j), )
            # "folder_cal.csv"
            csv_cal_path = os.path.join(self.root_path, "{}_cal.csv".format(folder_j), )
            # calculate if "folder_cal.csv" doesnt exist
            if not os.path.isfile(csv_cal_path):  # not os.path.isfile(csv_cal_path)
                # create a copy of data frame
                temp1 = self.readcsv(csv_ori_path).iloc[self.seq_timeseries_len * 5:, :]
                try:
                    # determine the compensation for confining pressure
                    # compensation = confining pressure * area of axial rod
                    temp1["Compensation"] = np.pi * 0.025 ** 2 / 4 * temp1["Pressure"].copy()
                    # determine the max axial load 
                    # max axial load = load - compensation
                    temp1["Max Axial Load"] = temp1["Load"].copy() - temp1["Compensation"].copy()
                    # determine the cyclic load
                    # cyclic load = 0.9 * max axial load
                    temp1["Cyclic Load"] = 0.9 * temp1["Max Axial Load"].copy()
                    # write files into csv files, and named with fold names
                    self.writecsv(pd.DataFrame(temp1), csv_cal_path)
                except:
                    # print hint
                    print("!!! Calculate Data: {}".format(folder_j))
                    KeyError
            # log duration
            self.log["CalDeform"].append(self.duration(self.start_time))
            count_j += 1
        # print hint
        print("Claculating Deformation done!")

    def calculate_mr(self, ):
        print("\nClaculating Mr...")
        count_k = 1
        for folder_k in self.file_names:
            # initiate start time
            self.set_time0()
            print("Processing[{}/{}]>> {}".format(count_k, self.folder_no, folder_k, ))
            # specimen number
            specimen_no = "{}_{}_{}_{}".format(*folder_k.split("_")[0: 4])
            # specimen height
            specimen_h = self.gauge_length.loc[specimen_no, "Height",]
            # "folder_cal.csv"
            csv_cal_path = os.path.join(self.root_path, "{}_cal.csv".format(folder_k), )
            # "folder_mr.csv"
            csv_mr_path = os.path.join(self.root_path, "{}_mr.csv".format(folder_k), )
            # make a copy of data frame as temp2
            temp2 = self.readcsv(csv_cal_path).copy()
            # initiate column name of Deformation
            temp2["Deformation"] = temp2["Time"].copy()
            temp3 = copy.deepcopy(self.mr_empty_dict)
            # calculate deformation and others
            for seq_i in self.sequence:
                for cycle_i in range(1, 101):
                    # range of cycle i
                    start_pt = self.seq_point[seq_i - 1] + (cycle_i - 1) * self.cycle_timeseries_len
                    end_pt = self.seq_point[seq_i - 1] + cycle_i * self.cycle_timeseries_len
                    cycle_range = list(range(start_pt, end_pt, ))
                    # calculate deformation
                    try:
                        # sum initial value of LVDT1 and LVDT2
                        initial_LVDT = temp2["Axial LVDT1"].iloc[start_pt].copy()\
                            + temp2["Axial LVDT2"].iloc[start_pt].copy()
                        temp2["Deformation"].iloc[cycle_range] = \
                            ((temp2["Axial LVDT1"].iloc[cycle_range].copy() + \
                                temp2["Axial LVDT2"].iloc[cycle_range].copy()) \
                                    - initial_LVDT) / 2
                    except:
                        # print hint
                        print("!!! Calculate Deformation: {}, {}, {}".format(
                            folder_k, seq_i, cycle_i
                        ))
                        KeyError
                    # calculate Mr related data
                    try:
                        # maximum deformation of the cycle
                        max_deform = temp2["Deformation"].iloc[cycle_range].max().copy()
                        # pressure
                        pressure = temp2["Pressure"].iloc[cycle_range].mean().copy()
                        temp3["Pressure"].append(pressure)
                        # cyclic axial load
                        cyclic_load = temp2["Cyclic Load"].iloc[cycle_range].max().copy()
                        temp3["Cyclic Axial Load"].append(cyclic_load)
                        # cyclic axial stress
                        cyclic_stress = cyclic_load / self.pad_area
                        temp3["Cyclic Axial Stress"].append(cyclic_stress)
                        # axial permanent deformation
                        permanent_deform = temp2["Deformation"].iloc[cycle_range].iloc[-1].copy()
                        temp3["Axial Permanent Deform"].append(permanent_deform)
                        # axial permanent strain = axial permanent deformation / gauge length
                        temp3["Axial Permanent Strain"].append(permanent_deform / specimen_h)
                        # axial resilient deformation
                        resilient_deform = max_deform - permanent_deform
                        temp3["Axial Resilient Deform"].append(resilient_deform)
                        # axial resilient strain = axial resilient deformation / gauge length
                        resilient_strain = resilient_deform / specimen_h
                        temp3["Axial Resilient Strain"].append(resilient_strain)
                        # axial resilient modulus = cyclic axial stress / axial resilient strain
                        resilient_modulus = cyclic_stress / resilient_strain / 1000
                        temp3["Axial Resilient Modulus"].append(resilient_modulus)
                        # sequence number
                        temp3["Sequence"].append(seq_i)
                        # cycle number
                        temp3["Acq Cycle"].append(cycle_i)
                    except:
                        # print hint
                        print("!!! Calculate Mr: {}, {}, {}".format(
                            folder_k, seq_i, cycle_i
                        ))
                        IndexError
            # log duration
            self.log["CalMr"].append(self.duration(self.start_time))
            count_k += 1
            # write files into csv files, and named with fold names
            self.writecsv(pd.DataFrame(temp3), csv_mr_path)

    def main(self, ):
        self.get_ori_data()
        self.calculate_data()
        self.calculate_mr()
        # write log
        self.writecsv(pd.DataFrame(self.log), "log.csv", )
    

# sampling interval
interval = 0.002

# get folder names
# Linux path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"
# Windows path = r"D:\TestData\RLTT\Mr"
root_path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"  # "Mr" = Resilient Modulus, "E" = Complex Modulus
for roots_i, folders_i, files_i in os.walk(root_path):
    if len(folders_i) > 0:
        file_folders = folders_i
        break

# Set particular folder name to analysis for tunning
# file_folders = ["RLTT_C16_95_1_1", "RLTT_C16_95_1_2", ]

# generate file abstract pathes
file_names = {}
for folder_i in file_folders:
    # generate folder abstract path
    file_folder = os.path.join(root_path, folder_i, )
    # generate file abstract pathes
    for roots_j, folers_j, files_j in os.walk(file_folder):
        # print hint
        print(file_folder, len(files_j))
        # generate file names according to length of file list
        if len(files_j) != 0:
            pass
            file_name = [
                r"MF015- Unbound Material Resilient Modulus AASHTO T307_trial_{}.csv".format(file_i) \
                    for file_i in range(1, len(files_j)+1)
            ]
            # generate file abstract pathes
            file_name = [os.path.join(file_folder, file_j) for file_j in file_name]
            # assign file abstract pathes into dictionary
            file_names[folder_i] = file_name


# write all files' abstract pathes into a csv file
# collect all files' abstract pathes into a list
temp = {"Filename": []}
for key_i in file_names:
    temp["Filename"] += file_names[key_i][:]
# print statistic result
print("SUM: #folders: {}, #CSV files:{}".format(
    len(file_folders), len(temp["Filename"])
))
# write filenames to filenames.csv
# pd.DataFrame(temp).to_csv("filenames.csv", index=False, sep=",", )

# generate gauge length
gauge_length = "gauge_length.csv"

demo = Calculation(root_path, file_names, gauge_length, interval, )
demo.main()
