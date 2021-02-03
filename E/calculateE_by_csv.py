#!/usr/bin/python3

# Goal: to calculate Resilient Modulus according data in csv files.

# import modules

import pandas as pd 
import numpy as np 
import statistics
import datetime as dt 
import os
import datetime

class Calculation:
    # initiate some variables
    def __init__(self, root_path, file_names, gauge_length, samp_interval, ):
        # load frequency
        self.freq = 0
        # sampling interval
        self.interval = samp_interval
        # offset summary
        self.offset = {"Folder": [], "Offset": [], }
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
        # empty dataframe for e data
        self.e_empty_dict = {
            "Sequence": [],
            "Acq Cycle": [],
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
        self.seq_timeseries = []
        # length of timeseries of a sequence (100 s)
        self.seq_timeseries_len = 0
        # timeseries of a cycle
        self.cycle_timeseries = []
        # length of timeseries of a cycle (100 s)
        self.cycle_timeseries_len = 0
        # sequence boundaries
        self.seq_point = []
        # start time
        self.start_time = datetime.datetime.now()
        # log dictionary
        self.log = {"Folder": [], }

    def cal_timeseries(self, folder, ):
        # get load frequency
        self.freq = int(folder.split("_")[4][0:-2])
        self.seq_timeseries = [
            time_i * self.interval for time_i in range(int(100 * (1/self.freq) / self.interval))
        ]
        # length of timeseries of a sequence (100 s)
        self.seq_timeseries_len = len(self.seq_timeseries)
        # timeseries of a cycle
        self.cycle_timeseries = [
            time_j * self.interval for time_j in range(int((1/self.freq) / self.interval))
        ]
        # length of timeseries of a cycle (100 s)
        self.cycle_timeseries_len = len(self.cycle_timeseries)
        # sequence boundaries
        self.seq_point = list(range(
            0, self.seq_timeseries_len * 15 + 1, self.seq_timeseries_len
        ))

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
        self.log["Combination"] = []
        # iterate fold names in keys of self.file_names
        count_i = 1
        for folder_i in self.file_names:
            # initiate start time
            self.set_time0()
            # calculate timeseries
            self.cal_timeseries(folder_i)
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
                            sort=False,
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

    def check_data(self, ):
        # print hint
        print("\nChecking Data...")
        self.log["Preprocessing"] = []
        count_j = 1
        for folder_j in self.file_names:
            self.set_time0()
            # calculate timeseries
            self.cal_timeseries(folder_j)
            # print hint
            print("Processing[{}/{}]>> {}".format(count_j, self.folder_no, folder_j, ))
            csv_ori_name = os.path.join(self.root_path, "{}_ori.csv".format(folder_j))
            try:
                pass
            except:
                # print hint
                print("!!! Preprocessing Data: {}".format(folder_j))
                KeyError
            pass
            # log duration
            self.log["Preprocessing"].append(self.duration(self.start_time))
            count_k += 1
        # print hint
        print("Preprocessing Data done!")
        pass

    def calculate_data(self, ):
        # print hint
        print("\nPreprocessing Data...")
        self.log["Preprocessing"] = []
        count_k = 1
        for folder_k in self.file_names:
            # initiate start time
            self.set_time0()
            # calculate timeseries
            self.cal_timeseries(folder_k)
            # print hint
            print("Processing[{}/{}]>> {}".format(count_k, self.folder_no, folder_k, ))
            csv_ori_name = os.path.join(self.root_path, "{}_ori.csv".format(folder_k), )
            # "folder_cal.csv"
            csv_cal_path = os.path.join(self.root_path, "{}_cal.csv".format(folder_k), )
            # calculate if "folder_cal.csv" doesnt exist
            if not os.path.isfile(csv_cal_path):  # not os.path.isfile(csv_cal_path)
                # create a copy of data frame
                temp1 = self.readcsv(csv_ori_name)
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
                    print("!!! Preprocessing Data: {}".format(folder_k))
                    KeyError
            # log duration
            self.log["Preprocessing"].append(self.duration(self.start_time))
            count_k += 1
        # print hint
        print("Preprocessing Data done!")

    def correct_offset(self, ):
        print("\nCorrecting Offset...")
        self.log["Correcting"] = []
        count_l = 1
        for folder_l in self.file_names:  # self.file_names, or ["RLTT_C16_100_2_1HZ_1"]
            # initiate start time
            self.set_time0()
            # calculate timeseries
            self.cal_timeseries(folder_l)
            self.offset["Folder"].append(folder_l )
            # print hint
            print("Processing[{}/{}]>> {}".format(count_l, self.folder_no, folder_l, ))
            # "folder_cal.csv"
            csv_cal_path = os.path.join(self.root_path, "{}_cal.csv".format(folder_l), )
            # make a copy of data frame as temp2
            temp2 = self.readcsv(csv_cal_path)
            # initiate an offsets list
            offset_list = []
            for seq_i in self.sequence:
                for cycle_i in range(1, 101):
                    # range of cycle i
                    start_pt = self.seq_point[seq_i - 1] + (cycle_i - 1) * self.cycle_timeseries_len
                    end_pt = self.seq_point[seq_i - 1] + cycle_i * self.cycle_timeseries_len
                    cycle_range = list(range(start_pt, end_pt, ))
                    # make a copy of present cycle data as temp3
                    temp3 = temp2["Cyclic Load"].iloc[cycle_range].copy().reindex()
                    # figure out max and min value location
                    min_loc = np.where(temp3 == temp3.min())[0][0]
                    max_loc = np.where(temp3 == temp3.max())[0][0]
                    # determine average offset of minimum load and maximum load
                    offset_list.append((min_loc + max_loc - self.cycle_timeseries_len/2) / 2)
            # determine the average offset of present expriment
            offset = int(round(statistics.mode(offset_list)))
            # adjust data postion
            # make a copy of the main data, and reset the index
            data_main = temp2.iloc[offset:, 1:].copy().reset_index(drop=True)
            # start point of the 15-seq 99-cycle
            seq_15_start = self.seq_timeseries_len * 14 + 99 * self.cycle_timeseries_len
            # make a copy of the interested data
            data_last = temp2.iloc[seq_15_start:seq_15_start+offset, 1:].copy().reset_index(drop=True)
            print("{}\n{}\n".format(data_main, data_last))
            data_main = pd.concat([data_main, data_last], ignore_index=True, )
            # assign some data of 14 sequence to the 15 sequence
            temp2.iloc[:, 1:] = data_main
            self.offset["Offset"].append(offset)
            # log duration
            self.log["Correcting"].append(self.duration(self.start_time))
            count_l += 1
            self.writecsv(temp2, csv_cal_path)

    def main(self, ):
        self.get_ori_data()
        # self.calculate_data()
        # self.correct_offset()
        # write log
        # self.writecsv(pd.DataFrame(self.log), "log.csv", )
    

# sampling interval
interval = 0.002

# get folder names
# Linux path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/E"
# Windows path = r"D:\TestData\RLTT\E"
root_path = r"D:\TestData\RLTT\E"  # "Mr" = Resilient Modulus, "E" = Complex Modulus
for roots_i, folders_i, files_i in os.walk(root_path):
    if len(folders_i) > 0:
        file_folders = folders_i

# Set particular folder name to analysis
# file_folders = ["RLTT_C13_95_1_3"]

# generate file abstract pathes
file_names = {}
for folder_i in file_folders:
    # generate folder abstract path
    file_folder = os.path.join(root_path, folder_i, )
    # generate file abstract pathes
    for roots_j, folers_j, files_j in os.walk(file_folder):
        # print hint
        print(file_folder, len(files_j), )
        # generate file names according to length of file list
        if len(files_j) != 0:
            file_name = [
                r"MF015- Unbound Material Resilient Modulus AASHTO T307_non-rest_trial_{}.csv".format(file_i) \
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
# # write filenames to filenames.csv
# pd.DataFrame(temp).to_csv("filenames.csv", index=False, sep=",", )

# generate gauge length
gauge_length = "gauge_length.csv"

demo = Calculation(root_path, file_names, gauge_length, interval, )
demo.main()
