#!/usr/bin/env python3

# import some necessary modules
import pandas as pd 
import numpy as np 
import os
import copy
import datetime
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import time


"""
=========================================================
Mr Models
=========================================================
[0]- Seed Model: Mr = k1 * (theta / pa) ^ k2
[1]- Uzan & Witczak Model: Mr = k1 * (theta / pa) ^ k2 * (sigmaD / pa) ^ k3
[2]- NCHRP Model: Mr = k1 * (theta / pa) ^ k2 * (tau / pa + 1) ^ k3
=========================================================
Function
=========================================================
Regression with last 5 seqs
"""

# define a class
class DoRegression:
    def __init__(self, root_path, file_names, ):
        # air pressure
        self.pa = 1
        # root path
        self.root_path = root_path
        # file names list
        self.file_names = file_names
        # folder number
        self.folder_no = len(self.file_names)
        # columns
        self.columns = ["Folder"] + \
                ["Seq{}".format(seq_i) for seq_i in list(range(1, 16))]
        # result dictionary
        self.result = {}
        for col_i in self.columns:
            self.result[col_i] = []
        # start time
        self.start_time = datetime.datetime.now()
        # log dictionary
        self.log = {
            "Folder": [],
            "Combination": [],
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

    def calculation(self, folder, data, no=5):
        # append folder
        self.result[self.columns[0]].append(folder)
        for seq_i in range(1, 16):
            cycle_range = list(range((seq_i)*100-no, (seq_i)*100))
            # append average value of each sequence
            self.result[self.columns[seq_i]].append(
                data["Axial Resilient Modulus"].iloc[cycle_range].mean()
            )

    def main(self, ):
        # print hint
        print("Combinating...")
        count_i = 1
        for folder_i in self.file_names:
            # print hint
            print("Processing[{}/{}]>> {}".format(count_i, self.folder_no, folder_i, ))
            # initiate start time
            self.set_time0()
            # log
            self.log["Folder"].append(folder_i)
            # generate mr file path and read file
            csv_mr_path = os.path.join(self.root_path, "{}_mr.csv".format(folder_i))
            temp1 = self.readcsv(csv_mr_path, ).copy()
            # input data columns = 
            # [Sequence, Acq Cycle, Cyclic Axial Load, Cyclic Axial Stress, Axial Permanent Deform, \
            # Axial Permanent Strain, Axial Resilient Deform, Axial Resilient Strain, Axial Resilient Modulus]
            # calculating some necessary parameters
            try:
                self.calculation(folder_i, temp1, )
            except:
                # print hint
                print("!!! Combinating: {}".format(folder_i, ))
            # log duration
            self.log["Combination"].append(self.duration(self.start_time))
            count_i += 1
        # output result
        self.writecsv(pd.DataFrame(self.result), "overall_last_5_cycle_ave_Mr.csv", )
        # write log
        self.writecsv(pd.DataFrame(self.log), "log.csv", )


# generate file names
# file path
# Linux path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"
# Windows path = r"D:\TestData\RLTT\Mr"
root_path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"
# get all files
for i, j, k in os.walk(root_path):
    if len(j) != 0:
        temp0 = k
# get all "*mr.csv" files
file_names = []
for file_i in temp0:
    if file_i.split("_")[-1].split(".")[0] == "mr":
        file_names.append(file_i.split("_mr")[0])

# Set particular folder name to analysis
# file_names = ["RLTT_C16_95_1_3"]

demo = DoRegression(root_path, file_names, )
demo.main()
