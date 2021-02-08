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
[1]- Uzan & Witczak Model: Mr = k1 * (theta / pa) ^ k2 * (tau / pa) ^ k3
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
        # rows to be dropped
        self.rows2drop = []
        # rows list for train
        self.rows_for_train = []
        for i in range(0, 95):
            self.rows_for_train += list(range(i, 1500+i, 100))
        # rows list for predict
        self.rows_for_predict = []
        for j in range(95, 100):
            self.rows_for_predict += list(range(j, 1500+j, 100))
        # input data columns
        # Sequence, Acq Cycle, Cyclic Axial Load, Cyclic Axial Stress, Axial Permanent Deform, \
        # Axial Permanent Strain, Axial Resilient Deform, Axial Resilient Strain, Axial Resilient Modulus
        self.parameters = {
            "Axial Resilient Modulus": "", "sigma3": "", "sigmaD": "",
            "theta": "", "tau": "", "theta / pa": "",
            "tau / pa": "", "tau / pa + 1": "",
        }
        # information about models
        self.model_xs = {
            # model: [[parameters], [coefficients], [intercepts], ]
            "Seed_Model": ["theta / pa"], 
            "Uzan&Witczak_Model": ["theta / pa", "tau / pa", ], 
            "NCHRP_Model": ["theta / pa", "tau / pa + 1", ], 
        }
        # regression result
        self.regr_result = {
            "coef_": [],
            "intercept_": [], 
            "MSE": [],
            "R2": [],
        }
        # result sum 
        self.result_sum = {
            "Folder": [],
            "Model": [],
            "coef_": [],
            "intercept_": [], 
            "MSE": [],
            "R2": [],
        }
        # start time
        self.start_time = datetime.datetime.now()
        # log dictionary
        self.log = {
            "Folder": [],
            "Regression": [],
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

    def calculation(self, mr, pressure, axial_cyclic_stress, ):
        # calculate some basic parameters
        paras = copy.deepcopy(self.parameters)
        paras["Axial Resilient Modulus"] = mr
        paras["sigma3"] = pressure
        paras["sigmaD"] = axial_cyclic_stress
        theta = axial_cyclic_stress + 3 * pressure
        paras["theta"] = theta
        tau = np.sqrt(2) / 3 * axial_cyclic_stress
        paras["tau"] = tau
        paras["theta / pa"] = theta / self.pa
        paras["tau / pa"] = tau / self.pa
        paras["tau / pa + 1"] = tau / self.pa + 1
        return paras

    def drop_row(self, x_, aim_rows, ):
        # drop row list
        row_list = []
        for seq_i in range(15):
            # aim range
            aim_range = range((seq_i+1)*100-aim_rows, (seq_i+1)*100)
            # aim average value
            x_ave = x_.iloc[aim_range].mean()
            for cycle_i in range(100):
                row_i = seq_i*100 + cycle_i
                # add row to drop list when the value is 
                # less than half of the aim
                # or more than twice of the aim
                if not 1/4*x_ave <= x_.iloc[row_i] <= 2*x_ave:
                    row_list.append(row_i)
        return row_list

    def regression(self, train_x, train_y, check_x, check_y):
        regr = linear_model.LinearRegression()
        # regression
        regr.fit(train_x, train_y)
        self.regr_result["coef_"] = regr.coef_
        self.regr_result["intercept_"] = regr.intercept_
        # determine and append k1
        # k1 = exp(intercept)/Pa
        k_temp.append(np.exp(regr.intercept_)/self.pa)
        # append k2 (k2 = coef_[0]) if it exists
        # append k3 (k3 = coef_[1])
        for k_2_3_i in regr.coef_:
            k_temp.append(k_2_3_i)
        # prediction
        predict_y = regr.predict(check_x)
        self.regr_result["MSE"] = mean_squared_error(check_y, predict_y)
        self.regr_result["R2"] = r2_score(check_y, predict_y)


    def main(self, ):
        # print hint
        print("Regressing...")
        count_i = 1
        for folder_i in self.file_names:
            # initiate start time
            self.set_time0()
            # log
            self.log["Folder"].append(folder_i)
            # generate mr file path and read file
            csv_mr_path = os.path.join(self.root_path, "{}_mr.csv".format(folder_i))
            temp1 = self.readcsv(csv_mr_path, ).copy()
            # get rows to be droppped
            self.rows2drop = self.drop_row(
                temp1["Axial Resilient Modulus"], 5, 
            )
            # calculating some necessary parameters
            try:
                # initiate input for calculation
                input_ = [
                    temp1["Axial Resilient Modulus"], 
                    temp1["Pressure"], temp1["Cyclic Axial Stress"],
                ]
                # temp2 = [theta, tau, theta / pa, sigmaD / pa, tau / pa + 1]
                temp2 = self.calculation(*input_)
                for model_i in self.model_xs:
                    # print hint
                    print("Processing[{}/{}]>> {} {}".format(count_i, self.folder_no, folder_i, model_i))
                    # generate reg file path
                    csv_reg_path = os.path.join(
                        self.root_path, "{}_{}_reg.csv".format(folder_i, model_i, ),
                    )
                    xs = self.model_xs[model_i]
                    # the number of x
                    x_no = len(xs)
                    # y train and predict dataframe
                    temp2["y"] = np.log10(temp2["Axial Resilient Modulus"]*1000)
                    y_train = temp2["y"].iloc[self.rows_for_train]
                    y_predict = temp2["y"].iloc[self.rows_for_predict]
                    # x train and predict empty dataframes
                    x_train = pd.DataFrame([])
                    x_predict = pd.DataFrame([])
                    # get x dataframes
                    for x_i in range(x_no):
                        temp2["x{}".format(x_i)] = np.log10(temp2[xs[x_i]])
                        x_train = pd.concat(
                            [x_train, temp2["x{}".format(x_i)].iloc[self.rows_for_train]], 
                            ignore_index=True, axis=1, 
                        )
                        x_predict = pd.concat(
                            [x_predict, temp2["x{}".format(x_i)].iloc[self.rows_for_predict]], 
                            ignore_index=True, axis=1, 
                        )
                    # drop rows
                    x_train = x_train.drop(self.rows2drop).to_numpy()
                    y_train = y_train.drop(self.rows2drop).to_numpy()
                    # do regression
                    self.regression(x_train, y_train, x_predict, y_predict, )
                    # append result
                    self.result_sum["Folder"].append(folder_i)
                    self.result_sum["Model"].append(model_i)
                    for item_i in self.regr_result:
                        self.result_sum[item_i].append(self.regr_result[item_i])
                    # output data used for regression 
                    self.writecsv(
                        pd.DataFrame(temp2), 
                        os.path.join(self.root_path, csv_reg_path, )
                    )
            except:
                # print hint
                print("!!! Regressing: {}".format(folder_i, ))
            # log duration
            self.log["Regression"].append(self.duration(self.start_time))
            count_i += 1
        self.writecsv(
            pd.DataFrame(self.result_sum), 
            os.path.join(self.root_path, "RLTT_reg_result.csv", )
        )
        # write log
        self.writecsv(
            pd.DataFrame(self.log), 
            os.path.join(self.root_path, "last_5_seq_regr_log.csv", )
        )


# generate file names
# file path
# Linux path = "/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"
# Windows path = r"D:\TestData\RLTT\Mr"
root_path = r"/mnt/c/Users/Chuanjun Lau/Documents/RLTT_Data/Mr"
# get all files
for i, j, k in os.walk(root_path):
    if len(j) != 0:
        temp0 = k
# get all "*mr.csv" files
file_names = []
for file_i in temp0:
    temp = file_i.split("_")
    if temp[-1].split(".")[0] == "mr":
        file_names.append(file_i.split("_mr")[0])

# Set particular folder name to analysis
# file_names = ["RLTT_C16_95_1_3_mr.csv"]

demo = DoRegression(root_path, file_names, )
demo.main()
