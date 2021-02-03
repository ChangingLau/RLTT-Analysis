#!/bin/python3

# import modules

import pandas as pd 
import numpy as np 
import datetime as dt 

class Combine:
    # initial some variables
    def __init__(self, file_name, force, pressure, ):
        self.csv_name = file_name  # csv file names list
        self.overall_data_set = pd.DataFrame(
            [], 
            columns=[
                "Time", "Actuator Displacement", "Load", 
                "Pressure", "Axial LVDT1", "Axial LVDT2", 
            ]
        )  # dataset for containing all data
        self.data_set = None  # dataset for containing data of (seq1 -seq15)
        self.seqs = list(range(16))  # sequences list of RLTT
        self.target_load = pd.DataFrame(force, columns=["TargetLoad"])   # target forces dataframe of RLTT
        # target confining pressure of RLTT
        self.target_pressure = pd.DataFrame(pressure, columns=["TargetPressure"])  
        self.compensation_load = pd.DataFrame([], columns=["Compensation"])
    
    def date_time(self):
        date = str(dt.datetime.now()).split()[0].replace("-","")
        time = str(dt.datetime.now()).split()[1].split(".")[0].replace(":","")
        return "{}{}".format(date, time)

    def readcsv(self, ori, sep=",", ):
        return pd.read_csv(ori, sep=sep, )  # read csv and return it as dataframe 

    def writecsv(self, content, target, index=False, ):
        # write csv with dataframe
        content.to_csv("{}_{}{}".format(target[:-4], self.date_time(), target[-4:]), index=index, )  
    
    def read_csvs(self):
        # access csv files and combine them into a single dataframe named self.overall_data_set
        for file in self.csv_name:
            temp0 = self.readcsv(file).iloc[1:, 2:]
            self.overall_data_set = pd.concat(
                [self.overall_data_set, temp0, ], ignore_index=True, 
            ).astype('float64')
        # assign the copy of data of (seq1-seq15), and reset the index
        self.data_set = self.overall_data_set.drop(list(range(250000))).reset_index(drop=True)
        
    def calculation(self):
        # determine the confining compensation load
        self.compensation_load["Compensation"] = np.pi * (0.025 / 2) ** 2 * self.data_set["Pressure"].copy()
        # determine target load according to:
        #   target load = target load from spec + confining compensation load
        self.target_load["TargetLoad"] = \
            self.target_load["TargetLoad"].add(self.compensation_load["Compensation"], fill_value=0)
        # concatenate confining compensation load and target load to self.data_set
        self.data_set = pd.concat(
            [
                self.data_set, self.compensation_load["Compensation"], 
                self.target_load["TargetLoad"], self.target_pressure["TargetPressure"], 
            ], 
            axis=1, 
        )
        # define a function to determine the relative error
        def cal_error(A, B, ):
            # A = current value
            # B = target value
            # determine load and pressure error, and concat them to data set
            load_error = self.data_set[B].sub(self.data_set[A], fill_value=0)
            return pd.DataFrame(
                load_error.div(self.data_set[B], fill_value=0) * 100,
                columns=["Rela{}Error".format(A)], 
            )
        self.data_set = pd.concat(
            [self.data_set, cal_error("Load", "TargetLoad"), ], axis=1, 
        )
        self.data_set = pd.concat(
            [self.data_set, cal_error("Pressure", "TargetPressure"), ], axis=1, 
        )
           
    def main(self):
        self.read_csvs()
        self.writecsv(self.data_set, "RLTT_1-15_origin.csv")
        self.calculation()
        # generate ouput index
        single_load = list(range(100)) + list(range(100,500,50))
        single_index = []
        for load_i in list(range(100)):
            single_index += [load_i*500+load_j for load_j in single_load]
        # self.writecsv(pd.DataFrame(single_index), "single_index.csv")
        overall_index = []
        for seq_j in list(range(15)):
            overall_index += [seq_j*50000+seq_k for seq_k in single_index]
        # overall_index.pop(-1)
        # self.writecsv(pd.DataFrame(single_index), "overall_index.csv")
        # output data for comparison
        self.writecsv(self.data_set.loc[overall_index,:], "RLTT_1-15_comparison.csv")



# define a function to generate cycle load list for different sequences
def generate_cycle_load_list(axial_load, timeseries, ):
    cycle_load_list = []
    # determine the coefficient between time(series) and theta
    theta_time_coef = (360 / 0.1) / 180 * np.pi
    for timeseq_i in range(len(cycle_timeseries)):
        # determine the target force according to haversine pulse shape when it is in haversine pulse
        if timeseq_i < 50:
            theta = cycle_timeseries[timeseq_i] * theta_time_coef
            target_force = ((1 - np.cos(theta)) / 2 * 0.9 + 0.1) * axial_load
            cycle_load_list.append(target_force)
        # set target force as contact load when it is in rest period
        else:
            rest_seqs = len(cycle_timeseries) - timeseq_i
            target_force = list([axial_load * 0.1] * rest_seqs) 
            cycle_load_list += target_force
            break
    return cycle_load_list


# generate file names and append them into a list named file_list
csv_list = range(1, 22)
file_partial_name = "MF015- Unbound Material Resilient Modulus AASHTO T307_trial_"
file_list = ["{}{}.csv".format(file_partial_name, csv, ) for csv in csv_list]
# define target confining pressure (1-15 sequences)
pressure = [20.7, 20.7, 20.7, 34.5, 34.5, 
            34.5, 68.9, 68.9, 68.9, 103.4, 103.4, 
            103.4, 137.9, 137.9, 137.9, 
]
# define target axial stress (1-15 sequences)
axial_stress = [20.7, 41.4, 62.1, 34.5, 68.9, 
        103.4, 68.9, 137.9, 206.8, 68.9, 103.4, 
        206.8, 103.4, 137.9, 275.8, 
]
# area of load pad = pi*d^2/4
area = np.pi * 0.1 ** 2 / 4
# define target axial load (1-15 sequences)
axial_load = [stress_i * area for stress_i in axial_stress]
time_gap=0.002  # time gap between data accessing in each cycle
cycle_duration = 1  # duratin of whole of pulse
# timeseries of a cycle
cycle_timeseries = [round(timeseq_j*time_gap, 3) for timeseq_j in range(int(cycle_duration/time_gap))]
# generate load and pressure for each sequence (1-15 sequences)
seq_load_list = []
seq_pressure_list = []
for seq_i in range(15):
    seq_load_list += list(generate_cycle_load_list(axial_load[seq_i], cycle_timeseries, ) * 100)
    seq_pressure_list += list([pressure[seq_i]] * (len(cycle_timeseries) * 100))


# Start combination
demo = Combine(file_list, seq_load_list, seq_pressure_list, )
demo.main()
# print(demo.overall_data_set.shape)
# print(demo.data_set.shape)
# print(demo.target_load.shape)
# print(demo.compensation_load.shape)
