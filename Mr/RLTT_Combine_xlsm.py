# import libraries
from os.path import split
import pandas as pd
import os

data_set = {}
sheet_set = ["Seq{}_data".format(seq_i) for seq_i in range(1, 16)]
file_dir = r"C:\Users\Chuanjun\Downloads\RLTT_0516"
file_set = []

for root, dirs, files in os.walk(file_dir):
   for name in files:
       if name.endswith("xlsm"):
           file_set.append(name)

for file_i in file_set:
    file_path = os.path.join(file_dir, file_i)
    data_temp = []
    for sheet_i in sheet_set:
        data_temp.append(
            pd.read_excel(
                file_path, sheet_i, header=None, 
                skiprows=73, nrows=1,
            ).iloc[0, 7]
        )
    data_set[file_i.split(".")[0]] = data_temp

output_data = pd.DataFrame(
    data_set, 
    index=[sheet_j.split("_")[0] for sheet_j in sheet_set]
)

output_data.to_csv(
    os.path.join(file_dir, "output.csv")
)