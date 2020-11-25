import os
import re
import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

# input_file = 'pssm_all/A0A0B4K7J2.pssm'
# output_file = 'pssm/A0A0B4K7J2.csv'

pssm = open(input_file).readlines()[3:-6]

# df = pd.DataFrame(columns=[i for i in range(0,20)])
# for line in pssm:
#     df.loc[len(df)] = [float(i) for i in line.split()[2:22]]

# df.to_csv(output_file,header=None,index=None,float_format='%.3f')

fout = open(output_file,"w")
for line in pssm:
    fout.write(f"{','.join([i for i in line.split()[22:42]])}\n")