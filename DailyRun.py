import numpy as np
import os 
import pandas as pd
import glob
from tqdm import tqdm 
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import sys
sys.path.append(r'C:\Users\meypau\OneDrive - Fraunhofer\Desktop\04_Programming\00_Functions')
sys.path.append(r'C:\Users\meypau\OneDrive - Fraunhofer\Desktop\04_Programming\00_Functions\NacelleLidar')
sys.path.append(r'C:\Users\meypau\OneDrive - Fraunhofer\Desktop\04_Programming\00_Functions\SonicProcessing')
from f_data2dat_functions import repitition_check, Delta_check, Range_check, rolling_median_despiking,IncrementDistributionFilter,  rollingDespiking, RollingIncrement
from f_PSD import plot_PSD
from Basics import time_axis
import logging
import pyarrow as pa
import pyarrow.csv as csv
from f_read import get_fresh_data, read_reference_Data


plt.close('all')

height = '110m'
#height2 = '55m'


from ReadSonicData_V4 import Sonic
plot = True
bp = r'Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\\'
basepath = bp + r'Data\post_processed\Daily10min\\'
#rawpath  = bp + r'Data\UpgradeData\ASCII\\'
rawpath  = bp + r'Data\UpgradeData\ASCII_MeyPau\\'
rawpath  = bp + r'Data\UpgradeData_corrected\ASCII\\'
#%%
sd110 = Sonic('Gill110', hor_rotation = -121.31, height = '110m', basepath = basepath)#-121.31-90
sd55 = Sonic('Gill55',   hor_rotation = -121.31, height = '55m', basepath = basepath)
sd25 = Sonic('Thies_25', hor_rotation = -122.04, height = '25m', basepath = basepath)#+90)<
#%%
import datetime
idx = pd.date_range(end = pd.to_datetime("today").date(), freq = '12h', periods = 10)
idx = pd.date_range(end = pd.to_datetime("today").date(), freq = '12h', start = pd.to_datetime('2022-07-08'))

#%%

for sd in [sd25, sd55,sd110]:#sd25,
    for start, end in zip(idx[:-1],idx[1:]):

        st,en =start.strftime('%Y_%m_%d_%H%S'),end.strftime('%Y_%m_%d_%H%S')
        #
        
        
        sd.read_raw(rawpath, start = st, end = en, overlap = '1h')
        plt.close('all')
        if 'data_raw' in dir(sd) and len(sd.data_raw)>0:
            sd.filter_raw_data(plot_reference=False)
            sd.calc_parameters()
            # if plot:
            #     sd.plot_things()
            
            sd.savedata(typ = 'filtered')
    sd.savedata(typ = 'Daily_Stats')
      


