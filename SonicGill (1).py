# -*- coding: utf-8 -*-
"""
Sonic 110m / 55m Gill
Original author: Dr. Pedro Santos

"""

#%% [Import module]
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import glob
#import pymysql
#from sqlalchemy import create_engine
import asyncio

#from datetime import datetime, timedelta, timezone
import datetime
from datetime import timedelta
from datetime import timezone


import os
import sys
module_path = os.path.abspath(r"C:\Users\hunlin\PythonConnector")
if module_path not in sys.path:
    sys.path.append(module_path)

from OneDasConnector import OneDasConnector

# nest_asyncio has bugs incompatible with Spyder 5.1.5 [2021.12.29]
# A quick fix: Downgrade to conda install spyder==4.2.5
import nest_asyncio
nest_asyncio.apply()
#%% [Sonic Offsets]
"""
Installation offsets in each Sonic anemometer. 
i.e. not necessary align with geographical North
"""
# offset for each sonic
# for the 55-m
offset55=121.84-90;
# for the 110-m
offset110=121.31-90;
#%% Input time of interest

infstart = datetime.datetime(2022, 3, 27, 0, 0, tzinfo=timezone.utc)
infend = datetime.datetime(2022, 3, 27, 0, 10, tzinfo=timezone.utc)

#%% [Sonics 110m raw data]

# settings and password
scheme = "https"
host = "onedas.iwes.fraunhofer.de"
port = 443
username = "ashim.giyanani@iwes.fraunhofer.de"
password = input("Please enter your password: ")

begin = infstart
end   = infend 
# must all be of the same sample rate
channel_paths = [
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_u/20 Hz",
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_v/20 Hz",
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_w/20 Hz",
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_SonicTempC/20 Hz",
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_SpeedOfSound/20 Hz"

]

connector = OneDasConnector(scheme, host, port, username, password) 

params = {
    "ChannelPaths": channel_paths,
}

data = asyncio.run(connector.load(begin, end, params))
sampleRate = 20 #Hz
time = [begin + timedelta(seconds=i/sampleRate) for i in range(len(data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_u/20 Hz"].values))]

filename_Sonic110 = "rawSonics110_"+begin.strftime('%Y%m%d%H')+"_to_"+end.strftime('%Y%m%d%H')

raw_Sonic110 = pd.DataFrame()
raw_Sonic110['u_usa'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_u/20 Hz"].values
raw_Sonic110['v_usa'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_v/20 Hz"].values
raw_Sonic110['w_usa'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_w/20 Hz"].values
raw_Sonic110['T'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_SonicTempC/20 Hz"].values 
raw_Sonic110['c_s'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION/gill_115_SpeedOfSound/20 Hz"].values


raw_Sonic110['time'] =time
raw_Sonic110['time']= pd.to_datetime(raw_Sonic110['time'], format=('%Y-%m-%d %H:%M:%S.%f'[:-2])) 

raw_Sonic110.set_index(['time'], inplace=True)
raw_Sonic110.index = raw_Sonic110.index.tz_convert(None)
raw_Sonic110['T'] = raw_Sonic110['T']+273.15
del scheme, host, port, username, password, params

#%% [function for calculation cross-variance ]
"""
Unbiased estimator of the cross-variance between A and B
Computing the fluctuations from a linear trend +
It accounts for the sample size on the calculation
- To check: scipy detrend may be different from MATLAB detrend
"""
from scipy import signal
def cross_variance_linear(a,b):
    c = np.mean(signal.detrend(a,type='linear')*signal.detrend(b,type='linear'))\
        *1/(1-1/len(a[~np.isnan(a)]))
    return c

#%% [flux]
VonKarman_const = 0.4; g =9.81
def gill_2r_fluxes_sample(s, offset=offset110, z_sonic=110.):
    col_name = ['U_horz', 'U_vec', 'wind_direction', 'inflow_angle',\
            'U', 'V', 'W', 'T ',\
            'u_max', 'v_max', 'w_max', 'T_max', \
            'u_min', 'v_min', 'w_min', 'T_min', \
            'cov_uu ', 'cov_uv ', 'cov_uw ',\
            'cov_vv', 'cov_vw', 'cov_ww',\
            'cov_uT ', 'cov_vT ', 'cov_wT ','cov_TT',\
            'U_horz_std', 'U_vec_std ',\
            'u_star', '1/L', 'zL']
    if s['u_usa'].isna().sum() < tlength*sampleRate:
        # UCAR ISFS
        Vaz = offset + 90
        #wd = (np.rad2deg(np.mean(np.arctan2(-s['v_usa'],s['u_usa']))) + Vaz)%360
        si = s['v_usa']/np.hypot(s['u_usa'],s['v_usa'])
        co = s['u_usa']/np.hypot(s['u_usa'],s['v_usa'])  
        wd_prime = ((180+Vaz-np.rad2deg(np.arctan2(np.sum(si), np.sum(co)))))%360
        wd = (wd_prime+180)%360
        # Wind vector Means 
        uh = np.mean(np.hypot(s['u_usa'],s['v_usa']))
        uvec = np.mean(np.sqrt(s['u_usa']**2+s['v_usa']**2+s['w_usa']**2))
        tilt = np.mean(np.arctan2(s['w_usa'], np.hypot(s['u_usa'],s['v_usa'])))*180/np.pi #inflow angle!
        # Wind vector STDs
        std_uh = np.std(np.hypot(s['u_usa'],s['v_usa']))
        std_uvec = np.std(np.sqrt(s['u_usa']**2+s['v_usa']**2+s['w_usa']**2)); 
        std_tilt = np.std(np.arctan2(s['w_usa'], np.hypot(s['u_usa'],s['v_usa'])))*180/np.pi #inflow angle!
        std_windvec = pd.DataFrame([std_uh,std_uvec]).transpose()
        std_windvec.columns =['U_horz_std', 'U_vec_std']
        
        # Output wind vector characteristics:
        df = pd.DataFrame([uh,uvec,wd,tilt]).transpose()
        df.columns = ['U_horz', 'U_vec', 'wind_direction', 'inflow_angle']
        
        # ----coordinate rotation---- #:
        alpha = np.arctan2(np.mean(s['v_usa']),np.mean(s['u_usa']));
        R01 = np.array([[np.cos(alpha), np.sin(alpha), 0], [-np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]])
        U1 = R01.dot(np.array([s['u_usa'], s['v_usa'], s['w_usa']])) # 1st rotation
        beta = np.arctan2(np.mean(U1[2,:]), np.sqrt(np.mean(U1[0,:])**2 + np.mean(U1[1,:])**2))
        R12 = np.array([[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]])
        U2 = R12.dot(U1); # 2nd Rotation
        U3 = U2.transpose()
        
        # ---- Rotated wind components ----#
        U3 = pd.DataFrame(U3)
        U3.columns = ['U', 'V', 'W']
        mean_U3 =pd.DataFrame(U3.mean()).transpose()
        mean_U3['T'] = np.mean(s['T'])
        df = pd.concat([df, mean_U3], axis=1)
        df['u_max'] = max(U3['U'])
        df['v_max'] = max(U3['V'])
        df['w_max'] = max(U3['W'])
        df['T_max'] = max(s['T'])
        df['u_min'] = min(U3['U'])
        df['v_min'] = min(U3['V'])
        df['w_min'] = min(U3['W'])
        df['T_min'] = min(s['T'])
        
        # ---- Reynolds stress tensor in the 2R coordinate system ----#
        df['cov_uu'] = cross_variance_linear(U3['U'],U3['U'])
        df['cov_uv'] = cross_variance_linear(U3['U'],U3['V'])
        df['cov_uw'] = cross_variance_linear(U3['U'],U3['W'])
        df['cov_vv'] = cross_variance_linear(U3['V'],U3['V'])
        df['cov_vw'] = cross_variance_linear(U3['V'],U3['W'])
        df['cov_ww'] = cross_variance_linear(U3['W'],U3['W'])
        
        # ---- Heat fluxes in the 2R coordinate system ----#
        df['cov_uT'] = cross_variance_linear(U3['U'],s['T'])
        df['cov_vT'] = cross_variance_linear(U3['V'],s['T'])
        df['cov_wT'] = cross_variance_linear(U3['W'],s['T'])
        df['cov_TT'] = cross_variance_linear(s['T'],s['T'])
        
        df = pd.concat([df, std_windvec], axis=1)
        # friction velocity u*
        df['u_star'] = np.power(df['cov_uw']**2+df['cov_vw']**2, 0.25)
        # Inverse of Obukhov length (1/L)
        df['1/L'] =VonKarman_const*(g/df['T'])*df['cov_wT']/(-df['u_star']**3)
        # zL considering the height of the sonic
        df['zL'] =z_sonic*df['1/L'] 
        df.columns=col_name
        df['time'] = s.index[0]
        df.set_index(['time'], inplace=True)
    else:
        df = pd.DataFrame([np.nan] * 31)
        df = df.transpose()
        df.columns=col_name
        df['time'] = s.index[0]
        df.set_index(['time'], inplace=True)
    return df


#%% [De-spike and fill vancancies]

# X-sec ensemble means
tlength = 600; #600 #[sec]
threshold_accept= 0.1;
timestamp = pd.date_range(start=pd.Timestamp(begin).tz_convert(None), end=pd.Timestamp(end).tz_convert(None), freq=f"{tlength}S")

quality_control = False; quantity_control = True
filt_sonic = []
if quality_control:
    for k in np.arange(len(timestamp)-1):
        s = raw_Sonic110.loc[(raw_Sonic110.index >= timestamp[k]) & (raw_Sonic110.index <timestamp[k+1] )].copy()
        
        if s['u_usa'].isna().sum()/(sampleRate*tlength) < threshold_accept :
            #print(f'accepted {timestamp[k]}')
            n_window = int(tlength*sampleRate)
            median = s.rolling(window=n_window,axis=0,min_periods=1).median()
            std = s.rolling(window=n_window,axis=0,min_periods=1).std()
            outliers = (s - median).abs() > 3*std
            s = s.where(~outliers, other=np.nan)
            s = s.fillna(median)
            if s.isna().sum().sum() >0:
                s = s.interpolate(method='linear', axis=0)
                s = s.fillna(method='ffill', axis=0).fillna(method='bfill', axis=0)
        else:
            avail = round(s['u_usa'].isna().sum()*100/(sampleRate*tlength), 2)
            print(f'{avail}% not pass availability test {timestamp[k]}')
            s.iloc[:] = np.nan
        filt_sonic.append(s)
    
elif quantity_control:
    for k in np.arange(len(timestamp)-1):
        s = raw_Sonic110.loc[(raw_Sonic110.index >= timestamp[k]) & (raw_Sonic110.index <timestamp[k+1] )].copy()    
        if len(s) > (sampleRate*tlength)*0.99 :
            n_window = int(tlength*sampleRate)
            median = s.rolling(window=n_window,axis=0,min_periods=1).median()
            std = s.rolling(window=n_window,axis=0,min_periods=1).std()
            outliers = (s - median).abs() > 3*std
            s = s.where(~outliers, other=np.nan)
            s = s.fillna(median)
            if s.isna().sum().sum() >0:
                s = s.interpolate(method='linear', axis=0)
                s = s.fillna(method='ffill', axis=0).fillna(method='bfill', axis=0)
        else:
            print(f'not pass availability test {timestamp[k]}')
            s.iloc[:] = np.nan
        filt_sonic.append(s)    


P = []
for j in np.arange(len(timestamp)-1):
    # Change the height and offset for 55/110m Gill sonic
    process_s = gill_2r_fluxes_sample(filt_sonic[j], offset=offset110, z_sonic=110.)
    #process_s = gill_2r_fluxes_sample(filt_sonic[j], offset=offset55, z_sonic=55.)
    print(j)
    P.append(process_s)

avg_Sonics = pd.concat(P, axis=0)



#%% [Compare results]
'''Comapre results with MATLAB processed sonic'''

# settings and password
scheme = "https"
host = "onedas.iwes.fraunhofer.de"
port = 443
username = "ashim.giyanani@iwes.fraunhofer.de"
password = input("Please enter your password: ")

begin = infstart
end = infend


# must all be of the same sample rate
channel_paths = [
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION_PROCESSED/gill_110m_U_horz/600 s",
    "/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION_PROCESSED/gill_110m_wind_direction/600 s"
 
]

# load data
connector = OneDasConnector(scheme, host, port, username, password) 
# without authentication: connector = OneDasConnector(scheme, host, port)

params = {
    "ChannelPaths": channel_paths,
}

data = asyncio.run(connector.load(begin, end, params))
sampleRate = 1/600 # 1 Hz (adapt to your needs)
time = [begin + timedelta(seconds=i/sampleRate) for i in range(len(data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION_PROCESSED/gill_110m_U_horz/600 s"].values))]


sonic100_MATLAB = pd.DataFrame()
sonic100_MATLAB['U_horz'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION_PROCESSED/gill_110m_U_horz/600 s"].values
sonic100_MATLAB['wind_direction'] = data["/AIRPORT/AD8_PROTOTYPE/METMAST_EXTENSION_PROCESSED/gill_110m_wind_direction/600 s"].values

sonic100_MATLAB['time'] =time
sonic100_MATLAB['time']= pd.to_datetime(sonic100_MATLAB['time'], format=('%Y-%m-%d %H:%M:%S.%f'[:-2])) 

sonic100_MATLAB.set_index(['time'], inplace=True)
sonic100_MATLAB.index = sonic100_MATLAB.index.tz_convert(None)



fig = plt.figure() 
ax = fig.add_subplot(111)
ax.plot(sonic100_MATLAB['wind_direction'], label='MATLAB')
ax.plot(avg_Sonics['wind_direction'], label='python')
ax.legend()


fig = plt.figure() 
ax = fig.add_subplot(111)
ax.plot(sonic100_MATLAB['U_horz'], label='MATLAB')
ax.plot(avg_Sonics['U_horz'], label='python')
ax.legend()


