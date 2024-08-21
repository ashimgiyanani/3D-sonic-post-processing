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
from f_data2dat_functions import repitition_check, Delta_check, Range_check, rolling_median_despiking,IncrementDistributionFilter,  rollingDespiking, RollingIncrement, find_outlier, find_mean_outlier, filter_thies_outlier_byTemperature, manual_filter_sonic
from f_PSD import plot_PSD
from Basics import time_axis
import logging
import pyarrow as pa
import pyarrow.csv as csv
from f_read import get_fresh_data, read_reference_Data
from scipy import signal
def cross_variance_linear(a,b):
    ## From Lin-Ya Hung / Pedro Santos, Fraunhofer IWES
    c = np.mean(signal.detrend(a,type='linear')*signal.detrend(b,type='linear'))\
        *1/(1-1/len(a[~np.isnan(a)]))
    return c

def cross_variance_linear_DataFrame(df,c1,c2):
    ## From Lin-Ya Hung / Pedro Santos, Fraunhofer IWES
    a = df[c1]
    b = df[c2]
    c = np.mean(signal.detrend(a,type='linear')*signal.detrend(b,type='linear'))\
        *1/(1-1/len(a[~np.isnan(a)]))
    return c

def hex2int(val):
    if val is val: ## not NAN
        if isinstance(val,str):
            return int(val.replace(';',''),16)
        elif isinstance(val,float) or isinstance(val,int):
            return val
    else:
        return np.nan

#%%

def met2hor(u,v):
    u_hor = np.sqrt(u**2+v**2)
    wd = (180+np.rad2deg(np.arctan2(u,v)))%360
    return u_hor,wd

def insert_AVG_STD(columns, string = 'AVG', sign = '['):
    retcols = []
    for c in columns:
        split = c.split(sign)
        retcols.append(split[0] + string +sign+ split [1])
    return retcols

class Sonic:
    #def __init__(self):
    
    def __init__(self,model = 'Gill115', height=None, hor_rotation = 0 , ver_rotation = 0,basepath = None):
        self.basepath = basepath
        self.height = height
        self.model = model
        self.hor_rotation = hor_rotation
        self.ver_rotation = ver_rotation
        if 'gill' in self.model.lower():
            self.device = 'Gill'
        elif 'thies' in self.model.lower():
            self.device = 'thies'
        else:
            raise ValueError('Unknown Device!')
        
        if basepath is not None:
            self._initDirectory()
        self._setupLogFile()
        
    def _setupLogFile(self):
        if self.basepath is not None:
            self.logfile = self.basepath + f'Processing_{self.height}.log'
            formatter = logging.Formatter('%(asctime)s -[%(levelname)s] %(funcName)s - %(message)s')
            handler = logging.FileHandler(self.logfile)    
            handler.setFormatter(formatter)
            logger = logging.getLogger(self.model)
            logger.setLevel(logging.DEBUG)
            logger.addHandler(handler)
            self.logger = logger
            ## from https://stackoverflow.com/questions/11232230/logging-to-two-files-with-different-settings
            
        else:
            self.logfile = None
            self.logger = logging.getLogger(self.logfile)
       
        
        
    
        
            
            
        # logging.basicConfig(filename=self.logfile,  filemode='w',
        #                     format='%(asctime)s -[%(levelname)s] %(funcName)s - %(message)s', datefmt='%Y-%m %d %H:%M:%S')
        # self.logger = logging.getLogger(self.logfile)
        # self.logger.setLevel(logging.DEBUG)        
        self.logger.info(f'##################\nCreated Log file for {self.model} ')
        
        
    def _initDirectory(self):
        print('Initializing Directory..')
        self.dpath = self.basepath + 'QC_Data\\' + self.height + '\\'
        self.figpath = self.basepath + 'Figures\\' + self.height + '\\'        
        
        for path in [self.dpath, self.figpath]:
            if not os.path.exists(path):
                print(f'{path} did not exist, will create one')
                os.mkdir(path)

        
        
    def read_raw(self, path, start = '2000_01_01_0000', end = '2100_01_01_0000',  overlap = '0min', columns = dict()):
        self.logger.info('___________________________')
        self.logger.info('Reading Raw data...')
        self.start = pd.to_datetime(start, format = '%Y_%m_%d_%H%S')
        self.end = pd.to_datetime(end, format = '%Y_%m_%d_%H%S')
        if overlap is None:
            overlap = '0min'
        
            
        self.read_offset = pd.to_timedelta(overlap)
        self.start_read = self.start- self.read_offset
        self.end_read = self.end+self.read_offset
        self.rename_columns = columns
        
        self.times = self.start_read.strftime('%Y%m%d%H%S') + '_'+ self.end_read.strftime('%Y%m%d%H%S') 
        self.logger.info(f'Reading Raw data...: {self.times}')
        #print(f'Reading from path: \n {path}')    
        #self.path = path
        
        self._getFileNames(path)
        self._readData()
        
            
        return self

    def from_raw(self,Sonic,  columns = None):
        print('----')
        print('Reading from Raw files')
        self.logger.info('Got Raw data from other height...')
        start,end = Sonic.start, Sonic.end
        self.read_offset = Sonic.read_offset
        #start,end = df.first_valid_index(), df.last_valid_index()
        df = Sonic.data_raw
        self.start, self.end = start, end
        self.times = start.strftime('%Y%m%d%H%S') + '_'+ end.strftime('%Y%m%d%H%S') 
        self.df_all_file = df
        
        self.data_raw = df.filter(like = self.device.lower()+'_'+self.height.replace('m',''))
    
    def _readData(self):
        dtypes = dict(TIMESTAMP = str, 
                      RECORD = object, 
                      gill_115_id = object,
                      gill_115_u = float, 
                      gill_115_v = float, 
                      gill_115_w = float, 
                      gill_115_unit = str, 
                      gill_115_SpeedOfSound = float, 
                      gill_115_SonicTempC = float, 
                      gill_115_status = float ,
                      gill_55_id = object,
                      gill_55_u = float, 
                      gill_55_v = float, 
                      gill_55_w = float, 
                      gill_55_unit = str, 
                      gill_55_SpeedOfSound = float, 
                      gill_55_SonicTempC = float, 
                      gill_55_status = float )
        kwargs = dict( skiprows = [0,2,3], na_values = ['NAN'], dtype = dtypes)
        if len(self.files)>0:
            

            df_all = pd.concat([pd.read_csv(file, **kwargs ) for file  in tqdm(self.files, desc = 'Reading data files: ')], axis = 0)
            print('Setting Timestamp..')
            try:
                df_all.index = pd.to_datetime(df_all['TIMESTAMP'], format = '%Y-%m-%d %H:%M:%S.%f')
            except:
                
                self.logger.debug('Format did not fit! Check')
                df_all.index = pd.to_datetime(df_all['TIMESTAMP'], errors='coerce')
                df_all = df_all.loc[df_all.index.dropna()]
        else:
            self.logger.debug('No files available to read...')
            df_all = pd.DataFrame()
            
        
        self.data_raw = df_all
        print('Finished Reading...')
        
    def _getFileNames(self,path):
        start, end = self.start_read.strftime('%Y_%m_%d_%H%M'), self.end_read.strftime('%Y_%m_%d_%H%M')
        print('Getting Filenames')
        try:
            files = self.all_files
        except:
            files = os.listdir(path)
            self.all_files = files
        self.files = [path + f for f in files if self.device.lower() in f and f[-19:-4] >= start and f[-19:-4] < end ]
        if len(self.files)>0:
            print('Found Files')
        else:
            print('Didn´t find any file!')
            
    def _correctColumnShift(self, plot = False):
        if 'Thies' in self.model:
            
            df = self.data_raw
            
            
            df['thies_ThiesStatus'] = df['thies_ThiesStatus'].apply(lambda val: hex2int(val))
            
            basecondition = ((df['thies_CheckSum'].isna()) & (df['thies_ThiesStatus'] >3) & (df["thies_AvTc"].between(-0.1,0.1)))
            
            if plot:
                figs, axe = plt.subplots(num = 3, clear= True)
                axe.plot(df['thies_Vy'], c = 'grey', alpha = 0.4, label = 'original')
                
                
                
            
            all_condition = np.zeros_like(df.iloc[:,0], dtype = bool)
            for i in range(10):
                
                c1 = ((df['thies_Vy']-df['thies_Vz'].shift()).abs() < df['thies_Vy'].diff().abs())# (~df["thies_AvTc"].shift(1).between(-0.1,0.1)))
                condition = basecondition & c1
            
                if condition.sum()>0:
                    
                    all_condition = all_condition | condition
            
                    self.logger.debug(f'Had to move {condition.sum()} columns!')
                    cols= [c for c in df.columns if not any(c in k for k in ["TIMESTAMP","RECORD"])]
                    
                    if plot:
                        axe.scatter(df.loc[condition].index, df.loc[condition, 'thies_Vy'], label = '')
                        
                        
                    
                    
                    df.loc[condition, cols[1:]] = df.loc[condition, cols[:-1]].values
                    df.loc[condition, cols[0]] = np.nan
                    
                    if plot:
                        axe.plot(df['thies_Vy'], c = 'grey', alpha = 0.4, label = 'Shifted')
                        axe.legend()
                    
                    
                
            #df[all_condition] = np.nan ## 
            self.logger.error(f'Had to remove {all_condition.sum()} datapoints, cannot reconstruct moving!')
            self.data_raw = df
            
            return all_condition.sum(), all_condition
            
        else:
            return 0, (self.data_raw.iloc[:,0].copy()*0).astype(bool)
        
    def filter_status(self, data):
        ## Filter for Status == 0 for Gill
        ## Filter exact bit values for Thies, see manuals:
            # https://www.thiesclima.com/db/dnl/4.383x.xx.xxx_US-Anemometer-3D_e.pdf   (Section 7.4.6.2)
            # http://gillinstruments.com/data/WindMaster/Manuals/1561-PS-0001%20WindMaster%20Windmaster%20Pro%20Manual%20Issue%2016.pdf (Section 11.4)
            # https://www.binary-code.org/bits/8?page=1
        status  = data['Status']
        if 'Thies' in self.model:
            statfill = status.fillna(255) ## will be filtered later, as 1111111
            status_binary = statfill.apply(lambda x: format(int(x), '08b')) ## from https://stackoverflow.com/questions/45018601/decimal-to-binary-in-dataframe-pandas
            
            ## eg: pd.Series([0,2,5,25]).apply(lambda x: format(int(x), '08b')).str.match('0....0..')
            # checkSum= data['CheckSum']
            strangeValues = [50,51,56,57,58,59] ## Status values, that cause strange measurements
            matchstr ='0.......'# '0....0..'
            condition = ((status_binary.str.match(matchstr)) & (~status.isin(strangeValues)) & (status.notna()))#& checkSum.notna() ## bit 0 and bit 5 indicate static and general malfunction
            ##~status_binary.str.match('....1...') 
            #condition = status_binary.str.match('0.......') #& checkSum.notna() ## bit 0 and bit 5 indicate static and general malfunction

        if 'Gill' in self.model:
            condition =  (status == 0)
        return condition
            
            
    def _getSignsNames(self):
        print('Cordinate System:')
        if 'Gill' in self.model:
            print('\t Gill Coordinate System...')
            signs = np.array([1,1,1,1]) ## mirror V component
            if self.height == '55m':
                rename_dict = dict(gill_55_u = 'u' , gill_55_v  = 'v', gill_55_w = 'w', gill_55_SonicTempC = 'T', gill_55_status = 'Status')
            elif self.height == '110m':
                rename_dict =  dict(gill_115_u= 'u' , gill_115_v = 'v', gill_115_w= 'w', gill_115_SonicTempC = 'T', gill_115_status = 'Status')
            else:
                raise ValueError('Unknown Height to get signs?')

        elif 'Thies' in self.model:
            print('\t Using Thies Coordinate System')
            signs = np.array([1,-1,1,1])
            rename_dict = dict(thies_Vx= 'v' , thies_Vy = 'u', thies_Vz= 'w', thies_AvTc = 'T', thies_ThiesStatus = 'Status', thies_CheckSum = 'CheckSum')
        else:
            raise ValueError('Unknown coordinate system')
        return signs, rename_dict
    
    
    def _makeIndexUnique(self, df):
        self.logger.debug('Index not unique, will remove first appearance')
        return df[~df.index.duplicated(keep = 'last')]
        
        
    def _makeIndexMonotonic(self, df):
        self.logger.debug('Index not monotonic increasing, will sort by datetime index')
        return df.sort_index()
        

    def filter_raw_data(self, nstds = [6,8],remove_faulty= 'all', plot_reference = True, reference = None, delta_check = True, 
                        Ranges = dict(u = [-45,45],v = [-45,45], w = [-6,10], T = [-20,60]), Deltas =dict(u = 10, v = 10, w = 10, T = 4),
                        plot = True, plot_incr_dist = False, plot_shift = False):
        ## from https://stackoverflow.com/questions/46964363/filtering-out-outliers-in-pandas-dataframe-with-rolling-median
        print('Filtering Raw data...')
        
        ## create df with rolling stats, upper and lower bounds
        self.logger.info('Filtering Raw Data...')
        
        
        try:
            self.filteredby
        except:
            self.filteredby = pd.DataFrame()
        
        
        signs, columns = self._getSignsNames()
        
        
        L_shift, shift_condition = self._correctColumnShift(plot = plot_shift)
        df = self.data_raw.copy()
        print(f'Time: {df.index[0]} -- {df.index[-1]}')

        
        
        
        if plot:
            fig, axes = plt.subplots(nrows = 4, sharex = True, figsize = (18,14), num = 1, clear = True)     
            axes[0].set(title =f'Period: {self.start.strftime("%Y-%m-%d_%H:%S")} -- {self.end.strftime("%Y-%m-%d_%H:%S")}')
        else:
            axes = [None,None,None,None]
        
        if plot_reference and axes [0] is not None:
            s,e = df.first_valid_index().tz_localize('UTC'), df.last_valid_index().tz_localize('UTC')
            ref = reference.loc[s:e]

            axes[0].plot(-ref['WS']*np.sin(np.deg2rad(ref['WD']+self.hor_rotation-90)), zorder = 20, label = 'Cup')
            axes[1].plot(-ref['WS']*np.cos(np.deg2rad(ref['WD']+self.hor_rotation-90)), zorder = 20, label = 'Cup')
            axes[2].plot(-ref['Z'], zorder = 20, label = 'Lidar')
            
        
        df_all_stat = df.rename(columns = columns)

        keep_condition = self.filter_status(df_all_stat)
        df_all_stat = df_all_stat[['u','v','w','T']]*signs
        df_qc = df_all_stat.copy()
        df_qc[~keep_condition] = np.nan
        
        filt_str = f'Filtered: All Status: {len(df_all_stat[~keep_condition])}, '
        
        self.filteredby.loc[df_qc.index[0],'ColumnShift'] = L_shift
        self.filteredby.loc[df_qc.index[0],'Status'] = len(df_all_stat[~keep_condition])
        
        if not df_qc.index.is_unique:
            df_qc = self._makeIndexUnique(df_qc)
        
        if not df_qc.index.is_monotonic:
            df_qc = self._makeIndexMonotonic(df_qc)
        
        
        
        
        if len(df_qc)>0:
            print('Filtering Columns:')
            manpath = r'Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData_corrected\SyncFiles\ManualRemovePeriods.xlsx'
            manuallist = pd.read_excel(manpath, sheet_name = 'Remove')
            manuallist = manuallist[manuallist['Height']== self.height]
            if len(manuallist)>0:
                manualcondition = manual_filter_sonic(df_qc, manuallist)
                df_qc[manualcondition] = np.nan
                filt_str += f' Manual Filter: {manualcondition.sum()}, '
                self.filteredby.loc[df_qc.index[0],'ManualFilter'] = manualcondition.sum()
            
            
            if 'thies' in self.model.lower():
                shiftcondition = filter_thies_outlier_byTemperature(df_qc['T'], ax = axes[3])
                df_qc[shiftcondition] = np.nan
                filt_str += f' Outlier_Temperatur: {shiftcondition.sum()}, '
                self.filteredby.loc[df_qc.index[0],f'Outlier_Temperatur'] = shiftcondition.sum()
            
            
            
            for i,(c, ax) in enumerate(zip(['u','v','w','T'], axes)):
                print(f'\t{c}...')
                filt_str += f'| {c}:'
                if ax is not None:
                    ax.plot( df_all_stat[c].dropna(), c = 'grey', alpha = 0.3)
                    if L_shift >0:
                        ax.scatter(df_all_stat[shift_condition].index, df_all_stat.loc[shift_condition, c],marker = '+',label = f'Shifted: {L_shift}')## next color
                    ax.scatter(df_all_stat[~keep_condition].index, df_all_stat.loc[~keep_condition, c],label = f'Status: {len(df_all_stat[~keep_condition])}')## next color
                    
                #if c in ['w','T']:
                    # try:
                # df_qc[c], L_outlier = find_mean_outlier(df_qc[c], ax = ax)
                # filt_str += f' Outlier: {L_outlier},'
                # self.filteredby.loc[df_qc.index[0], f' {c} Outlier'] = L_outlier                    
                
                for j,nstd in enumerate(nstds):
                    df_qc[c],L_roll = rolling_median_despiking(df_qc[c],  nstd = nstd, ax = ax)
                    filt_str += f' RollingMedian{j}: {L_roll}, '
                    self.filteredby.loc[df_qc.index[0],f'{c}_RollingMedian{j}'] = L_roll
                    


                for incr in [1,2]:
                    df_qc[c], L_increment = RollingIncrement(df_qc[c], ax = ax, increment = incr)
                    filt_str += f' Increment{incr}: {L_increment}, '    
                    self.filteredby.loc[df_qc.index[0], f' {c}_Increment{incr}'] = L_increment

                        
                    
                    
                df_qc[c], L_inc2 = IncrementDistributionFilter(df_qc[c], logger = self.logger, ax = ax, plotown = plot_incr_dist)
                self.filteredby.loc[df_qc.index[0],f'{c}_Increment_dist'] = L_inc2
                filt_str += f' Increment_Dist: {L_inc2}, '
                
                df_qc[c],L_rep = repitition_check(df_qc[c], c= 10, ax = ax)
                df_qc[c],L_ran = Range_check(df_qc[c], Ranges[c], ax = ax)
                
                self.filteredby.loc[df_qc.index[0],f'{c}_Repitition'] = L_rep
                self.filteredby.loc[df_qc.index[0],f'{c}_Range'] = L_ran
                
                
                if delta_check:
                    df_qc[c],L_delta = Delta_check(df_qc[c], dmax = Deltas[c], ax = ax)
                    filt_str += f' Delta: {L_delta}'
                    self.filteredby.loc[df_qc.index[0],f'{c}_Delta'] = L_delta
                    
                filt_str += f' Repitition: {L_rep}, Range: {L_ran}, kept: {len(df_qc[c].dropna())}'
                
                if ax is not None:
                    
                    try:
                        ax.set(ylim = [df_qc[c].min()-1, df_qc[c].max()+1], ylabel = c)
                    except:
                        print('Nothing given here')
                
            
                
            
            df_qc = df_qc[(df_qc.index >=self.start)&( df_qc.index < self.end)] ## remove offset
            
            
            if remove_faulty == 'all':
                df_qc[df_qc.isna().any(axis = 1)] = np.nan
                
            for i,(c, ax) in enumerate(zip(['u','v','w','T'], axes)):
                if ax is not None:
                    ax.plot(df_qc[c].dropna(), label = f'{c} Final: {df_qc[c].notna().sum()}', c = 'k')
                    ax.legend(loc = 'lower left')
                
                
            df_qc['T'] += 273.15 ## degree to Kelvin
            if ax is not None:
                time_axis(ax)
                self.qcfig = fig
            
            print(filt_str)
            self.logger.info(filt_str)
            self.data = df_qc.copy()
            df_qc = ''
            print('Finished Filtering...')
        else:
            self.logger.error(f'Not enough data for period {self.times}')
            
            
       
    def calc_parameters(self, avgtime= '10min', **kwargs):
        if len(self.data)>0:
            self._calcMeteorologicalComponents()
            self.add_AVG_Streamwise(avgtime= avgtime)
            self.calculateStatistics(avgtime = avgtime, **kwargs)
        else:
            self.logger.debug('No data available!')
        
    def add_AVG_Streamwise(self, avgtime = '10min'):
        # if not self.data.index.is_unique:
        #     self.logger.error(__name__+' Index is not unique!')
        
        resmean = self.data[['u_met','v_met']].resample(avgtime).mean()
        uvag, wdavg = met2hor(resmean['u_met'], resmean['v_met'])
        
        try:
            self.data['10min WD [deg]'] = wdavg.copy().reindex(self.data.index, method = 'ffill')
        except:
            self.logger.error(__name__+' Index is not unique!')
            self.data['10min WD [deg]'] =wdavg.sort_index().copy().reindex(self.data.sort_index().index, method = 'ffill')
        
        rotvec = R.from_euler('z', (180-self.data['10min WD [deg]'])%360, degrees = True)
        self.data[['10minAVG streamwise [m/s]','10minAVG crosswise [m/s]', '10minAVG vertical [m/s]']] =  rotvec.apply(self.data[['v_met', 'u_met' , 'w_met']])


    def check_index(self, method, data):
        if not data.index.is_unique:
            self.logger.error(method+' Index is not unique!')
        
        
        
    def _calcMeteorologicalComponents(self):
        print('Calculating Meteorological Coordinate System')
        rotvec = R.from_euler('z',[self.hor_rotation-90], degrees = True) # z angle defined positive anticlockwise! -90 offset, so U points towards east
        self.rotvec = rotvec
        metcoords = rotvec.apply(self.data[['u','v','w']])
        self.data[['u_met','v_met','w_met']] = metcoords
        self.data['U_hor'], self.data['WD']    = met2hor(self.data['u_met'], self.data['v_met'])  #( np.rad2deg(np.arctan2(self.data['u_met'],self.data['v_met']))+180 )%360
        self.data.drop(['u','v','w'],axis = 1,inplace = True)
        #self.rotated_raw = self.data.copy()



    
    
    def calculateStatistics(self, avgtime= '10min', min_count_percent = 40, dt = 1/20, inclusive_variances = False ):
        columns = ['T','u_met','v_met','w_met','U_hor',#'WD',
                   '10minAVG streamwise [m/s]','10minAVG crosswise [m/s]','10minAVG vertical [m/s]']
        
        newcolumns = ['T [K]','u [m/s]','v [m/s]','w [m/s]','U_hor [m/s]',#'WD [deg]',
                      'Streamwise [m/s]','Crosswise [m/s]','Vertical [m/s]']
        
                  
        data = self.data[columns].copy()
        data.columns = newcolumns
        res = data.dropna().resample(avgtime, label = 'left', closed = 'left')
        avg = res.mean()
        std = res.std()
        std.columns = insert_AVG_STD(std.columns, string = 'STD ')#['['.join(c.split['[']   c+' STD' for c in std.columns]
        avg.columns = insert_AVG_STD(avg.columns, string = 'AVG ')# [c+ ' AVG' for c in avg.columns]
        
        dt_avg = pd.to_timedelta(avgtime).total_seconds()
        data10 = pd.concat([avg, std], axis=1)
        #self.data10 = data10
        data10['U_hor Vec AVG [m/s]'],data10['WD Vec AVG [deg]'] = met2hor(data10['u AVG [m/s]'], data10['v AVG [m/s]'])
        ## Only consider periods, where more than mincount data
        data10['Datapoints'] = res.count().iloc[:,0].astype(int)
        data10[data10['Datapoints'] < min_count_percent/100 * dt_avg / dt] = np.nan
        #data10 = data10.loc[]
        if inclusive_variances:
            data10 = self._calculateVariances(data10, data)
        try:
            
            self.data10 = pd.concat([self.data10, data10], axis = 0)
            print('Adding to existing 10min')
        except AttributeError:
            try:
                self.data10
                print('I have it')
                self.logger.error('Concatenating Statistic Data did not work. Maybe first time?')
            except:
                self.logger.info('Created Statistic file')
                print('Not here, need to create it')
                self.data10 = data10
                
    def _calculateVariances(self,data10,):
        cols = ['T [K]','u [m/s]','v [m/s]','w [m/s]']
        from itertools import combinations_with_replacement
        for c1,c2 in combinations_with_replacement(['T [K]','u [m/s]','v [m/s]','w [m/s]'],2):
            name = f'cov{c1.split()[0]}{c2.split()[1]}'
            data10[name] = cross_variance_linear_DataFrame(data10, c1,c2)
        
            
        
        
    
    
    def savedata(self,   typ= 'filtered'):
        

        if 'filtered' in typ and not typ == 'filteredby':
            
            path = self.dpath
            df_qc = self.data.round(2)#.copy()
            
            df_qc.drop([ '10minAVG vertical [m/s]','U_hor','WD'], axis = 1, inplace = True)#U_hor = "U_hor [m/s]",WD = "WD [deg]"
            df_qc.rename(columns = dict(T= 'T [K]',u_met='u [m/s]',v_met='v [m/s]', w_met = 'w [m/s]'), inplace = True)
            
            if typ == 'filtered':
                df_qc = df_qc[['u [m/s]','v [m/s]','w [m/s]','T [K]']]
            
            for hour in tqdm(np.unique(df_qc.index.hour), desc = 'Writing files:'):
                ret  = df_qc[df_qc.index.hour == hour]
                if not ret.index.is_unique:
                    self.logger.error('UniqueError: Index is not unique but doubled! Check')
                if not ret.index.is_monotonic_increasing:
                    self.logger.error('Monotonic Error: Index is not monotonic increasing! Check')
                ## write file...
                if len(ret)>100:
                    retp = ret.copy()#.reindex()
                    retp['Time'] = retp.index.astype(str)
                    cs = ['Time']+[c for c in retp.columns if c !='Time' ]
                    df_pa_table = pa.Table.from_pandas(retp[cs],preserve_index=False)
                    
                    fname = f'{self.model}_{ret.index[0].strftime("%Y%m%d_%H%S")}.dat'
                    self.logger.info( f'Writing file:  {fname}')

                    csv.write_csv(df_pa_table, path + fname)

        elif typ == 'All_Stats':
            
            self.logger.info('Writing All STATS file')
            path = self.basepath
            fname = f'{path}{self.model}_AVG.dat'
            ret = self.data10.round(3)
            if not ret.index.is_unique:
                ret = self._makeIndexUnique(ret)
            
            if not ret.index.is_monotonic:
                ret = self._makeIndexMonotonic(ret)
                
                
            ret.to_csv(fname)
            
        elif typ == 'Daily_Stats':
            self.logger.info('Writing Daily Stats data')
            path = self.basepath+ 'Statistics\\'
            
            for date in tqdm(np.unique(self.data10.index.date)):
                ret = self.data10[self.data10.index.date == date].round(3)
                fname = f'{path}{self.model}_AVG_{ret.index[0].strftime("%Y%m%d")}.dat'
                ret.to_csv(fname)#f'Testfeld_{self.model}_AVGdata_{ret.index[0].strftime("%Y%m%d")}_{ret.index[-1].strftime("%Y%m%d")}.dat')
            
        elif typ == 'filteredby':
            self.logger.info('Writing FilteredBy file')
            path = self.basepath
            ret = self.filteredby
            ret.index.name = 'TIMESTAMP'
            ret.to_csv(f'{path}{self.model}_FilteredBy{self.start.strftime("%Y%m%d")}-{self.end.strftime("%Y%m%d")}.dat')#f'Testfeld_{self.model}_AVGdata_{ret.index[0].strftime("%Y%m%d")}_{ret.index[-1].strftime("%Y%m%d")}.dat')
            
        else:
            raise ValueError(f'Unkknown typ -{typ}-')
        self.logger.info('Finished Writing file')
            
    def plot_things(self, plot_PSD = False):
        savepath = self.figpath
        self.logger.info('Saving figures...')
        try:
            if plot_PSD:
                plot_PSD([self.data[['10minAVG streamwise [m/s]', '10minAVG crosswise [m/s]','10minAVG vertical [m/s]']].dropna()], PSD_kwargs = dict(T_window_sec = 60*2, sf = 20),PSD_type ='fS')
            if savepath is not None:
                print('Saving figures...')
                if plot_PSD:
                    plt.gcf().savefig(savepath + f'PSD_{self.times}.jpg')
                    # plt.gcf().clear()
                    # plt.close(plt.gcf())
                self.qcfig.savefig(savepath + f'QC_{self.times}.jpg')
            # self.qcfig.clear()
            # plt.close(self.qcfig)
                
        except Exception as e:
            self.logger.error('Saving Figures did not work')
            self.logger.error(e)
        
    # def plot_raw_stats(self, col = 'gill_115_u'):
        
    #     plt.figure()
    #     plt.plot(self.data_raw[col], label = 'raw')
    #     plt.plot(self.data[col], label = 'filtered 4$\sigma$')
    #     plt.plot(self.mov_med[col], label = 'Moving Median')
    #     plt.plot(self.data10[col], label = 'averaged')
    #     plt.ylim(-10,25)
    #     plt.legend()

if __name__ == '__main__':
    import time
    s = time.time()
    mast, turbine,LiDAR_mast = read_reference_Data() 
    
    plt.close('all')
    #path = r'C:\Users\meypau\OneDrive - Fraunhofer\Desktop\07_Testfeld\SonicTry\Data\Raw\\'
    #logfile = r'C:\Users\meypau\OneDrive - Fraunhofer\Desktop\01_EMUwind\00_Project Files\Deliverabe C.8 (Testfeld Data)\Scripts\log.log'
    height = '110m'
    height2 = '55m'
    import gc

    plot = True
    basepath = r'Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\post_processed\EMUwindQC\\'
    rawpath = r'Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\\'
    # rawpath = r'D:\Sonic\55m\\'
    # basepath = r'D:\Sonic\\'

    #idx = pd.date_range('2021-11-01', end = '2022-02-01', freq = '12h')
    idx = pd.date_range('2021-03-01', end = '2021-06-01', freq = '12h')
    
    #idx = pd.date_range('2021-11-04', end = '2021-11-06', freq = '12h')
    # idx = pd.date_range('2021-11-01', end = '2021-11-01 12:00', freq = '12h')
    # idx = pd.date_range('2021-11-01', end = '2021-11-20', freq = '12h')
    

    
    
    sd110 = Sonic('Gill110', hor_rotation = -121.31, height = '110m', basepath = basepath)#-121.31-90
    sd55 = Sonic('Gill55', hor_rotation = -121.31, height = '55m', basepath = basepath)
    alle = range(len(idx)-1)
    
    reference = pd.concat([mast[['WS_Mast_SE_115m','WD_Mast_111m']], LiDAR_mast['Z_LiDAR_115m']], axis = 1).dropna()
    reference.columns = ['WS','WD','Z']
    
    for i in alle:
        sd110.logger.info(f'Period {i}/{len(idx)}')
        st,en = idx[i].strftime('%Y_%m_%d_%H%S'),idx[i+1].strftime('%Y_%m_%d_%H%S')
        print('-------------------------------------------')
        print(f'{(i/len(alle)*100):.1f} %')
        print(f'Time since start: {round((time.time()-s)/60,2)}min')
        print(f'Height: {height},  Period {i+1}/{len(alle)} ')
        print(f'Time: {st} -- {en}')

        sd110.read_raw(rawpath, device = 'gill_115', start = st, end = en, overlap = '1h')
        if len(sd110.data_raw)>0:
            sd110.filter_raw_data(plot_reference= True, reference =reference )
            sd110.calc_parameters()
            if plot:
                sd110.plot_things()
            sd110.savedata()      
        
            print(f'Height {height2}')
            sd55.from_raw(sd110, device = 'gill_55')
            sd55.filter_raw_data(plot_reference= True, reference =reference)
            sd55.calc_parameters()
            if plot:
                sd55.plot_things()
            sd55.savedata()
        else:
            print('No data available')
            
        # sd110.data_raw = ''
        # sd110.data = ''
        # sd55.data_raw= ''
        # sd55.data= ''
        # gc.collect()
        # # except Exception as e:
        # #     print(e)
        # #     sd.logger.error('Major Error! Skipped files')
        # #     sd.logger.error(e)
            
    sd110.savedata(typ ='stats')
    sd55.savedata(typ ='stats')
    
            
    
