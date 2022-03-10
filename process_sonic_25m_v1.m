%% Process sonic data BHV
% Author: Pedro Santos
% contact: pedro.santos@iwes.fraunhofer.de
% ###################################################

clear all;clc;

% Enter the period in UTC
periodini=datenum(2021,12,1,0,0,0); % including this
periodend=datenum(2021,12,2,0,0,0); % not including this
%periodini=datenum(2021,4,1,0,0,0); % including this
%periodend=datenum(2021,4,1,12,0,0); % not including this

period=datenum(periodini);
freq=20; % 20Hz sampling frequency
ND=20; % number of decades for log averaging spectra
% Choose time length
tlength=600; % 10min ensemble means
%tlength=1800; % 30min ensemble means

% offset for each sonic
% for the 25-m
offset25=122.04+90;


index=1;
%vector_sonics_geo=[]; vector_sonics=[]; spectra_sonics={}; timestamp=[];
vector_sonics25=[]; timestamp=[];

%%
tic

while period<datenum(periodend)

% Catch datevec
    aa=datevec(period);
    
    try
        % Locate the hourly file
        filename=dir(strcat('Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\*thies*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        %filename=dir(strcat('Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\*gill*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        %filename=dir(strcat('D:\DATA_IWES\sonic\thies\*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        %filename=dir(strcat('D:\DATA_IWES\sonic\gill\*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        %filename=dir(strcat('C:\Users\sanped\OneDrive - Fraunhofer\Dokumente\0.projects\2.Tauern\1.Campaign\sonic\20Hz_Mast69_SN15670\',num2str(aa(1)),num2str(aa(2),'%02d'),num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        if length(filename)==2
            filename=filename(1);
        else
            if length(filename)==0
            vector_sonics25(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
            % STAMP AT THE START
            timestamp(index:index+(3600/tlength)-1,1)=[period:(tlength/(24*3600)):period+(3600/(25*3600))];
            period=period+(3600/(24*3600)); % jump to the next hourly file
            index=index+(3600/tlength);
            end
        end
    CATCH
         vector_sonics25(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
         % STAMP AT THE START
         timestamp(index:index+(3600/tlength)-1,1)=[period:(tlength/(24*3600)):period+(3600/(25*3600))];
         period=period+(3600/(24*3600)); % jump to the next hourly file
         index=index+(3600/tlength);
    end
    if datetime(aa)==datetime(datevec(period))
        source=strcat(filename.folder,'/',filename.name);
        % read the hourly file 
        A=readmatrix(source,'Delimiter',',','OutputType','char','NumHeaderLines',4);
        time=datenum(table2array(cell2table(A(:,1))));
        B=readmatrix(source,'Delimiter',',','Whitespace','ABCDEF','NumHeaderLines',4);
        %loop over each 10min or 30min period
        for k=1:(3600/tlength)
          % Find 10min or 30min lines of data
            F=find(time>=period+(k-1)*(tlength/(24*3600)) & time<period+k*(tlength/(24*3600)));
            id=B(F,2);
            % Catch X, Y, Z and T components in sonic coordinates
            % for the 25-m sonic
            X=B(F,3);
            Y=B(F,4);
            Z=B(F,5);
            T=B(F,6)+273.15;
            % Stamped in the START of the period
            timestamp(index,1)=period+(k-1)*(tlength/(24*3600));
            % Stamped in the END of the period
            %timestamp(index,1)=period+(k)*(tlength/(24*3600));
            
            if length(X)>freq*tlength*.99 % accept only 99% of samples
                %Despike and fill outliers for sonic data
                X=filloutliers(X,'nearest','movmedian',1200);
                Y=filloutliers(Y,'nearest','movmedian',1200);
                Z=filloutliers(Z,'nearest','movmedian',1200);
                T=filloutliers(T,'nearest','movmedian',1200);
                % Fill remaining NaNs
                X=fillmissing(X,'linear','EndValues','nearest');
                Y=fillmissing(Y,'linear','EndValues','nearest');
                Z=fillmissing(Z,'linear','EndValues','nearest');
                T=fillmissing(T,'linear','EndValues','nearest');
                %process sonics 20Hz
                % Process sonic in 2-rotation coordinates
                [MF]=thies_2r_fluxes_BHV_sample(X,Y,Z,T,offset25,freq,ND); 
                vector_sonics25(index,:)=MF; % 28   
                index=index+1;
            else
                vector_sonics25(index,:)=ones(1,28).*nan; % 28
                index=index+1;
            end
            
            
        end
      
    end
      
codetime=toc;
disp(['Sonic data processed for',' ',datestr(period,'dd.mm.yyyy HH:MM'),' at ',num2str(codetime),'s'])        
period=datenum(datetime(datevec(period))+hours(1));
clear filename

end
disp(['Fast sonic post-processing: DONE'])


%% Create full time series

% means on un-rotated coordinate system

vars={'U_horz [m/s]', 'U_vec [m/s]', 'wind direction [deg]', 'inflow angle [deg]',...
    'u [m/s]', 'v [m/s]', 'w [m/s]', 'T [K]', ...
    'u_max [m/s]', 'v_max [m/s]', 'w_max [m/s]', 'T_max [K]', ...
    'u_min [m/s]', 'v_min [m/s]', 'w_min [m/s]', 'T_min [K]', ...
    'cov_uu [m2/s2]', 'cov_uv [m2/s2]', 'cov_uw [m2/s2]',...
    'cov_vv [m2/s2]', 'cov_vw [m2/s2]', 'cov_ww [m2/s2]',...
    'cov_uT [mK/s]', 'cov_vT [mK/s]', 'cov_wT [mK/s]','cov_TT [K2]',...
    'U_horz std [m/s]', 'U_vec std [m/s]',...
    'u_star []', '1/L [1/m]', 'zL [-]'
    };

% ustar with along and cross-wind components
vector_sonics25(:,29)=(sqrt(vector_sonics25(:,15).^2+vector_sonics25(:,17).^2)).^0.25;

% Inverse of Obukhov length (1/L)
vector_sonics25(:,30)=((-(vector_sonics25(:,29).^3).*(vector_sonics25(:,4)))./(0.4*9.81*vector_sonics25(:,21))).^(-1);


% zL considering the height of the sonic
vector_sonics25(:,31)=25*vector_sonics25(:,30);

TT25=timetable(datetime(datevec(timestamp)),vector_sonics25(:,1),...
    'VariableNames',{'U_horz [m/s]'});


for i=2:length(vars)

    TT25=addvars(TT25,vector_sonics25(:,i),'NewVariableNames',{vars{i}});
    
end

