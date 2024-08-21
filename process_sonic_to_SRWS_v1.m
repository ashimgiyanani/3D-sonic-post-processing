% Read 1-min SRWS data
clear all;clc;

% Enter the period in UTC
periodini=datenum(2021,3,27,0,0,0); % including this
periodend=datenum(2022,3,27,0,10,0); % not including this

period=datenum(periodini);
freq=20; % 20Hz sampling frequency
% Choose time length
tlength=60; % 1min ensemble means
ND=10;
% offset for each sonic
% for the 110-m
offset110=121.31;

% Projected time series onto SRWS beams
vars={'v1 [m/s]', 'v2 [m/s]', 'v3 [m/s]', ...
    'v1_max [m/s]', 'v2_max [m/s]', 'v3_max [m/s]', ...
    'v1_min [m/s]', 'v2_min [m/s]', 'v3_min [m/s]', ...
    'cov_v1v1 [m2/s2]', 'cov_v2v2 [m2/s2]', 'cov_v3v3 [m2/s2]'
    };

%%
index=1;
%vector_sonics_geo=[]; vector_sonics=[]; spectra_sonics={}; timestamp=[];
timestamp=[];
vector_sonics110=[]; vector_sonics110_LOS=[];

tic

while period<datenum(periodend)

% Catch datevec
    aa=datevec(period);
    
    try
        % Locate the hourly file
        folder = 'Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII';
        %filename=dir(strcat('C:\transfer\*gill*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        %filename=dir(strcat('\\iwes.fraunhofer.de\Data\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\*gill*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        filename=dir(strcat('Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\*gill*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
        if length(filename)==2
            filename=filename(1);
        else
            if length(filename)==0
            vector_sonics55(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
            vector_sonics110(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
            % STAMP AT THE START
            timestamp(index:index+(3600/tlength)-1,1)=[period:(tlength/(24*3600)):period+(3600/(25*3600))];
            period=period+(3600/(24*3600)); % jump to the next hourly file
            index=index+(3600/tlength);
            end
        end
    CATCH
         vector_sonics55(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
         vector_sonics110(index:index+(3600/tlength)-1,:)=ones(3600/tlength,28).*nan; 
         % STAMP AT THE START
         timestamp(index:index+(3600/tlength)-1,1)=[period:(tlength/(24*3600)):period+(3600/(25*3600))];
         period=period+(3600/(24*3600)); % jump to the next hourly file
         index=index+(3600/tlength);
    end
    if datestr(aa)==datestr(datevec(period))
        filename.folder = folder;
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
            % for the 55-m sonic
            X=B(F,12);
            Y=B(F,13);
            Z=B(F,14);
            T=B(F,17)+273.15;
            % Stamped in the START of the period
            timestamp(index,1)=period+(k-1)*(tlength/(24*3600));
            % Stamped in the END of the period
            %timestamp(index,1)=period+(k)*(tlength/(24*3600));
            
            % for the 110-m sonic
            X=B(F,4);
            Y=B(F,5);
            Z=B(F,6);
            T=B(F,9)+273.15;
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
                [MF,MP]=gill_fluxes_SRWS(X,Y,Z,T,offset110,freq,ND); 
                vector_sonics110(index,:)=MF; % 28   
                vector_sonics110_LOS(index,:)=MP; % 12
                index=index+1;
                else
                vector_sonics110(index,:)=ones(1,28).*nan; % 28
                vector_sonics110_LOS(index,:)=ones(1,12).*nan; % 12
                index=index+1;
               end
            
        end
      
    end


TTLOS=timetable(datetime(datevec(timestamp)),vector_sonics110_LOS(:,1),...
    'VariableNames',{'v1 [m/s]'});
for i=2:length(vars)

    TTLOS=addvars(TTLOS,vector_sonics110_LOS(:,i),'NewVariableNames',{vars{i}});
end

% Save the 1-min data in Z
%writetimetable(TT55,strcat('\\172.29.13.76\daten\raw\DB_AD8_METMAST_EXTENSION\DATA\processed_sonics\gill_55m_20Hz\',datestr(time(1),'YYYYmmddhh'),'00_gill_55m_20Hz'));

clear TTLOS


codetime=toc;
disp(['Sonic data processed for',' ',datestr(period,'dd.mm.yyyy HH:MM'),' at ',num2str(codetime),'s'])        
period=datenum(datetime(datevec(period))+hours(1));
clear filename

end
disp(['Fast sonic post-processing: DONE'])

