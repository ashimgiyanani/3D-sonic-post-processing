%% Process sonic data BHV
% TO DO:
% Add spectra to the processing

clear all;clc;

% Enter the period in UTC
periodini=datenum(2022,2,1,0,0,0); % including this
periodend=datenum(2022,3,1,0,0,0); % not including this
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

%%
tic

while period<datenum(periodend)

% Restart the output array every hour
vector_sonics25=[]; timestamp=[];
% Restart the index every hour
index=1;
    
% Catch datevec
    aa=datevec(period);
    
    try
        % Locate the hourly file
        filename=dir(strcat('\\iwes.fraunhofer.de\Data\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\*thies*',num2str(aa(1)),'_',num2str(aa(2),'%02d'),'_',num2str(aa(3),'%02d'),'_',num2str(aa(4),'%02d'),'*.dat'));
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
                [MF,~,~]=thies_2r_fluxes_BHV(X,Y,Z,T,offset25,freq,ND); 
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



% ustar with along and cross-wind components
vector_sonics25(:,29)=(sqrt(vector_sonics25(:,15).^2+vector_sonics25(:,17).^2)).^0.25;
% Inverse of Obukhov length (1/L)
vector_sonics25(:,30)=((-(vector_sonics25(:,29).^3).*(vector_sonics25(:,4)))./(0.4*9.81*vector_sonics25(:,21))).^(-1);
% zL considering the height of the sonic
vector_sonics25(:,31)=25*vector_sonics25(:,30);
% Make a timetable
TT25=timetable(datetime(datevec(timestamp)),vector_sonics25(:,1),...
    'VariableNames',{'U_horz [m/s]'});
for i=2:length(vars)
    TT25=addvars(TT25,vector_sonics25(:,i),'NewVariableNames',{vars{i}});
end
%Save the hourly data
%writetimetable(TT25,strcat('.\data\thies_25m\',datestr(timestamp(1),'YYYYmmddhh'),'00_thies_25m'));
% Save the hourly data into OneDAS
writetimetable(TT25,strcat('\\172.29.13.76\daten\raw\DB_AD8_METMAST_EXTENSION\DATA\processed_sonics\thies_25m\',datestr(timestamp(1),'YYYYmmddhh'),'00_thies_25m'));

clear TT25

end
disp(['Fast sonic post-processing: DONE'])


%% Create full time series

%% Load full time series (already gap-filled)



%% compare with met mast

scheme          = 'https';
host            = 'onedas.iwes.fraunhofer.de';
port         	= 443;
% Edit your username here
username        = 'pedro.santos@iwes.fraunhofer.de';
% Enter your password here
password        = 'Mmlp1001p!p!p';

% Choose the period of interest
dateTimeBegin 	= datetime(2021, 12, 1, 0, 0, 0, 'TimeZone', 'UTC');
dateTimeEnd 	= datetime(2021, 12, 2, 0, 0, 0, 'TimeZone', 'UTC');

% must all be of the same sample rate
channelPaths = { ...
    %'/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0030_V4/1 s_mean'
    %'/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0100_D4/1 s_mean_polar'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0030_V4/600 s_mean'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0100_D4/600 s_mean_polar'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0000_V1/600 s_mean'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0010_V2/600 s_mean'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0070_D1/600 s_mean_polar'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0020_V3/600 s_mean'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0040_V5/600 s_mean'
    '/AIRPORT/AD8_PROTOTYPE/GENERAL_DAQ/M0050_V6/600 s_mean'
    };

% load data
connector = OneDasConnector(scheme, host, port, username, password);
% without authentication: connector = OneDasConnector(scheme, host, port)
params.ChannelPaths = channelPaths;
data                = connector.Load(dateTimeBegin, dateTimeEnd, params);

% Save data in a structure
for k=1:length(channelPaths)
    AA{k,1}=data(channelPaths{k});
    vars{k,1}=AA{k,1}.Description;
end

% Create a time table
% Make a timestamp array
%sampleRate  = 1; % 1/600 Hz or 10-min (adapt to your needs)
sampleRate  =  1/600; %Hz or 10-min (adapt to your needs)
dt          = 1 / sampleRate / 86400;
time        = (datenum(dateTimeBegin) : dt : datenum(dateTimeEnd) - dt).';

% Compile the time table
TT=timetable(datetime(datevec(time)),AA{1}.Values,...
    'VariableNames',{AA{1}.Description});

for i=2:length(channelPaths)
    TT=addvars(TT,AA{i}.Values,'NewVariableNames',{AA{i}.Description});
end

%% Combine time series from met mast and sonics

TTfinal=synchronize(TT,TT25);

figure
subplot(2,1,1)
plot(TTfinal.Time,TTfinal.("U_horz [m/s]"))
hold on
plot(TTfinal.Time,TTfinal.("Wind speed 5 / 25 m"))

subplot(2,1,2)
plot(TTfinal.Time,TTfinal.("wind direction [deg]"))
hold on
plot(TTfinal.Time,TTfinal.("Wind direction 4 / 23.5"))


%% Get a stability code for the lowest level (for now is at 55m)
%Classify stability

ttemp=datevec(TTfinal.Time);
FF=find(TTfinal.("U_horz [m/s]")<3);
% LMO in terms of ZL
LMO=TTfinal.("zL [-]");
LMO(FF,:)=nan;

sonicno=1;
neutralth=0.1;
nearth=0.5;
%veryth=0.5;
% Unstable = 1
stabilityLMO{1,1}=find(LMO(:,sonicno)<=-nearth);
% Near unstable =2
stabilityLMO{2,1}=find(LMO(:,sonicno)<=-neutralth & LMO(:,sonicno)>-nearth);
% Neutral = 3
stabilityLMO{3,1}=find(abs(LMO(:,sonicno))<neutralth);
% Near stable = 4
stabilityLMO{4,1}=find(LMO(:,sonicno)>=neutralth & LMO(:,sonicno)<nearth);
% Stable = 5
stabilityLMO{5,1}=find(LMO(:,sonicno)>=nearth);

%% Add stability code on TT
scode=ones(length(TTfinal.Time),1)*nan;
for i=1:5
    scode(stabilityLMO{i,1})=i;
end


%%
line=0:0.1:25;
Y=TTfinal.("U_horz [m/s]");
X=TTfinal.("Wind speed 5 / 25 m");
Z=TTfinal.("wind direction [deg]");
W=TTfinal.("Wind direction 4 / 23.5");

%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270) | W>330 | W<90));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>0 & W<60)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>100 & W<120)));
FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>0 & W<360)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>0 & W<60) | (W>100 & W<120) | (W>210 & W<220)));



figure;
plot(X(FF,1),Y(FF,1),'.','color',[0.5 0.5 0.5]);
set(get(gca,'XLabel'),'Fontsize',14,'Interpreter','latex','String','cup [m s$$^{-1}$$]')
set(get(gca,'YLabel'),'Fontsize',14,'Interpreter','latex','String','sonic [m s$$^{-1}$$]')
title('Cup 55 m vs Sonic 55 m ($190^\circ$ to $220^\circ$)','Fontsize',14,'Interpreter','latex');
%set(gca,'Fontsize',12)
p=polyfit(X(FF,1),Y(FF,1),1);
rho=corr(X(FF,1),Y(FF,1));


hold on
plot(line,line,'k:','LineWidth',1.5);
grid on
slope=num2str(p(1,1));
inter=num2str(p(1,2));
rho2=num2str(round(rho^2,3));

text(15,8,['y=',slope(1,1:4),'x + ',inter(1,1:4)],'Fontsize',14,'Interpreter','latex');
text(15,6,['R$$^2$$=',rho2(1,1:6)],'Fontsize',14,'Interpreter','latex');
text(15,3,['N=',num2str(length(X(FF)))],'Fontsize',14,'Interpreter','latex');
xlim([0 max(line)])
ylim([0 max(line)])

%export_fig cup_sonic_55m.pdf -painters

%%
line=0:0.1:25;
Y=TTfinal.("U_horz [m/s]_TT110");
X=TTfinal.("Wind speed 1 / 115.9 m");
Z=TTfinal.("wind direction [deg]_TT110");
W=TTfinal.("Wind direction  1 / 110.8 m");

%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270) | W>330 | W<90));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>0 & W<60)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>100 & W<120)));
%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>190 & W<220)));
FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>0 & W<60) | (W>100 & W<120) | (W>210 & W<220)));



figure;
plot(X(FF,1),Y(FF,1),'.','color',[0.5 0.5 0.5]);
set(get(gca,'XLabel'),'Fontsize',14,'Interpreter','latex','String','cup [m s$$^{-1}$$]')
set(get(gca,'YLabel'),'Fontsize',14,'Interpreter','latex','String','sonic [m s$$^{-1}$$]')
title('Cup 115.9 m vs Sonic 110 m ($150^\circ$ to $270^\circ$)','Fontsize',14,'Interpreter','latex');
%set(gca,'Fontsize',12)
p=polyfit(X(FF,1),Y(FF,1),1);
rho=corr(X(FF,1),Y(FF,1));


hold on
plot(line,line,'k:','LineWidth',1.5);
grid on
slope=num2str(p(1,1));
inter=num2str(p(1,2));
rho2=num2str(rho^2);

text(15,8,['y=',slope(1,1:4),'x + ',inter(1,1:4)],'Fontsize',14,'Interpreter','latex');
text(15,6,['R$$^2$$=',rho2(1,1:4)],'Fontsize',14,'Interpreter','latex');
text(15,3,['N=',num2str(length(X(FF)))],'Fontsize',14,'Interpreter','latex');
xlim([0 max(line)])
ylim([0 max(line)])

%export_fig cup_sonic_110m.pdf -painters


%%

line=0:0.1:360;
Y=TTfinal.("U_horz [m/s]_TT55");
X=TTfinal.("Wind speed 4 / 55 m");
Z=TTfinal.("wind direction [deg]_TT55");
W=TTfinal.("Wind direction 4 / 23.5");

%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270) | W>330 | W<90));
FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) );%& ((W>150 & W<270)));


figure;
plot(W(FF,1),Z(FF,1),'.','color',[0.5 0.5 0.5]);
set(get(gca,'XLabel'),'Fontsize',14,'Interpreter','latex','String','vane [deg]')
set(get(gca,'YLabel'),'Fontsize',14,'Interpreter','latex','String','sonic [deg]')
title('Vane 23.5 m vs Sonic 55 m','Fontsize',14,'Interpreter','latex');
%set(gca,'Fontsize',12)
p=polyfit(X(FF,1),Y(FF,1),1);
rho=corr(X(FF,1),Y(FF,1));


hold on
plot(line,line,'k:','LineWidth',1.5);
grid on
slope=num2str(p(1,1));
inter=num2str(p(1,2));
rho2=num2str(rho^2);
%slope2=num2str(linearfit);
%text(15,8,['y=',slope(1,1:4),'x + ',inter(1,1:4)],'Fontsize',14,'Interpreter','latex');
%text(15,6,['R$$^2$$=',rho2(1,1:4)],'Fontsize',14,'Interpreter','latex');
%text(15,3,['N=',num2str(length(X(FF)))],'Fontsize',14,'Interpreter','latex');
xlim([0 max(line)])
ylim([0 max(line)])

%export_fig cup_sonic_55m.pdf -painters

%%

line=0:0.1:360;
Y=TTfinal.("U_horz [m/s]_TT110");
X=TTfinal.("Wind speed 1 / 115.9 m");
Z=TTfinal.("wind direction [deg]_TT110");
W=TTfinal.("Wind direction  1 / 110.8 m");

%FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) & ((W>150 & W<270) | W>330 | W<90));
FF=find(~isnan(X) & ~isnan(Y) & ~isnan(Z) & ~isnan(W) );%& ((W>150 & W<270)));


figure;
plot(W(FF,1),Z(FF,1),'.','color',[0.5 0.5 0.5]);
set(get(gca,'XLabel'),'Fontsize',14,'Interpreter','latex','String','vane [deg]')
set(get(gca,'YLabel'),'Fontsize',14,'Interpreter','latex','String','sonic [deg]')
title('Vane 110 m vs Sonic 110 m','Fontsize',14,'Interpreter','latex');
%set(gca,'Fontsize',12)
p=polyfit(X(FF,1),Y(FF,1),1);
rho=corr(X(FF,1),Y(FF,1));


hold on
plot(line,line,'k:','LineWidth',1.5);
grid on
slope=num2str(p(1,1));
inter=num2str(p(1,2));
rho2=num2str(rho^2);
%slope2=num2str(linearfit);
%text(15,8,['y=',slope(1,1:4),'x + ',inter(1,1:4)],'Fontsize',14,'Interpreter','latex');
%text(15,6,['R$$^2$$=',rho2(1,1:4)],'Fontsize',14,'Interpreter','latex');
%text(15,3,['N=',num2str(length(X(FF)))],'Fontsize',14,'Interpreter','latex');
xlim([0 max(line)])
ylim([0 max(line)])

%export_fig cup_sonic_55m.pdf -painters

%% Catch 

V=horzcat(TTfinal.("Wind speed 6 / 10 m"),TTfinal.("Wind speed 5 / 25 m"),TTfinal.("Wind speed 4 / 55 m"),TTfinal.("Wind speed 3 / 85 m"),TTfinal.("Wind speed 1 / 115.9 m"));
height=[10 25 55 85 116];

%% Dimensionless shear with Universal Function 
z1=55;z2=110;
zmed=(z2-z1)/log(z2/z1); % reference height from Pena (2008) BLM paper
ZL=zmed.*TTfinal.x1_L_1_m__TT55(FF);
K=0.4;
fi=(K./TTfinal.u_star___TT55(FF)).*(px{1,1}(:,2)+2*px{1,1}(:,1)*log(zmed)); % Calculate phi with the fitting
%Mean Flux-Profile
%
for i=1:70   % bins for average analysis
clear F; F=find(ZL>(-5+(i-1)*0.1)-0.05 & ZL<(-5+(i-1)*0.1+0.05));
ZLmed(i,:)=nanmean(ZL(F,:));

fimed(i,:)=nanmean(fi(F,:));
fistd(i,:)=nanstd(fi(F,:))./length(fi(F));
end
%}
% Theorerical Profile from Hogstrom (1988)
phit=[];
ZLt=[-2:0.1:2];
a=12;b=4.7;p=-1/3;
% Unstable
phit(find(ZLt<0))=(1-a.*ZLt(ZLt<0)).^(p);
%Neutral
phit(find(ZLt==0))=1;
%Stable
phit(find(ZLt>0))=(1+b.*ZLt(ZLt>0));


%% Dimensionless shear based on 2nd order poly on ln(z) as in Hogstrom (1988)



% ALL
for j=1:4 % for each stability
for k=1:size(FAA{j,1},1); % for each profile
px{j,1}(k,:)=polyfit(log(height),VAA{j,1}(k,:),2); % fitting Phi with a polynomial as in Hogstrom (1988)
end;
for i=1:5 %for each height
fifit{j,1}(:,i)=(K./ustar(FAA{j,1})).*(px{j,1}(:,2)+2*px{j,1}(:,1).*log(height(i))); % Calculate phi with the fitting
end
end

