% function [vector_out vector_V]=gill_fluxes_SRWS(X,Y,Z,T,offset,freq,ND)
% [vector_out third_order]=gill_2r_fluxes(X,Y,Z,T,status,offset)
% Full version had: [vector_out despiking_vector flag_sonic_run std_p third_order]
% 20Hz processing: Applies Double Rotation + calculates all relevant
% turbulent fluxes for GILL sonics.
% INPUT:
% [X,Y,Z]: 3D wind coordinates already filtered by STATUS code
% offset: angle from N to boom positive clock-wise
% freq: aquisition frequency, in order to calculate the spectra
% ND: number of decades to average the spectra
% CAUTION #1: Crosswind temperature corection is inside the sonic firmware,
% not in the code
% CAUTION #2: No flow distortion correction (only for METEK).
% CAUTION #3: No detrending, see gill_2r_fluxes_detrend
% CAUTION #4: No despiking, see gill_2r_fluxes_despike
% OUTPUT:
% vector_out []
% third_order []
% fCo_LMS_log [5xND]
% TODO:
% Implement despiking
% Implement filter for repeated values
    sonics = import_sonics('Z:\Projekte\109797-TestfeldBHV\30_Technical_execution_Confidential\TP3\AP2_Aufbau_Infrastruktur\Infrastruktur_Windmessung\02_Equipment\01_Wartung_Messmast_GE-NET_DWG_20190226\Data\UpgradeData\ASCII\TOA5_gill_115m_55m11604_2022_01_26_1800.dat', 5, 72004);
    X = sonics.u115;
    Y = sonics.v115;
    Z = sonics.w115;
    T = sonics.T115;
    offset=121.31;
    freq = 20;
    ND=10;

     %%%%%%%%%%%% Despike data %%%%%%%%%%%%%%%%%%%%%
    [X ~]=despiking_number(X);
    [Y ~]=despiking_number(Y);
    [Z ~]=despiking_number(Z);
    [T ~]=despiking_number(T);
    
    % Calculate WS and WD with offset
    % UCAR ISFS
    Vaz=offset+90; % for Gill
    %Vaz=offset-90; % for METEK
    cosWD=cos(atan2(-X(1),-Y(1)));
    sinWD=sin(atan2(-X(1),-Y(1)));
    dir=mod((180/pi)*atan2(nanmean(sinWD),nanmean(cosWD))+Vaz,360); 
    Spd=(hypot(X,Y));
    % Get wind speed components
    U=-Spd.*sind(dir);
    V=-Spd.*cosd(dir);        
    W=Z;
    
    WD=mod(nanmean(180/pi * (atan2(-X,-Y))+Vaz),360);
    Uh=nanmean(hypot(X,Y));
    Uvec=nanmean(sqrt(X.^2+Y.^2+Z.^2)); 
    tilt=nanmean(atan2(Z,hypot(X,Y)))*180/pi; % inflow angle!
    % Wind vector STDs
    Uhstd=nanstd(hypot(X,Y));
    Uvecstd=nanstd(sqrt(X.^2+Y.^2+Z.^2)); 
    % Output wind vector characteristics: mean and std
    windvec=[Uh,Uvec,WD,tilt];
    windvecstd=[Uhstd,Uvecstd];
    
    %%%---compare radial speeds
    [V1,V2,V3] = to_rad_speed_SRWS(U,V,W,0,0); % angles of beams inside of to_rad_speed
    
    mU=nanmean(U);
    mV=nanmean(V);
    mW=nanmean(W);
    mT=nanmean(T);

    media = [mU, mV, mW, mT];
    maxia = [max(U), max(V), max(W), max(T)];
    minia = [min(U), min(V), min(W), min(T)];
    
    mV1=nanmean(V1);
    mV2=nanmean(V2);
    mV3=nanmean(V3);
    
    mediaV = [mV1, mV2, mV3];
    maxiaV = [max(V1), max(V2), max(V3)];
    miniaV = [min(V1), min(V2), min(V3)];
    variaV = [cross_variance_linear(V1,V1), cross_variance_linear(V2,V2), cross_variance_linear(V3,V3)];
    
    vector_V=[mediaV,maxiaV,miniaV,variaV];
    
    % Reynolds stress tensor in the geographical coordinate system
    Reynolds_stress=[cross_variance_linear(U,U),cross_variance_linear(U,V),cross_variance_linear(U,W),cross_variance_linear(V,V),cross_variance_linear(V,W),cross_variance_linear(W,W)];
    
    % Heat fluxes in the sonic coordinates (aligned with gravity)
    %heat_fluxes=[cross_variance_linear(X,T),cross_variance_linear(Y,T),cross_variance_linear(Z,T),cross_variance_linear(T,T)];
     
    % Heat fluxes in the 2R coordinate system
    heat_fluxes=[cross_variance_linear(U,T),cross_variance_linear(V,T),cross_variance_linear(W,T),cross_variance_linear(T,T)];
    

    vector_out=[windvec,media,maxia,minia,Reynolds_stress,heat_fluxes,windvecstd];
