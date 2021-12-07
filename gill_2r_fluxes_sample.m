function [vector_out]=gill_2r_fluxes_sample(X,Y,Z,T,offset,freq)
% Author: Pedro Santos
% contat: pedro.santos@iwes.fraunhofer.de
% [vector_out]=gill_2r_fluxes_sample(X,Y,Z,T,status,offset,ND)
% 20Hz processing: Applies Double Rotation + calculates all relevant
% INPUT:
% [X,Y,Z]: 3D wind coordinates already filtered by STATUS code
% offset: angle from N to boom positive clock-wise
% ND: number of decades to average the spectra
% OUTPUT:
% vector_out []

    Xm=mean(X);
    Ym=mean(Y);
    Zm = mean(Z);

    % UCAR ISFS
    Vaz=offset+90; % for Gill
    WD=mod(nanmean(rad2deg(atan2(-Y,X))+Vaz),360);
    % Wind vector Means 
    Uh=nanmean(hypot(X,Y));
    Uvec=nanmean(sqrt(X.^2+Y.^2+Z.^2)); 
    tilt=nanmean(atan2(Z,hypot(X,Y)))*180/pi; % inflow angle!
    % Wind vector STDs
    Uhstd=nanstd(hypot(X,Y));
    Uvecstd=nanstd(sqrt(X.^2+Y.^2+Z.^2)); 
    tiltstd=nanstd(atan2(Z,hypot(X,Y)))*180/pi; % inflow angle!
    % Output wind vector characteristics:
    windvec=[Uh,Uvec,WD,tilt];
    windvecstd=[Uhstd,Uvecstd];
    %%--- coordinate rotation:
    alpha=atan2(nanmean(Y),nanmean(X));
    R01=[cos(alpha) sin(alpha) 0;-sin(alpha) cos(alpha) 0;0 0 1];
    U1=R01*[X Y Z]'; % 1 rotation
    beta=atan2(nanmean(U1(3,:)),sqrt(nanmean(U1(1,:))^2+nanmean(U1(2,:))^2));
    R12=[cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)];
    U2=R12*U1; % Double Rotation
    U3=U2;
    %%-- Rotated wind components
    U=U3(1,:)';mU=mean(U3(1,:)');
    V=U3(2,:)';mV=mean(U3(2,:)');
    W=U3(3,:)';mW=mean(U3(3,:)');
    mT=mean(T);
    
    media=  [mU, mV, mW, mT];
    maxia = [max(U), max(V), max(W), max(T)];
    minia = [min(U), min(V), min(W), min(T)];

    % Reynolds stress tensor in the 2R coordinate system
    Reynolds_stress=[cross_variance_linear(U,U),cross_variance_linear(U,V),cross_variance_linear(U,W),cross_variance_linear(V,V),cross_variance_linear(V,W),cross_variance_linear(W,W)];
    
    % Heat fluxes in the 2R coordinate system
    heat_fluxes=[cross_variance_linear(U,T),cross_variance_linear(V,T),cross_variance_linear(W,T),cross_variance_linear(T,T)];        

    vector_out=[windvec,media,maxia,minia,Reynolds_stress,heat_fluxes,windvecstd];
