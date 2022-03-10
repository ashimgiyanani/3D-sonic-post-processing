function [vector_out]=thies_2r_fluxes_tauern(X,Y,Z,T,offset,freq,ND)
% [vector_out]=thies_2r_fluxes(X,Y,Z,T,status,offset)
% Full version had: [vector_out despiking_vector flag_sonic_run std_p third_order]
% 20Hz processing: Applies Double Rotation + calculates all relevant
% INPUT:
% [X,Y,Z]: 3D wind coordinates already filtered by STATUS code
% offset: angle from N to boom positive clock-wise
% freq: aquisition frequency, in order to calculate the spectra
% ND: number of decades to average the spectra
% CAUTION #3: No detrending, see thies_2r_fluxes_detrend
% CAUTION #4: No despiking, see thies_2r_fluxes_despike
% OUTPUT:
% vector_out []

    %%%%%%%%%%%% Despike data %%%%%%%%%%%%%%%%%%%%%
    %[X ~]=despiking_number(X);
    %[Y ~]=despiking_number(Y);
    %[Z ~]=despiking_number(Z);
    %[T ~]=despiking_number(T);
    %%%%%%%%%%%%%  data correction  %%%%%%%%%%%%%%%%%%
    % Remove 2D correction from METEK
    %[X,Y,Z]=Remove2DCorrection(X,Y,Z);
    % Apply correction on T
    %T = T+0.00248139*(0.75*(X.^2+ Y.^2)+0.5*Z.^2);
    % Apply 3D correction on METEK
    %[X,Y,Z]=Metek3DCorrection(X,Y,Z);
    

    Xm=mean(X);
    Ym=mean(Y);
    Zm = mean(Z);

    % UCAR ISFS
    %Vaz=offset+90; % for Gill
    %Vaz=offset-90; % for METEK
    Vaz=180-offset; % for Thies
    WD=mod(nanmean(rad2deg(atan2(-Y,X))+Vaz),360);
    % Alternative
    %WD=mod(rad2deg(atan2(-Xm,-Ym))+Vaz,360);
    % Alfredo
    %cosWD=cos(atan2(Y,X));
    %sinWD=sin(atan2(Y,X));        
    %WD=mod((180/pi)*atan2(nanmean(sinWD),nanmean(cosWD))+180+offset,360); %%%--- 20deg line offset!
    % Ebba
    % hor_direction=mod(atan2(Ym,Xm)*180/pi + 180 + 15,360); %%%--- 15deg line offset!
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
    % Triple Rotation (optional, not reccomended):
    %gama=0.5*atan(2*mean((U2(2,:)-mean(U2(2,:))).*(U2(3,:)-mean(U2(3,:))))/(mean((U2(2,:)-mean(U2(2,:))).^2)-mean((U2(3,:)-mean(U2(3,:))).^2)));
    %R23=[1 0 0;0 cos(gama) sin(gama);0 -sin(gama) cos(gama)];
    %U3=R23*U2;
    U3=U2;
    %%--- Rotation by Ebba
    %
    %     vx=Xm/norm([Xm Ym]);
    %     vy=Ym/norm([Xm Ym]);
    %     Uxy=X*vx+Y*vy;
    %     V=X*vy-Y*vx;
    % 
    %     psi = atan2(Zm,mean(Uxy));
    %     U = cos(psi)*Uxy + sin(psi)*Z;
    %     W = -sin(psi)*Uxy + cos(psi)*Z;
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
    
    % Heat fluxes in the sonic coordinates (aligned with gravity)
    %heat_fluxes=[cross_variance(X,T),cross_variance(Y,T),cross_variance(Z,T),cross_variance(T,T)];
    
    % Heat fluxes in the 2R coordinate system
    heat_fluxes=[cross_variance_linear(U,T),cross_variance_linear(V,T),cross_variance_linear(W,T),cross_variance_linear(T,T)];
    
    % Corrected heat flux according to Liu et. al. (2001), in the sonic
    % coordinte system (aligned with gravity)
    %heat_fluxes(3)=heat_fluxes(3)+(2*(mT+273.15)/(331.3+0.606*mT)^2)*(Reynolds_stress(3)*Xm*0.75+Reynolds_stress(5)*Ym*0.75);           
    

    vector_out=[windvec,media,maxia,minia,Reynolds_stress,heat_fluxes,windvecstd];
