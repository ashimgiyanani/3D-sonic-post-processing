function [vector_out third_order fCo_LMS_log]=thies_2r_fluxes_tauern(X,Y,Z,T,offset,freq,ND)
% [vector_out third_order]=thies_2r_fluxes(X,Y,Z,T,status,offset)
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
% third_order []
% fCo_LMS_log [5xND]

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Third-order moments
    third_order(1)=mean((U-mU).^3);                  %uuu
    third_order(2)=mean((U-mU).^2.*(V-mV));          %uuv
    third_order(3)=mean((U-mU).^2.*(W-mW));          %uuw
    third_order(4)=mean((U-mU).^2.*(T-mT));          %uut
    
    third_order(5)=mean((U-mU).*(V-mV).^2);          %uvv
    third_order(6)=mean((U-mU).*(V-mV).*(W-mW));     %uvw
    third_order(7)=mean((U-mU).*(V-mV).*(T-mT));     %uvt
    
    third_order(8)=mean((U-mU).*(W-mW).^2);          %uww
    third_order(9)=mean((U-mU).*(W-mW).*(T-mT));     %uwt
    
    third_order(10)=mean((U-mU).*(T-mT).^2);         %utt
    
    third_order(11)=mean((V-mV).^3);                 %vvv
    third_order(12)=mean((V-mV).^2.*(W-mW));         %vvw
    third_order(13)=mean((V-mV).^2.*(T-mT));         %vvt
    third_order(14)=mean((V-mV).*(W-mW).^2);         %vww
    third_order(15)=mean((V-mV).*(W-mW).*(T-mT));    %vwt
    third_order(16)=mean((V-mV).*(T-mT).^2);         %vtt
    
    third_order(17)=mean((W-mW).^3);                 %www
    third_order(18)=mean((W-mW).^2.*(T-mT));         %wwt
    third_order(19)=mean((W-mW).*(T-mT).^2);         %wtt
    
    third_order(20)=mean((T-mT).^3);                 %ttt
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copectra calculation
    %%%--- In general for the FFT ---%%%
    N=length(X);                %%%--- Length of the series
    %N=10800; % 18 HZ always, even if some values are filtered out
    T=N/freq;                 %%%--- Total time series
    dt=T/N;
    time=0:dt:T-dt;
    f=(0:N-1)/(N*dt);
    df=1/T;
    %%%--- for logartithmically averaging the spectra
    fs=f(2:floor(N/2));
    x1=logspace(log10(fs(1)),log10(fs(end)),ND);
    psi_log=zeros(5,ND-1);
    x2=zeros(1,ND-1);
 
 
%%%--- FFT work on sonic data
        uf=fft(U'-mU')/N;
        vf=fft(V'-mV')/N;
        wf=fft(W'-mW')/N;
        Couu_LMS=2*conj(uf).*uf/df;
        Couw_LMS=2*conj(uf).*wf/df;
        Covv_LMS=2*conj(vf).*vf/df;
        Covw_LMS=2*conj(vf).*wf/df;
        Coww_LMS=2*conj(wf).*wf/df;
        %%%--- logarithmic average of the spectra
        for j=1:length(x1)-1
            [~,~,bin]=histcounts(fs,[x1(j) x1(j+1)]);
            rows=find(bin==1);
            psi_log(1,j)=mean(Couu_LMS(1,rows));
            psi_log(2,j)=mean(Covv_LMS(1,rows));
            psi_log(3,j)=mean(Coww_LMS(1,rows));
            psi_log(4,j)=mean(Couw_LMS(1,rows));
            psi_log(5,j)=mean(Covw_LMS(1,rows));
            x2(j)=mean(fs(rows));
        end
        fCo_LMS_log(1,:)=x2.*psi_log(1,:);
        fCo_LMS_log(2,:)=x2.*psi_log(2,:);
        fCo_LMS_log(3,:)=x2.*psi_log(3,:);
        fCo_LMS_log(4,:)=x2.*psi_log(4,:);
        fCo_LMS_log(5,:)=x2.*psi_log(5,:);
        

    vector_out=[windvec,media,maxia,minia,Reynolds_stress,heat_fluxes,windvecstd];
