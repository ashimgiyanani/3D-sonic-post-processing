function [y I]=despiking3_number(x);
% despiking_number3(x): despike single points and linear interpolate voids
% This function is required for despike_number
% INPUT: x, signal to be processed
% OUTPUT: y, despiked signal. I, indexes of spikes
%Source: Ebba >> SEGALINI,KTH MECHANICS, 2011

    y=x;
    x1=x(1:end-1)-x(2:end);
    % if flag==0
        sigma=std(x1);
    % else
    %     sigma=sqrt(median((x1-median(x1)).^2));
    % end
    z=find(abs(x1)>6*sigma);
    I=[];
    if length(find(diff(z)==1))>0
        z1=z(find(diff(z)==1)+1);
        for i=1:length(z1)
            y(z1(i))=(x(z1(i)+1)+x(z1(i)-1))/2;
        end
        I=z1;
    end
