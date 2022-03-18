function [y I]=despiking_number(x);
% despiking_number(x): despike data and linear interpolate voids
% INPUT: x, signal to be processed
% OUTPUT: y, despiked signal. I, indexes of spikes
%Source: Ebba >> SEGALINI,KTH MECHANICS, 2011

%This part removes all the spikes happening in a single point
    [y I]=despiking3_number(x);

    %This part removes all the spikes happening in two points
    l_even=[1:floor(length(x)/2)]*2; l_odd=l_even-1;
    [y(l_even) I1]=despiking3_number(y(l_even)); I=union(I,l_even(I1));
    [y(l_odd) I1]=despiking3_number(y(l_odd)); I=union(I,l_odd(I1));

    [y I1]=despiking3_number(y);     I=union(I,I1);
    [y I1]=despiking3_number(y);     I=union(I,I1);