function c=cross_variance_linear(a,b)
% Cross_variance_linear (a,b)
% Unbiased estimator of the cross-variance between A and B
% Computing the fluctuations from a linear trend +
% It accounts for the sample size on the calculation
    c=nanmean((detrend(a)).*(detrend(b)))*1/(1-1/length(a(~isnan(a))));
