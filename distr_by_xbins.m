function [y_median,y_quantiles,xbins,nvals] = distr_by_xbins(x,y,plot_quantiles,nbins)
% For paired vectors x and y, compute the median and quantiles of the
% y corresponding to the binned ranges of x. The bins of x are determined
% based on the distribution of values in x, so that each bin has an equal
% number of datapoints.
% Nate Zuk (2020)

if nargin < 3 || isempty(plot_quantiles)
    plot_quantiles = [0.1 0.9]; % two values indicating the quantiles of the distributions to plot
end

if nargin < 4 || isempty(nbins)
    nbins = 20;
end

% order values of x, and determine the bin ranges
srtx = sort(x);
srtx(isnan(srtx)) = []; % skip NaN values
bin_idx = round(linspace(1,length(srtx),nbins+1)); % even spacing of indexes
xbins = [srtx(bin_idx(1:end-1)); srtx(bin_idx(end))+min(diff(srtx(bin_idx)))];
    % add extra on the end to include the largest element of x

y_quantiles = NaN(nbins,2);
y_median = NaN(nbins,1);
nvals = NaN(nbins,1);
for n = 1:nbins
    idx = x>=xbins(n) & x<xbins(n+1);
    y_quantiles(n,1) = quantile(y(idx),plot_quantiles(1));
    y_quantiles(n,2) = quantile(y(idx),plot_quantiles(2));
    y_median(n) = median(y(idx));
    nvals(n) = sum(idx);
end