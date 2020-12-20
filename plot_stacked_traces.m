function [yt,h] = plot_stacked_traces(xvals,traces,cvals)
% Plot each column of traces stacked along the y axis using the x-values
% specified by xvals. The size and coloring of each trace varies based on
% the values provided in cvals.
% Nate Zuk (2020)

max_clr_rgb = [0 0 1]; % for the color value corresponding to the maximum cval
min_clr_rgb = [0.5 0.5 0.5]; % for the color value corresponding to the minimum cval
max_line_width = 4; % line width for maximum cval
min_line_width = 1; % line width for minimum cval

ntraces = size(traces,2);

% Position the traces so that they are two standard deviations apart, on
% average
ysep = 3*std(reshape(traces,[numel(traces),1]));
yt = ysep*(0:ntraces-1);

% Compute the gradient of colors by linearly interpolating between
% max_clr_rgb and gray ([0.5 0.5 0.5])
clr_vals = NaN(ntraces,3);
line_widths = NaN(ntraces,1); % width is largest for 
for ii = 1:size(traces,2)
    c_prop = (cvals(ii)-min(cvals))/(max(cvals)-min(cvals)); 
        % get the proportion of the current cvalue relative to min and max
        % cvalues
    clr_vals(ii,:) = (max_clr_rgb-min_clr_rgb)*c_prop+min_clr_rgb;
    line_widths(ii) = (max_line_width-min_line_width)*c_prop+min_line_width;
end

% Create the plot
h = figure;
hold on
for ii = 1:size(traces,2)
    plot(xvals,traces(:,ii)+yt(ii),'Color',clr_vals(ii,:),'LineWidth',line_widths(ii));
end
set(gca,'FontSize',14,'YTick',yt,'YTickLabel',1:ntraces,'YLim',[yt(1)-2*ysep yt(end)+2*ysep]);