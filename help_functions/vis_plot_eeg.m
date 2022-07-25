function h = vis_plot_eeg(t,data,opts)
% h = plot_eeg(t,data,opts)
% Plots signal of each channel in data. 
%
% INPUTS:
%       t: time vector (length = # of time points)
%       data: matrix of ieeg signals (# channels x # time points)
%       opts (optional): structure containing options (all optional)
%           .offset: offset between plotted signals (default: 200)
%           .plot_labels: true if want to plot electrode names on y-axis
%               (default: false)
%           .plot_numbers: true if want to plot electrode numbers on y-axis
%               (default: false)
%           .labels: cell array of label names; need if plotting labels on
%               y-axis or if want to color lines by electrode groups
%           .clrs: # channels x 3 matrix of colors of lines
%           .clr_groups: false to make each electrode signal a different
%               color (default); true to make each group of electrodes
%               a different color
%           .group_clrs: 3 x # electrode groups specifying color for each
%               group of electrodes (i.e., different color for each
%               grid/strip/shank).
%           .xlab (optional): label for x axis
%           .linewidth: linewidth of time series (default: 1)
%
% OUTPUTS:
%       h: handle to plot
%       
%
% Modified from plotEEG function by Yujiang Wang
%
% Gabrielle M. Schroeder
% Newcastle University School of Computing Science
% 30 June 2017
%
% To do:
% add option for specifying yaxis label font size
% add option for only plotting some labels
% add option for changing beta

% number of channels
n_chan=size(data,1);

% set default options
if nargin<3
    opts.offset=200;
elseif ~isfield(opts,'offset')
    opts.offset=200;
end

% y-axis labels
if ~isfield(opts,'plot_numbers'); opts.plot_numbers=false; end
if ~isfield(opts,'plot_labels'); opts.plot_labels=false; end
if ~isfield(opts,'linewidth'); opts.linewidth=1; end

if opts.plot_numbers && opts.plot_labels
    error('Can label channels with either names or numbers, not both (only one of opts.plot_numbers and opts.plot_labels can be true)')
end

% check that length of labels is correct
if isfield(opts,'labels')
    if length(opts.labels)~=size(data,1)
        error('The number of channel labels does not match the number of channels in data')
    end
end

% colors

if ~isfield(opts,'clr_groups'); opts.clr_groups=false; end

if isfield(opts,'clrs') && (isfield(opts,'group_clrs') || opts.clr_groups==1)    
    error('Ask to plot either group colors or individual signal colors, not both.')

% use specified group colors to determine color of each line
elseif isfield(opts,'group_clrs')    
    [groups,idx]=get_electrode_groups(opts.labels);    
    opts.clrs=zeros(n_chan,3);    
    for i=1:length(groups)
        opts.clrs(idx{i},:)=repmat(opts.group_clrs(i,:),[length(idx{i}) 1]);
    end
    
    % make every other color slightly lighter or darker
    beta=-0.4;
    if beta>0
        opts.clrs(2:2:n_chan)=opts.clrs(2:2:n_chan).^(1-beta);
    else
        opts.clrs(2:2:n_chan)=opts.clrs(2:2:n_chan).^(1/(1+beta)); 
    end
    
% create group colors and assign appropriate one to each line
elseif opts.clr_groups
    [groups,idx]=get_electrode_groups(opts.labels);
    opts.group_clrs=lines(length(groups));
    opts.clrs=zeros(n_chan,3);
    for i=1:length(groups)
        opts.clrs(idx{i},:)=repmat(opts.group_clrs(i,:),[length(idx{i}) 1]);
    end
    
    % make every other color slightly lighter or darker
    beta=-0.4;
    if beta>0
        opts.clrs(2:2:n_chan)=opts.clrs(2:2:n_chan).^(1-beta);
    else
        opts.clrs(2:2:n_chan)=opts.clrs(2:2:n_chan).^(1/(1+beta)); 
    end
    
end

% change offset of channels
offset_data=data;
for i=1:n_chan
    offset_data(i,:)=data(i,:)-nanmean(data(i,:));      % first de-mean so centered at zero
    offset_data(i,:)=offset_data(i,:)-opts.offset*(i-1);    % offset
end

% plot
h=plot(t,offset_data','LineWidth',opts.linewidth);

% change colors
if isfield(opts,'clrs')
    for i=1:n_chan; set(h(i),'color',opts.clrs(i,:)); end
end
axis tight

% plot labels
if opts.plot_labels
%     set(gca,'ytick',(opts.offset*(n_chan-1)*-1):opts.offset:0,'yticklabel',flipud(opts.labels));
    ax=ancestor(h(1),'axes');
    ax.YTick=(opts.offset*(n_chan-1)*-1):opts.offset:0;
    ax.YTickLabel=flipud(opts.labels);
    ax.YAxis.FontSize=5;
% otherwise, plot numbers
elseif opts.plot_numbers
    n_ticks=round(n_chan/5,0)-1;
    ax=ancestor(h(1),'axes');
    ax.YTick=(opts.offset*((n_ticks*5)-1)*-1):(opts.offset*5):(opts.offset*-4);
    ax.YTickLabel=(n_ticks*5):-5:5;
    ax.YAxis.FontSize=8;
else
    ax=ancestor(h(1),'axes');
    ax.YTick=[];
end

% label x-axis
if isfield(opts,'xlab')
    xlabel(opts.xlab)
end
