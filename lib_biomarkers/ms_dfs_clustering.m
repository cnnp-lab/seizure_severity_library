function [tbl_out,parset] = ms_dfs_clustering(meta_tbl,df_tbl,opts)
% takes arbitryary amount of segments and clusters the channels
% of each segment based on a distance measure based on their
% dominant frequencies into clusters. Outputs a variety of measures. Within
% the function, a lot of time series of measures are computed, but these
% are then boiled down to a single value for the table. In general, dfs is
% short for dominant frequency similarity.
% inputs:
%   - meta_tbl: contains meta data about the patients
%   - df_tbl: contains the dominant frequencies and the time points of
%       their measurement 
%   - opts: for more info look into arguments below
% outputs:
%   - tbl_out: table with following measures in the rows: 
%       'segment_id','dfs_allchns','ncl_dfs','pchns_largcl','dfs_largcl',
%       'dfs_dur','df_cl_med','dfs_cl_std','nchns_largcl','dfs_cluster'

arguments
    meta_tbl
    df_tbl
    
    % distance measure
    opts.a = 1
    
    % clustering
    opts.perc_sim = 0.1;
    opts.min_s_sim = 5;
    opts.max_s_sim = 5;
    opts.cluster_dist = 'average'
    
    % plotting
    opts.fig_ind = 10
    opts.plot_rel = false
end
distfct = @(x,y) distfct_a(x,y,opts.a);
    
    
%% main 
% align the inputs by the segment ids
meta_tbl = [array2table((1:size(meta_tbl,1))'),meta_tbl];
tbl = tbl_join(meta_tbl,df_tbl);
tbl = sortrows(tbl,1);
tbl = tbl(:,2:end);

% create variables to efficiently generate the output table
nsegs = size(tbl,1);
dfs_allchns_cell = zeros(nsegs,1);
dfs_ncl_cell = zeros(nsegs,1);
dfs_largcl_cell = zeros(nsegs,1);
nchns_largcl_cell = zeros(nsegs,1);
pchns_largcl_cell = zeros(nsegs,1);
dfs_dur_cell = zeros(nsegs,1);
df_lcl_med_cell = zeros(nsegs,1);
df_lcl_std_cell = zeros(nsegs,1);
dfs_cluster_cell = cell(nsegs,1);

% main loop
for s = 1:nsegs
    
    % preprocess data
    df_data = tbl.df{s};
    t = tbl.t_wndw{s};
    ind = 0 <= t & t < tbl.duration(s);
        % data points that are included into the distance measure, i.e. as
        % of now only the ictal period
    df_data = df_data(:,ind);
    t = t(ind);  
    no_sec_insync = min(opts.max_s_sim,max(opts.min_s_sim,tbl.duration(s)*opts.perc_sim));

    % calc initial clustering 
    nchns = size(df_data,1);
    Z = linkage(df_data,opts.cluster_dist,distfct);
    cutoff_low = 0.00001;
    fs_freq = numel(t) / (t(end)-t(1));
    cutoff = max(cutoff_low,sum(ind) - no_sec_insync*fs_freq);
    if cutoff == cutoff_low
        warning('cutoff is basically 0')
    end
    T = cluster(Z,'cutoff',cutoff,'criterion','distance');
    no_cluster = length(unique(T));
    T_final = T;


    % put all single clusters into one cluster
    no_mult_cl = 0;
    clusters = cell(no_cluster+1,1);
    single_ind = zeros(nchns,1);
    for c = 1:no_cluster
        ind = find(T == c);
        if length(ind) <= 1
            T_final(ind) = -1;
            single_ind(ind) = 1;
        else
            no_mult_cl = no_mult_cl + 1;
            T_final(ind) = no_mult_cl;
            clusters{no_mult_cl} = ind';
        end
    end
    T_final(T_final == -1) = no_mult_cl + 1;  %evtl without max(2,
    if no_cluster > no_mult_cl
        clusters{no_mult_cl+1} = find(single_ind);
    end
    clusters(no_mult_cl+2:end) = [];


    % specify return values
    ndata = length(t);
    cluster_rel_sync_ts = zeros(no_mult_cl,ndata); 
        % normed by nchns in cluster
    cluster_abs_sync_ts = zeros(no_mult_cl,ndata);
        % normed by nchns in segment
    for c = 1:no_mult_cl

        nchns_cl = length(clusters{c});
        res = ms_simfreq_noaverg(df_data(clusters{c},:));
        cluster_rel_sync_ts(c,:) = res / (nchns_cl*(nchns_cl-1) /2);
        cluster_abs_sync_ts(c,:) = res / (nchns*(nchns-1) /2);
    end
    
    % calc dfs all chns
    dfs_allchns_ts = ms_simfreq_noaverg(df_data) / (nchns*(nchns-1) /2);


    % calc med_freq and std_freq in each cluster
    cluster_med_freq_cl = nan(no_mult_cl,1);        % TODO: returns 0 if no cluster?
    cluster_std_freq_cl = nan(no_mult_cl,1);
    for cl = 1:no_mult_cl
        data_cl = df_data(clusters{cl},:);
        cluster_med_freq_cl(cl) = median(median(data_cl,1));
        cluster_std_freq_cl(cl) = std(median(data_cl,1));
    %     cluster_std_freq(cl) = mean(std(data_cl,[],1));
    end


    % reorder clusters and output vals by their relative sync val
    if no_mult_cl > 0

        % determine order
    %     order_val = mean(cluster_rel_sync,2);
        order_val = zeros(no_mult_cl,1);
        for c = 1:no_mult_cl
            order_val(c) = numel(clusters{c});
        end
        [~,I] = sort(order_val,'descend');

        % reorder existing outpus
        cluster_rel_sync_ts = cluster_rel_sync_ts(I,:);
        cluster_abs_sync_ts = cluster_abs_sync_ts(I,:);
        cluster_med_freq_cl = cluster_med_freq_cl(I,:);
        cluster_std_freq_cl = cluster_std_freq_cl(I,:);
        I = [I;no_mult_cl+1];
        clusters = clusters(I);
        for c=1:no_mult_cl
            T_final(clusters{c}) = c;
        end
    end
    
    % calc measures from timeseries above
    no_cl = size(clusters,1)-1;
    dfs_allchns = mean(dfs_allchns_ts);
    if size(cluster_rel_sync_ts,1) == 0
        cluster_med_freq = nan;
        cluster_std_freq = nan;
        dfs_lcl = nan;
        perc_chn = 0;
    else
        cluster_med_freq = cluster_med_freq_cl(1);
        cluster_std_freq = cluster_std_freq_cl(1);
        perc_chn = numel(clusters{1}) / numel(tbl.segment_channel_labels{1});
        dfs_lcl = mean(cluster_rel_sync_ts(1,:));
    end
    
    % calc SWD duration
    if no_cl > 0 && tbl.duration(s) > 0
        ind_ictal = 0 < t & t <= tbl.duration(s);
        t_ictal = t(ind_ictal);
        [val,ind_max] = max(cluster_rel_sync_ts(1,:));
        ind_lwr = find(cluster_rel_sync_ts(1,1:ind_max) < val/2,1,'last');
        ind_hgr = ind_max + find(cluster_rel_sync_ts(1,ind_max+1:end) < val/2,1,'first');
        if isempty(ind_lwr)
            ind_lwr = 1;
        end
        if isempty(ind_hgr)
            ind_hgr = numel(t_ictal);
        end

        SWD_dur = t_ictal(ind_hgr) - t_ictal(ind_lwr);
    else
        SWD_dur = 0;
    end
    
    
    % save calculations
    dfs_allchns_cell(s) = dfs_allchns;
    dfs_ncl_cell(s) = no_cl;
    dfs_largcl_cell(s) = dfs_lcl;
    dfs_dur_cell(s) = SWD_dur;
    df_lcl_med_cell(s) = cluster_med_freq;
    df_lcl_std_cell(s) = cluster_std_freq;
    dfs_cluster_cell{s} = clusters;
    nchns_largcl_cell(s) = numel(clusters{1});
    pchns_largcl_cell(s) = perc_chn;
end
    
% create outputs
% tbl_out = table(tbl.segment_id,dfs_allchns_cell{:},dfs_ncl_cell{:},pchns_largcl_cell{:},...
%     dfs_largcl_cell{:},dfs_dur_cell{:},df_cl_med_cell{:},dfs_cl_std_cell{:},nchns_largcl_cell{:},dfs_cluster_cell,...
%     'VariableNames',{'segment_id','dfs_allchns','ncl_dfs','pchns_largcl','dfs_largcl','dfs_dur','df_cl_med',...
%     'dfs_cl_std','nchns_largcl','dfs_cluster'});
tbl_out = table(tbl.segment_id,dfs_allchns_cell,dfs_ncl_cell,pchns_largcl_cell,...
    dfs_largcl_cell,dfs_dur_cell,df_lcl_med_cell,df_lcl_std_cell,nchns_largcl_cell,dfs_cluster_cell,...
    'VariableNames',{'segment_id','dfs_allchns','ncl_dfs','pchns_largcl','dfs_largcl','dfs_dur','df_lcl_med',...
    'df_lcl_std','nchns_largcl','dfs_cluster'});
parset = opts;


%% plotting

if opts.fig_ind > -1 && nsegs == 1
    if ~strcmp('segment_data',tbl.Properties.VariableNames)
        warning('Need segment_data for plotting results.')
        return
    end
    
    if opts.fig_ind == 0
        figure()
    elseif opts.fig_ind >0.5
        figure(opts.fig_ind)
    end


    % generate colors by clusters
    ncl = size(clusters,1);
    colors = [distinguishable_colors(ncl-1); [.85 .85 .85]];
    nchns = size(tbl.segment_data{1}, 1);
    T = cell2array_cl(clusters,nchns);
    color_mat = colors(T,:);

    % create plot options
    clear plot_opts
    %     plot_opts.labels = Sz_preproc.segment_channel_labels{s_pat};
    plot_opts.labels = size(tbl.segment_data{1},1):-1:1;
    plot_opts.offset=600;
    plot_opts.plot_labels = true;
    plot_opts.linewidth = 0.5;
    plot_opts.clrs = color_mat;

    % plot time series
    ax1 = subplot(4,1,1:3);
    str = 'clustering by dominant frequencies';
    plot_sz(tbl,'fig_ind',0.5,'pltopts',plot_opts,'str',str);
    xlim([-10,tbl.duration(1)+10])

    % plot dominant freq   
    ax2 = subplot(4,1,4);
    hold on
    for c = 1:nchns
        plot(t,df_data(c,:),'Color',color_mat(c,:))
    end
    xline(0,'LineWidth',2)
    xline(tbl.duration(1), 'LineWidth',2)
    xlabel('Time (s)')
    ylabel('dominant frequency (Hz)')
    
    
%     % plot time series (experimental, didnt debug this)
%     ax2 = subplot(2,1,2);
%     plot(t,[glob_sync_allchns; cluster_abs_sync_ts]')
%     str = 'dfs for allchns and just of the clusters';
%     standard_title(tbl(s,:),'str',str);
%     legend_cell = cell(no_cl,1);
%     for i = 1:no_cl
%         legend_cell(i) = {['cl' num2str(i)]};
%     end
%     legend([{'dfs_allchns'};legend_cell],'Interpreter','none')

    linkaxes([ax1,ax2],'x')
    
end


%% internal functions
function [D,mat] = distfct_a(x,y,a)
% calculates the distance between 1-by-nwndws frequency vector x and all
% frequency vectors in rows of n-by-nwndws matrix y. Also works for scalar
% values

    mat = (1 - (x > 0 & y > 0) .* exp(-a*(x-y).^2)); 
    D = sum(mat,2);

%     mat = (1 - (x > 0 & y > 0) .* exp(-a*(log(1+x)-log(1+y)).^2)); 
%     mat(isnan(mat)) = 0;
%     D = sum(mat,2, 'omitnan');
%     D(isnan(D)) = 0;
end


% copied from calc.ms, execpt for not averaging dist_arr!
function out = ms_simfreq_noaverg(freq)
% calculates sync in between all freq-time graphs, i.e. calculates sync
% sync value of all the provided channels in freq for each time point
    nchns_inp = size(freq,1);
    out = zeros(1,size(freq,2));
    for i = 1:nchns_inp-1
        [~,mat] = distfct(freq(i,:), freq(i+1:end,:));
        out = out + size(mat,1)- sum(mat,1);
    end 
end

% transforms cell cluster representation into array
function T = cell2array_cl(cl,nchns)
    T = nan(nchns,1);
    for c=1:numel(cl)
        T(cl{c}) = c;
    end
end
end