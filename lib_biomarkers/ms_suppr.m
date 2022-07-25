function [tbl,parset] = ms_suppr(meta_tbl,vals,t,opts)
% calculates PGES duration, partial suppression duration and its strength 
% based on the input data. Input data can be e.g. line length or some
% amplitude measure. Algo will determine extreme points of the input data 
% against preictal period and then calculate the output values based on the
% occurence of extreme points in the post ictal period.
% 
% inputs:
%   - meta_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - val_tbl: cell array with data matrix chns-by-time for each segment
%   - t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below

    arguments
        meta_tbl
        vals
        t

        % determine extreme points
        opts.use_z = 0          % 1: zscoring, 0: percentile
        opts.thrhld_prctile = 0.05
        opts.z_trhld = -3
        
        % output creatoin
        opts.buffer_post_s = 1 %1sec buffer after seizure to ignore
        opts.pges_th = 0.8
        opts.partialsupp_th = 0.1
        opts.consec_len_s = 2.5
        
        % plotting
        opts.fig_ind = 0            % -1: no plotting, 0:new figure, k>0: figure with index k
        
    end
    buffer_post_s = opts.buffer_post_s;
    thrhld_prctile = opts.thrhld_prctile;
    consec_len_s = opts.consec_len_s;
    pges_th = opts.pges_th;
    partialsupp_th= opts.partialsupp_th;
    fig_ind = opts.fig_ind;
    z_trhld = opts.z_trhld;
    if size(meta_tbl,1) ~= numel(vals) || size(meta_tbl,1) ~= numel(t) || ...
        numel(vals) ~= numel(t)
        error('All inputs must have same amount of rows.')
    end
    
    nsegs = size(meta_tbl,1);
    tbl = cell(1,nsegs);
    for s = 1:nsegs

        % specify variables for algorithm
        data_wndw = vals{s};
        t_wndw = t{s};
        ind_pre = 1:time2wnd_ind(t_wndw,0);
        data_pre = data_wndw(:,ind_pre);
        nchns = size(data_wndw,1);

        % determine extreme values
        if opts.use_z
            % v2: zscoring with mad
            ind_mat = zscore_vs_preictal(log(1+data_wndw), ind_pre(end),'is_mad',1) < z_trhld;
        else
            % v1: percentile approach
            thrhld = prctile(data_pre,thrhld_prctile*100, 2);
            ind_mat = data_wndw < thrhld;
        end
        
        
        % final class
        wndw_len_s = t_wndw(2)-t_wndw(1);
        consec_len_w = ceil(consec_len_s / wndw_len_s);
        d = ind_mat;
        tl=time2wnd_ind(t_wndw,meta_tbl.duration(s)+buffer_post_s);
        data_post = d(:,tl:end);
        time_post = t_wndw(tl:end);
        data_consec=movsum(data_post,consec_len_w,2);
        %plot(sum(ms==recruitment_threshold))
%         imprint=(ms==recruitment_threshold);

        
        supp_ind = data_consec > 0;
        
        supp_pct=sum(supp_ind,1)/nchns;
        
        PGES_ind=supp_pct>pges_th;
        PGES_dur_s=sum(PGES_ind)*wndw_len_s;%assumes no overlap in windows, which is true in current implementation.
        
        partialsupp_ind=supp_pct<=pges_th & supp_pct>partialsupp_th;
        part_suppr_dur_s=sum(partialsupp_ind)*wndw_len_s;%assumes no overlap in windows, which is true in current implementation.
        part_suppr_val=median(supp_pct(partialsupp_ind));
        
        % create output values
        tbl{s} = table(meta_tbl.segment_id(s),PGES_dur_s, part_suppr_dur_s, part_suppr_val,'VariableNames',...
            {'segment_id','PGES_dur','part_suppr_dur','suppr_strength'});
        parset = opts;
        
        %% old
% 
%         dur = zeros(nchns,1);
%         for c = 1:nchns
% 
%             % find ind where PGES ends
%             ind_PGES_end = find(~supp_ind(c,:),1) - 1;
%             if isempty(ind_PGES_end)
%                 ind_PGES_end = size(data_consec,2);
%             end
% 
%             % save variables
%             dur(c) = (ind_PGES_end)*wndw_len_s;
%             ind_dur_mat(c,ind_predurbuffer+1:ind_predurbuffer + ind_PGES_end) = 1;
% 
%         end
% 
%         % calc meta_tbl duration
%         PGES_dur_s = (sum(mean(ind_dur_mat,1) > mean_th)) * wndw_len_s;
% 
%         % calc partial outputs
%         yval_part_suppr = mean(ind_dur_mat,1);
%         ind_post = time2wnd_ind(t_wndw,meta_tbl.duration(s)+buffer_post_s):length(t_wndw);
%         % [rel_suppr, ind_max] = max(yval_rel_suppr(ind_post)); 
%             % ind_max will always be first post time point (after buffer) so no
%             % need to calc!
%         part_suppr_val = yval_part_suppr(ind_post(1));
%         part_suppr_dur = find(yval_part_suppr(ind_post) <= part_suppr_val /2,1) - 1;
%         if isempty(part_suppr_dur)
%             part_suppr_dur = numel(ind_post);
%         end
%         part_suppr_dur_s = (part_suppr_dur) *wndw_len_s;


%         % create output values
%         tbl{s} = table(meta_tbl.segment_id(s),PGES_dur_s, part_suppr_dur_s, part_suppr_val,'VariableNames',...
%             {'segment_id','PGES_dur','part_suppr_dur','suppr_strength'});
%         parset = opts;
    end
    tbl = vertcat(tbl{:});
    
    % plotting
    if nsegs == 1 && fig_ind > -1
        if fig_ind == 0
            figure()
        else
            figure(fig_ind)
        end
        ax1 = subplot(12,1,1:5);
        plot_sz(meta_tbl(s,:),'fig_ind',0.5)
        standard_title(meta_tbl(s,:),'str', 'amplitude suppression')

        ax2 = subplot(12,1,6:9);
        imagesc(supp_ind,'XData',time_post)
        ylabel('low amplitude')


        ax4 = subplot(12,1,10:12);
        plot(time_post, supp_pct,'-k');
        yline(pges_th,'r');
        yline(partialsupp_th,'b');
        hold on;
        a = area(time_post, supp_pct > pges_th, 'FaceColor', [.5 .5 .5 ]);
        a.FaceAlpha = 0.35;
        hold off;
        ylabel('PGES meta_tbl')
        xlabel('Time (s)')

        linkaxes([ax1,ax2,ax4],'x')
        %     linkaxes([ax3,ax2], 'y')
        xlim([t_wndw(1) t_wndw(end)])

    end
    
end