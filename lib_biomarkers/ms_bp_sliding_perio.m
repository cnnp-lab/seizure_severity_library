function tbl = ms_bp_sliding_perio(meta_tbl,bp_tbl,opts)
    % aims to detect regular bandpower fluctuations in all channels 
    %
    % takes bandpower matrix of a seizure, sums it up along the channels,
    % computes windowed periodograms of this 1D function, checks the relative
    % power contained between 1.5-2 Hz with a half-mode-estimator in each 
    % window and extracts the maximum value of all windows
    %
    % inputs:
    %   - Sz: seizure table, provide segment data for nicer plots
    %   - bp_tbl: table with columns segment_id,bp,t_wndw
    %   - opts.fig_ind: -1 no plots, 0 new figure, >0 figure with that
    %       index
    % 
    % potential todos:
    %   - are we trying to detect burst-suppression patterns?
    %   - high amplitude channel dominate the resut. Normalize amplitude or similar?

    arguments
        meta_tbl
        bp_tbl
        
        opts.fig_ind = 0
    end
    
    % parse params
    bp_all = bp_tbl.bp;
    t_wndw_all = bp_tbl.t_wndw;
    nsz = size(meta_tbl,1);
    if size(meta_tbl,1) ~= size(bp_tbl,1)
        error('Segment and calculation table must have same number of rows.')
    end
    
    % calculation loop
    val = zeros(nsz,1);
    wndw = zeros(nsz,2);
    for s = 1:nsz
        bp = sum(bp_all{s},1);
        t_wndw = t_wndw_all{s};
    
        % set params
        t_start = 0;
        t_end = meta_tbl.duration(1);
        fs = floor(( (t_wndw(end) - t_wndw(1)) / size(bp,2) )^(-1));
        % ind = 0 < t_wndw & t_wndw < Sz_preproc.duration(s);
        ind = time2wnd_ind(t_wndw,t_start):time2wnd_ind(t_wndw,t_end);
        bp_ictal = bp(ind);
        t_ictal = t_wndw(ind);

        % sliding window periodogram
        freq_pw = linspace(0.5,4,30);
        if meta_tbl.duration(1) < 4
            val = 0;
        else
            if meta_tbl.duration(1) <= 8
                block_ind = [1 length(bp_ictal)];
                t_val = 0;
            else
                [block_ind, ~, t_val] = sliding_window(length(bp_ictal), 8 * fs, 'wndw_shift',floor(1*fs),...
                    'pre',0,'fs',fs);
            end

            nwndws = length(t_val);
            pxx_all = cell(nwndws,1);
            for w = 1:nwndws
                [pxx,f] = periodogram(bp_ictal(block_ind(w,1):block_ind(w,2)),[],freq_pw,fs);
                pxx_all{w} = pxx;
            end
            pxx_all = cell2mat(pxx_all);
        %     [pval, freq] = pwelch(bp(ind),24,4,freq_pw,fs);


            % calc measure
            f_bnd = [1.5 2.5];
            ind_ms = f_bnd(1) < f & f<f_bnd(2);

            % * v1: rel power in f_bnd
            % rel_pw_wndw = sum(pxx_all(:,ind_ms),2) ./ sum(pxx_all,2);
            % [val,ind_max] = max(rel_pw_wndw);
            % val
            % ff = [t_val', rel_pw_wndw];


            % * v2: power in half maximum thing where spike in f_bnd
            [val_all,ind_max_all] = max(pxx_all(:,ind_ms),[],2);
            strength = zeros(nwndws,1);
            for w = 1:nwndws

                % set local variables
                data = pxx_all(w,ind_ms);
                val = val_all(w);
                ind_max = ind_max_all(w);

                % calc half mode estimator
                ind_lwr = find(data(1:ind_max) < val/2,1,'last')+1;
                ind_hgr = ind_max + find(data(ind_max+1:end) < val/2,1,'first')-1;
                if isempty(ind_lwr)
                    ind_lwr = 1;
                end
                if isempty(ind_hgr)
                    ind_hgr = size(data,2);
                end

                % calc relative area under half mode estimator
                strength(w) = sum(data(ind_lwr:ind_hgr)) / sum(pxx_all(w,:));
            end
            [val_tmp,ind_max] = max(strength);
            val(s) = val_tmp;
            
            % find time window
            ind = block_ind(ind_max,1):block_ind(ind_max,2);
            wndw(s,:) = [t_ictal(ind(1)) t_ictal(ind(end))];
        end
    end
        
    % output table
    tbl = table(meta_tbl.segment_id,val',wndw,'VariableNames',...
        {'segment_id','ms_slidperio','ms_sp_timepoint'});
    
    % plotting
    if nsz == 1 && opts.fig_ind >= 0
        
        if ~ismember('segment_data',meta_tbl.Properties.VariableNames)
            warning('Provide segment data for analysis plots')
        else
            % create figure
            if opts.fig_ind == 0
                figure()
            else
                figure(opts.fig_ind)
            end
            
            % actual plotting
            ax1 = subplot(5,1,1);
            plot_sz(meta_tbl(s,:),'fig_ind',0.5)
            xlabel('time (s)')

            ax2 = subplot(512);
            plot(t_wndw,bp)
            % xlim([-120,370])
            xlabel('time (s)')
            title('summed bandpower')

            ax3 = subplot(513);
            plot_sz(meta_tbl(s,:),'fig_ind',0.5)
            title('time series of chosen window')
            ind = block_ind(ind_max,1):block_ind(ind_max,2);
            xlim([t_ictal(ind(1)) t_ictal(ind(end))])
            xlabel('time (s)')

            ax4 = subplot(514);
            plot(t_ictal(ind),bp_ictal(ind));
            title('bandpower of chosen window')
            xlim([t_ictal(ind(1)) t_ictal(ind(end))])
            xlabel('time (s)')

            ax5 = subplot(515);
            plot(freq_pw,pxx_all(ind_max,:))
            title('psd of bandpower of chosen window')
            xlabel('frequency (Hz)')
            
            linkaxes([ax1,ax2],'x')
            linkaxes([ax3,ax4],'x')
        end
    end
    
    
    
end