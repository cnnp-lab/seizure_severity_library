classdef DF_db < Abstract_calc_db
% This database aims to detect SWD. It recreastes the method used in "Automatic 
% Detection of Spike and Wave Discharges in the EEG of Genetic Absence Epilepsy 
% Rats From Strasbourg" from Hese et al. 2009 (dx.doi.org/10.1109/TBME.2008.2008858)
%
% for more information about the calculation performed see method
% detect_DF_chn below
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params()
    
    % name of calculation performed by database
    name_calc = 'df';     
    
    % parameters 
    default_params.wndw_len = 256 * 2.5;        % suggestion from paper: 1/4 fs
    default_params.wndw_overlap = 256*2.5 - 64; % suggestion from paper: floor(fs / opts.freq_fund(1))
    
    % frequency params
    default_params.freq_fund = [1 12.5];        % boundaries for df
    default_params.freq_cutoff = [0.5 49];     % range that harmonics are searched in
    default_params.freq_step = 0.1;
    
    % harmonic analysis params (best not change)
    default_params.beta_s = 8;
    default_params.alpha = 0.05;
    default_params.delta_bnd = [0.1,0.5];
    
    % classification and gamma
    default_params.gamma = 0;       % whether to include factor for accounting for sleep spindles or not
    default_params.class = 1;       % whether to classify each time point as SWD (points where not SWD classified automatically have df 0
    
    % window for classification
    default_params.class_buffer = 1;

    % example version
    vers_id = 'v2';
    vers_cmmt = 'classifies only time series from seizure start to finish, thrhld algo on normalized E_final';
    
end

function tbl = calc_seg(data_tbl,params,fig_hndl) 

    % calcuation
    nchns = size(data_tbl.segment_data{1},1);
    df_chns = cell(nchns,1);
    for c = 1:nchns
        res =  DF_db.calc_chn(data_tbl,params,c,[]);
        df_chns{c} = res.df;
    end
    
    % create output tbl
    df_mat = cell2mat(df_chns);
    tbl = table({df_mat},{res.t_wndw},'VariableNames',{'df','t_wndw'});
    
    % Your plotting here (instead of the warnings if available)
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        % plot segment analysis
        ax1 = subplot(4,1,1:2);
        plot_sz(data_tbl,'str','time series and dominant frequencies','fig_ind',0.5);
        ax2 = subplot(4,1,3:4);
        plot(res.t_wndw, df_mat)
        ylabel('dom. freq. (Hz)')
        xlabel('Time (s)')
        linkaxes([ax1 ax2],'x')
        xlim([res.t_wndw(1) res.t_wndw(end)])
    end
    
end

function tbl = calc_chn(data_tbl,params,chn,fig_hndl) 
    % inputs:
    %   - data_tbl: table with one row for the channel data
    %   - params: struct with all parameters in fields. For more info about
    %       params and their default values look at get_algo_params
    %   - plot_flag: bool. If true, plots result plot into current plot
    % outputs:
    %   - tbl_out: table with values for dominant freq

    % NOTES: 
    % - no classification possible for first half window_len at data begin,
    % no classification possible at end for cutoff because of uneven shifts

    % ---------------------------
    % STEP 0: parameters and data preparation
    % ---------------------------
    s = 1;
    if ~istable(data_tbl) || size(data_tbl,1) ~= 1
        error('Data input has to be table with exactly one row')
    end
    wndw_overlap = params.wndw_overlap;
    wndw_len = params.wndw_len;
    freq_fund = params.freq_fund;
    freq_cutoff = params.freq_cutoff;
    freq_step = params.freq_step;
    beta_s = params.beta_s;
    alpha = params.alpha;
    delta_bnd = params.delta_bnd;
    data = data_tbl.segment_data{s};
    wndw_shift = wndw_len - wndw_overlap;

    % t_eeg = (1:length(data))/data_tbl.segment_fs - data_tbl.segment_pre(s);
    % ictal_ind = 0 <= t_eeg & t_eeg < data_tbl.duration(s);
    % data = data(chn,ictal_ind);

    data = data(chn,:);
    fs = data_tbl.segment_fs(s);

    % notation as in paper
    L = wndw_len;
    N = L; % size of N makes no difference as constant in all E_prelim
    S = wndw_shift;

    % assert inputs
    if floor(freq_cutoff(1) / freq_step) * freq_step ~= freq_cutoff(1)
        error('Low frequency cutoff needs to be devisible by frequency step!');
    end
    if freq_cutoff(2) > 2*fs
        freq_cutoff(2) = floor(fs/2);
        warning(['freq_cutoff lowered to ' num2str(freq_cutoff(2)) ' because of low sample frequency'])
    end
    % if freq_fund(2) * 2 > freq_cutoff(2)
    %     error('Frequency cutoff must be at least twice fundamental frequency')
    % end
    if beta_s < 0
        error('Beta_s needs to be larger than 0!');
    end
    if floor(L/2) * 2 ~= L
        error('Window length must be even!')
    end
    if floor(S/2) * 2 ~= S
        error('Windows step must be even!')
    end


    % ---------------------------
    % STEP 1: calculate spectrogram
    % ---------------------------
    % presentation step 2

    freq = freq_cutoff(1):freq_step:freq_cutoff(2);
    wndw_overlap = wndw_len - wndw_shift;
    [stft,~,t_wndw] = spectrogram(data, hamming(wndw_len), wndw_overlap, freq,fs,'yaxis'); 
        % zero padding performed in matlabs spectrogram
    spectr = abs(stft).^2;
    ff = [t_wndw - data_tbl.segment_pre(s);spectr];   % for debugging


    % ---------------------------
    % STEP 2: estimate background spectrum
    % ---------------------------
    % presentation step 3.1
    spectr_bg = median(spectr,2)/log(2);


    % ---------------------------
    % STEP 2.5: confine window for classification to some time interval
    % ---------------------------
    % params
    start_class_s = 0;                        % in seconds
    end_class_s = data_tbl.duration(s);     % in seconds
    wndw_s = wndw_len / fs;
    buffer = params.class_buffer;     % buffer around the specified time

    % calc ind
    t_eeg = (1:size(data_tbl.segment_data{s},2))/data_tbl.segment_fs(s) - data_tbl.segment_pre(s);
    t_wndw = t_wndw - data_tbl.segment_pre(s);
    % ind_eeg = -buffer +(-wndw_s/2 +start_s) <= t_eeg & t_eeg < (end_s+wndw_s/2) +buffer;    % wndw_s to ensure that last wndw totally in data
    % ind_wndw = -buffer +(-wndw_s/2 +start_s) <= t_wndw & t_wndw < (end_s) +buffer;         % buffer to ensure that classification for all ictal time points
    ind_wndw = -buffer +(-wndw_s/2 +start_class_s) <= t_wndw & t_wndw < (end_class_s +wndw_s/2) +buffer;         % buffer to ensure that classification for all ictal time points
    t_wndw = t_wndw(ind_wndw);
    ind_eeg = t_wndw(1)-wndw_s/2 <= t_eeg & t_eeg < t_wndw(end) + wndw_s/2;
    t_eeg = t_eeg(ind_eeg);


        % stuff in brackets ensures that whole specified period classified,
        % then add the buffer in both directions
    spectr = spectr(:,ind_wndw);
    data = data(:,ind_eeg);


    % ---------------------------
    % STEP 3: Harmonic analysis
    % ---------------------------
    % presentation step 3.2
    % calc spectrogram prime 
    sz = size(spectr);
    [f_ind, w_ind] = find(spectr >= spectr_bg);    
        % idea: f_ind and w_ind are indices of potential candidates that get
        % reduced until only very last one left

    % % contained in paper: only include points with local maximum in frequency
    % less_l = spectr(sub2ind(sz, max(1,f_ind-1), w_ind)) <= spectr(sub2ind(sz,f_ind, w_ind));
    % less_r = spectr(sub2ind(sz,f_ind,w_ind)) >= spectr(sub2ind(sz, min(f_ind+1,sz(1)), w_ind));
    % ind = less_l & less_r; % determine if local max
    % f_ind = f_ind(ind);
    % w_ind = w_ind(ind);

    spectr_prime = sparse(f_ind,w_ind,spectr(sub2ind(sz,f_ind,w_ind)),size(spectr,1),size(spectr,2));


    % extract reasonable candidate frequencies
    % * confine to fundamental freq band
    f_fund_low_ind = (sum(freq < freq_fund(1))+1);
    f_fund_high_ind = (sum(freq < freq_fund(2)));
    ind = (f_fund_low_ind <= f_ind) & (f_ind <= f_fund_high_ind); 
    f_ind = f_ind(ind);
    w_ind = w_ind(ind);


    % * calculate combs (mainly for E_prelim later) (debugcheck: passed)
    f_ind_uniq = unique(f_ind);
    max_no_combs = ceil(freq_cutoff(2) / freq_fund(1));
    no_combs = zeros(max(f_ind_uniq),1);
    combs_logind = cell(max(f_ind_uniq),max_no_combs);
        % stores the combs as logical pattern by frequency_index and comb_index

    for ind = 1:length(f_ind_uniq)

        f_ind_curr = f_ind_uniq(ind);
        freq_curr = freq(f_ind_curr);

        harm = freq_curr:freq_curr:freq_cutoff(2);
        no_combs(f_ind_curr) = length(harm);

        harm_ind_lwr = time2wnd_ind(freq, harm * (1-alpha),'op','less')+1;
        harm_ind_hgr = time2wnd_ind(freq, harm * (1+alpha),'op','lesseq');

        for i = 1:length(harm_ind_hgr)
            harm_ind = zeros(size(freq));
            harm_ind(harm_ind_lwr(i) : harm_ind_hgr(i)) = 1;
            combs_logind{f_ind_curr,i} = logical(harm_ind);
        end

    end


    % * check the 1st and 2nd harmonic conditions
    ind = false(size(f_ind));
    for cand = 1:length(f_ind)

        f = f_ind(cand);
        w = w_ind(cand);

        ind(cand) = ...
            any(spectr_prime(combs_logind{f,1},w) > beta_s * spectr_bg(combs_logind{f,1})) || ...
            any(spectr_prime(combs_logind{f,2},w) > beta_s * spectr_bg(combs_logind{f,2}));

    end
    f_ind = f_ind(ind);
    w_ind = w_ind(ind);


    % calc prelim evidence of reasonable candidate freq (TODO: error warning if no candidates)
    % presentation step 4
    % * main calc loop of E_prelim
    V_w = sum(hamming(wndw_len).^2)/N;
    prelim_val = zeros(length(f_ind),1);
    prelim_val_nogmm = zeros(length(f_ind),1);
    for cand = 1:length(f_ind)  

        % set variables
        f_ind_curr = f_ind(cand);
        w_ind_curr = w_ind(cand);
        no_combs_curr = no_combs(f_ind_curr);

        % calc comb values
        h_m0_all = nan(no_combs_curr,1);
        for comb = 1:no_combs_curr
            h_m0_all(comb) = sum(spectr_prime(combs_logind{f_ind_curr,comb},w_ind_curr));
        end
        sum_h_m0 = sum(h_m0_all);

        % prelim evidence
        psi = @(x) (x <= 0.1) .* x/0.1 + (x>0.1);
        phi = @(x) (x <= 0.7) + (~(x<=0.7)) .* (x-1)/(0.7-1);
        gamma = psi(h_m0_all(1) / sum_h_m0) * prod(phi(h_m0_all / sum_h_m0));
        V_xk = N/L/V_w * sum(spectr(:,w_ind_curr));

        % calc preliminary evidence
        if params.gamma
            prelim_val(cand) = gamma * sum_h_m0 / V_xk;
        else
            prelim_val(cand) = sum_h_m0 / V_xk;
        end
    end
    E_prelim = zeros(size(spectr)); % needs to be full for convolution below anyways
    E_prelim(sub2ind(size(E_prelim),f_ind,w_ind)) = prelim_val;


    % * smooth E_prelim (not same smooth as in presentation! Maybe test to smooth
    %   here instead of in the very end)
    weight = [1 4 6 4 1] / 16;
    cutoff_weight = floor(length(weight)/2);
    for f_ind_curr = 1:size(E_prelim,1)
        tmp = conv(weight, E_prelim(f_ind_curr,:));
        E_prelim(f_ind_curr,:) = tmp(1+cutoff_weight:end-cutoff_weight);
    end


    % calc final candidate (presentation step 5.1 & 6)
    [E_final, ind_E_final] = max(E_prelim,[],1);
    dom_freq = freq(ind_E_final) .* (E_final > 0);

    % ---------------------------
    % STEP 4: Classification
    % ---------------------------
    % presentation step 5.2
    % thrhld determined on with normed E_final, but after that has to be at
    % least 0.1 in actual scale. Clustering then performed with this final
    % threshold

    if params.class
        E_final_norm = E_final / max(E_final);
        % E_final_norm = E_final_gmm;

        % generate threshhold data
        thrhld_cand = linspace(delta_bnd(1),delta_bnd(2),100);
        detected = nan(length(thrhld_cand),1);
        for t=1:length(thrhld_cand)
            detected(t) = sum(...
                classify_data(E_final_norm, thrhld_cand(t), L, S, fs, size(data,2))) ...
                / size(data,2);
        end

        % find plateau with HSM
        left = 1;
        right = length(detected);
        pointer = left;
        thrhld_scaled = -1;
        while thrhld_scaled == -1
            if right - left == 2
                thrhld_scaled = thrhld_cand(left+1);
            elseif right - left == 1
                thrhld_scaled = thrhld_cand(left);
            else
                dist_ind = floor((right - left) / 2);
                min_dist = Inf;
                for j = 0:ceil((right-left)/2)-1
                    d = detected(left+j) - detected(left+j + dist_ind);
                    if d < min_dist
                        min_dist = d;
                        pointer = j;
                    end
                end
                right = left + pointer + dist_ind;
                left = left + pointer;         
        %         disp(['left: ' num2str(thrhld_cand(left)) ', right: ' num2str(thrhld_cand(right))])
            end
        end

        thrhld_scaled = min(1, max( 0.1/max(E_final), thrhld_scaled));
        % thrhld within 0.1 and 1
        [class_time,class_wndw] = classify_data(E_final_norm, thrhld_scaled, L, S, fs, size(data,2));
        thrhld = thrhld_scaled * max(E_final);

        % drop dominant frequencies where time points where not classified as SWD
        % 1st half of presentation step 7
        dom_freq = dom_freq .* class_wndw;

    else
        thrhld = 1;
    end

    % ---------------------------
    % generate correct outputs
    % ---------------------------
    % smooth final frequencies
    % 2nd half of presentation step 7
    mm_half_wndw = floor(wndw_len/2 / wndw_shift);
    dom_freq = movmedian(dom_freq,[mm_half_wndw,mm_half_wndw]);

    % make nice output table
    tbl = table(chn,dom_freq,t_wndw,'VariableNames',{'chn','df','t_wndw'});

    % ---------------------------
    % plot (intermediate-) results
    % ---------------------------
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        set(groot, 'DefaultAxesTickLabelInterpreter', 'none');

        % feedback plotting
        tiledlayout(7,1)
        ax1 = nexttile;
            plot(t_eeg,data)
            str = ['SWD detection, chn ' num2str(chn)];
            freq = freq_cutoff;
            standard_title(data_tbl,'str',str,'freq',freq)
            hold on;
            a = area(t_eeg, (max(data) - min(data)) * class_time + min(data), 'basevalue',min(data),'FaceColor', [.5 .5 .5 ]);
            a.FaceAlpha = 0.35;
            hold off;
            ylabel('classification')

        ax2 = nexttile;
            imagesc(log(flip(spectr,1)),'XData',[t_eeg(1) t_eeg(end)], 'YData', [freq_cutoff(2) freq_cutoff(1)])
            ylabel('log(spectr)')
            colorbar

        ax5 = nexttile;
            imagesc(log(spectr_prime),'XData',[t_eeg(1) t_eeg(end)], 'YData',[freq_cutoff(1) freq_cutoff(2)])
            ylabel('log(S_prime)','Interpreter','none')
            colorbar

        ax6 = nexttile;
            freq = freq_cutoff(1):freq_step:freq_cutoff(2);
            f_ind_fund = freq_fund(1) <= freq & freq <= freq_fund(2);
            E_p_plot = E_prelim(f_ind_fund,:);
            imagesc(E_p_plot,'XData',[t_eeg(1) t_eeg(end)], 'YData',[freq_fund(1) freq_fund(2)])
            ylabel('E_prelim','Interpreter','none')
            colorbar

        ax3 = nexttile;
            bar(t_wndw,E_final );
            yline(thrhld)
            ylabel('E_final ','Interpreter','none')

        ax4 = nexttile;
            plot(t_wndw, dom_freq)
            ylabel('freq')

        ax7 = nexttile;
            plot(t_wndw,movmedian(dom_freq.*class_wndw,[mm_half_wndw,mm_half_wndw]));
            ylabel('freq movmedian')
            xlabel('Time (s)')


        linkaxes([ax1 ax2 ax3 ax4 ax5 ax6, ax7],'x')
        xlim([t_eeg(1) t_eeg(end)]);


        % further plot for classification analysis
        if params.class
            figure()
            tiledlayout(2,1)
            nexttile
                plot(thrhld_cand, detected)
                hold on
                xline(thrhld_scaled)
                hold off
                xlabel('Threshold')
                ylabel('percentage of windows classified as SWD')
                title('HSM analysis')
    
            nexttile
                plot(t_eeg,data)
                hold on;
                a = area(t_eeg, (max(data) - min(data)) * class_time + min(data), 'basevalue',min(data),'FaceColor', [.5 .5 .5 ]);
                a.FaceAlpha = 0.35;
                hold off;
        end



    end


    % local functions
    function [class_smpl, class_wndw] = classify_data(E_final, thrhld, L, S, fs, no_smpls)
    % returns the classification in form of 0-1 array given a certain thrhld,
    % E_final and thrhld need to be same scale (i.e. no normalization
    % performed)


        % classify segments
        class_ind = E_final > 0;
        isOne = false;
        hasLargerEvidence = false;
        ptr = -1;
        segments_wndw = nan(length(class_ind),2);
        for i = 1:length(class_ind)
             if isOne && ~class_ind(i)
                 if hasLargerEvidence 
                     segments_wndw(i-1,:) = [ptr i-1];
                 end
                 isOne = false;
             elseif ~isOne && class_ind(i)
                 isOne = true;
                 ptr = i;
                 hasLargerEvidence = E_final(i) >= thrhld;
             else
                hasLargerEvidence = hasLargerEvidence || E_final(i) >= thrhld;
             end
        end
        if isOne && class_ind(i) == 1
            if hasLargerEvidence
                segments_wndw(i,:) = [ptr i];
            end
        end
        segments_wndw = rmmissing(segments_wndw);   


        % transform segments to sample domain
    %     segments_smpl = [-S/2 + segments_wndw(:,1)*S + L/2, S/2 + 1 + segments_wndw(:,2)*S + L/2];
    %         % formular from paper
        segments_smpl = [(-ceil(S/2) + (segments_wndw(:,1)-1)*S + L/2), (ceil(S/2) + (segments_wndw(:,2)-1)*S + L/2)];

        % remove segments short than 1s
        for i=1:size(segments_smpl,1)
            if segments_smpl(i,2)-segments_smpl(i,1) < fs
                segments_smpl(i,:) = [nan nan];
                segments_wndw(i,:) = [nan nan];
            end
        end
        segments_smpl = rmmissing(segments_smpl);
        segments_wndw = rmmissing(segments_wndw);

        % merge segments with distance smaller than 1s
        for i=1:size(segments_smpl,1)-1
            if segments_smpl(i+1,1)-segments_smpl(i,2) < fs
                segments_smpl(i+1,1) = segments_smpl(i,1);
                segments_smpl(i,:) = [nan nan];
                segments_wndw(i+1,1) = segments_wndw(i,1);
                segments_wndw(i,:) = [nan nan];
            end
        end
        segments_smpl = rmmissing(segments_smpl);
        segments_wndw = rmmissing(segments_wndw);

        % create 0-1 array on sample domain
        class_smpl = zeros(1,no_smpls);
        class_wndw = zeros(size(E_final));
        nwndws = size(class_wndw,2);
        if ~isempty(segments_smpl)
            for i=1:size(segments_smpl,1)           
                class_smpl(segments_smpl(i,1) : min(no_smpls,segments_smpl(i,2))) = 1;
                class_wndw(segments_wndw(i,1) : min(nwndws,segments_wndw(i,2))) = 1;
            end
        end


end
end

% Don't change anything here
function [tbl,parset] = plot(data,varargin)
    % for documentation look at Abstract_calc_db.plot
    varargin = [varargin {'subclass_str'} {mfilename('class')}];
    [tbl,parset] = Abstract_calc_db.plot(data,varargin{:});
end

end
end

