function [data_out, t_wndw] = sync_coh(Sz_preproc,s,wndw_len,wndw_shift,opts)
% calculates the coherence between all channels in sliding window
% approach, returns sum of co_h leading eigenvalues per window together
% with time point the window is associated with

%% args processing
arguments
    Sz_preproc
    s
    wndw_len
    wndw_shift
    
    opts.fig_ind = 2
    opts.coh = -1
    opts.freq = -1
    
    opts.co_h = 0.05;
    opts.co_l = 0.60;
end
fig_ind = opts.fig_ind;
if opts.freq == -1
    freq = [0.5,120];
end

fs = Sz_preproc.segment_fs(s);
pre = Sz_preproc.segment_pre(s);

%% calculation 
% calc coherence if not given
if opts.coh == -1
    coh = bivarsw_coherence(Sz_preproc.segment_data{s}, Sz_preproc.segment_fs(s), ...
        wndw_len, wndw_len - wndw_shift, freq);
else
    warning('Coherence data provided. Ensure that window length and overlap correctly passed to function!')
    coh = opts.coh;
end

% calc eigenvalues:
nchns = size(Sz_preproc.segment_data{s},1);
nwndws = size(coh,2);
eigvals = nan(nchns, nwndws);
for w = 1:nwndws
    coh_mat = squareform(coh(:,w)) + eye(nchns);   % adding eye or not makes no difference (just makes eigvals one larger)
    eigvals(:,w) = sort(eig(coh_mat),'descend');
end
eigvals_norm = eigvals./sum(eigvals,1);

% sum leading eigvals
if opts.co_h > 1
    ind_up = 1:opts.co_h;
else
    ind_up = max(1,floor(opts.co_h*nchns));
end
ind_lw = min(nchns,nchns-floor(opts.co_l*nchns)+1);
sum_eig_hgr = sum(eigvals_norm(1:ind_up,:),1);
sum_eig_lwr = sum(eigvals_norm(ind_lw:end,:),1);
% disp(['Coh: High eigvals from 1 to ' num2str(ind_up)]);
% disp(['Coh: Low eigvals from ' num2str(ind_lw) ' to ' num2str(nchns)]);
data_out = sum_eig_hgr;


%% plotting
[~,~,t_wndw] = sliding_window(size(Sz_preproc.segment_data{s},2),...
    wndw_len,'wndw_overlap',wndw_len-wndw_shift,'fs',fs,'pre',pre,'ass',0.5);
if fig_ind > 0
    % plot time series
    figure(fig_ind)
    t_eeg = (1:size(Sz_preproc.segment_data{s},2))/Sz_preproc.segment_fs(s) - Sz_preproc.segment_pre(s);
    ax1 = subplot(4,1,1:2);
    standard_plot_timeSeries(Sz_preproc,s);
    standard_title(Sz_preproc,s,'str','coherence');

    % plot zscored eigvals
    ax3 = subplot(413);
        imagesc(zscore(eigvals_norm,1,2), 'XData', t_wndw)

    % plot synchrony function
    figure(fig_ind)
    ax2 = subplot(414);
        plot(t_eeg,interp1(t_wndw,sum_eig_hgr,t_eeg))
        hold on;
        plot(t_eeg,interp1(t_wndw,sum_eig_lwr,t_eeg))
        hold off;
        legend('hgr','lwr')


    % synchronize axes
    linkaxes([ax1 ax2 ax3],'x')
    figure(fig_ind)
    xlim([t_eeg(1) t_eeg(end)])
end