function [data_out,t_wndw] = sync_corr(Sz_preproc,s,wndw_len, wndw_shift,opts)
% calculates the correlation between the channels in sliding window
% approach, returns sum of co_h leading eigenvalues per window together
% with time point the window is associated with

%% assert data
arguments
    Sz_preproc
        % needs fields: segment_data,segment_fs,segment_pre,chn_label,...
    s
    wndw_len
    wndw_shift
    
    opts.fig_ind = 1
    opts.co_l = 0.6
    opts.co_h = 0.05
    opts.offset = 300;
end
fig_ind = opts.fig_ind;
co_l = opts.co_l;
co_h = opts.co_h;


% standard
data = Sz_preproc.segment_data{s};
fs = Sz_preproc.segment_fs(s);
pre = Sz_preproc.segment_pre(s);

%% calc corr
% sliding wndw eigval
nchns = size(data,1);
[block_ind, ~, t_wndw] = sliding_window(size(data,2),wndw_len,...
    'wndw_shift',wndw_shift,'fs',fs,'pre',pre,'ass',0.5);   
    % 0.5 ass because imagesc plots symmetrically around the time point
nwndws = length(t_wndw);
eigvals = nan(nchns, nwndws);
for w = 1:nwndws
    data_wndw = data(:,block_ind(w,1) : block_ind(w,2));
    data_wndw = zscore(data_wndw,1,2);
    corr_mat = 1/wndw_len * (data_wndw * data_wndw');
    eigvals(:,w) = sort(eig(corr_mat),'descend');
end
eigvals_norm = eigvals / nchns; % because symm and trace is one

% sum leading eigvals
if co_h >= 1
    ind_up = 1:co_h;
else
    ind_up = max(1,floor(co_h*nchns));
end
ind_lw = min(nchns,nchns-floor(co_l*nchns)+1);
sum_eig_hgr = (sum(eigvals_norm(1:ind_up,:),1));
sum_eig_lwr = (sum(eigvals_norm(ind_lw:end,:),1));
% disp(['Corr: High eigvals from 1 to ' num2str(ind_up)]);
% disp(['Corr: Low eigvals from ' num2str(ind_lw) ' to ' num2str(nchns)]);
data_out = sum_eig_hgr;


    
%% plotting
if fig_ind > 0
    % plot time series
    t_eeg = (1:size(Sz_preproc.segment_data{s},2))/Sz_preproc.segment_fs(s) - Sz_preproc.segment_pre(s);

    figure(fig_ind)
    ax1 = subplot(4,1,1:2);
        standard_plot_timeSeries(Sz_preproc,s,'offset',opts.offset);
        standard_title(Sz_preproc,s,'str','correlation');

    % plot zscored eigvals
    ax2 = subplot(413);
        imagesc(zscore(eigvals_norm,1,2), 'XData', t_wndw)

    % plot synchrony function
    figure(fig_ind)
    ax3 = subplot(414);
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


end
