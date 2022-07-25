function [mpc_vals, t_wndw] = sync_plv(Sz_preproc,s,wndw_len,wndw_shift,opts)
% plv = phase lock value

arguments
    Sz_preproc
    s
    wndw_len 
    wndw_shift 
    
    opts.fig_ind = 6
end



%% calc mean phase coherence (mpc)
tic 
data = Sz_preproc.segment_data{s};
fs = Sz_preproc.segment_fs(s);
pre = Sz_preproc.segment_pre(s);
[block_ind,~,t_wndw] = sliding_window(size(data,2), wndw_len, ...
    'wndw_overlap',wndw_len - wndw_shift,'pre',pre,'fs',fs);
mpc_vals = zeros(size(t_wndw));
for w = 1:length(t_wndw)
    disp([num2str(w) '/' num2str(length(t_wndw))])

    % calc mean phase coherence matrix
    nchns = size(data,1);
    mpc_mat = zeros(nchns);
    for i = 1:nchns-1
        for j = i+1:nchns
%                 disp([num2str(i) ',' num2str(j)]);
            mpc_mat(i,j) = calc_mpc_R(data([i j],block_ind(w,1):block_ind(w,2)));
        end
    end

    % compress matrix into one value
    mpc_vals(w) = sum(sum(mpc_mat)) / (nchns*(nchns-1)/2);
end
toc

%% plot
figure(opts.fig_ind)
t_eeg = (1:size(Sz_preproc.segment_data{s},2))/Sz_preproc.segment_fs(s) - Sz_preproc.segment_pre(s);
ax1 = subplot(2,1,1);
    standard_plot_timeSeries(Sz_preproc,s)
    standard_title(Sz_preproc,s,'str','correlation, leading eigenval')
ax2 = subplot(2,1,2);
    plot(t_eeg,interp1(t_wndw,mpc_vals,t_eeg))
    %     hold on; yline(3,'r'); hold off;
    %     hold on; yline(-3,'r'); hold off;


% synchronize axes
linkaxes([ax1 ax2],'x')
xlim([t_eeg(1) t_eeg(end)])





% ------------------------
function R = calc_mpc_R(data)
    if size(data,1) ~= 2
        error('Data needs to be 2-by-something!');
    end

    analy_sign = hilbert(data')';
%     real_sign = real(analy_sign);
%     stilde = imag(analy_sign);
%     rel_phase = atan((stilde(1,:).*real_sign(2,:) - real_sign(1,:).*stilde(2,:)) ./ ...
%         (real_sign(1,:).*real_sign(2,:) - stilde(1,:).*stilde(2,:)));
    phase = atan(imag(analy_sign) ./ real(analy_sign));
    rel_phase = phase(1,:) - phase(2,:);
    R = 1/size(data,2) * abs(sum(exp(1i*(rel_phase))));
end
    


%% 
% analy_sign = hilbert(data')';
% hilbert_phase = atan(imag(analy_sign) ./ real(analy_sign));
% 
% s = real(analy_sign);
% stilde = imag(analy_sign);
% alt_ind = atan((stilde(1,:).*s(2,:) - s(1,:).*stilde(2,:)) ./ ...
%     (s(1,:).*s(2,:) - stilde(1,:).*stilde(2,:)));
% 
% rph_h = hilbert_phase(1,:) - hilbert_phase(2,:);
% rph_t = alt_ind;
% rel_phase = rph_t;

end