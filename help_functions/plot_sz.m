function plot_sz(sz_data, opts)
arguments
    sz_data                 % table with one row with segment data
    
    opts.fig_ind = 0        % figure to plot in, 0 new figure, 0.5 current
    opts.pltopts            % plot opts as in vis_plot_eeg
    opts.plt_numbers = 1;   % plot numbers instead of channel labels
    opts.str = ''           % additional string for title
end
    
    % input handling
    if ~istable(sz_data) || size(sz_data,1) ~= 1
        error('Data has to be table with exaclty 1 row.')
    end
    
    % prepare params
    s = 1;
    if ~isfield(opts,'pltopts')
        pltopts.offset=800;
        pltopts.plot_labels = true;
        if opts.plt_numbers
            pltopts.labels = size(sz_data.segment_data{s},1):-1:1;
        else
            pltopts.labels = sz_data.segment_channel_labels{s};
        end
    else
        pltopts = opts.pltopts;
    end

    % figure
    if opts.fig_ind > 0.5
        figure(opts.fig_ind);
    elseif opts.fig_ind == 0 
        figure();
    end
    
    % actual plotting
    eeg=sz_data.segment_data{s}';
    fs = sz_data.segment_fs(s);
    t=[1:max(size(eeg))]/fs - sz_data.segment_pre(s);
    nch=size(eeg,2);
    vis_plot_eeg(t,eeg',pltopts);
    hold on
    plot([0 0],[-pltopts.offset*nch 200],'-k','LineWidth',2)
    plot([sz_data.duration(s) sz_data.duration(s)],[-pltopts.offset*nch 200],'-k','LineWidth',2)
    hold off
    
    % title
    set(gca,'FontSize',10)
    standard_title(sz_data,'str',opts.str)
    
%     % plot PSDs (for whole time series)
%     vis_psd(Sz_preproc, s)
%     figure(2)
%     title(['Subject ' subjectID ', Seizure ' Sz_preproc.segment_id{s}], 'Interpreter', 'none')
%     fi=[1:0.1:150];
%     pxx=pwelch(eeg,2*fs,fs,fi,fs);
%     loglog(fi,pxx)
%     title(['Subject ' subjectID ', Seizure ' Sz_preproc.segment_id{s}], 'Interpreter', 'none')

end