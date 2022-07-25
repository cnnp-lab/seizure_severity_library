classdef Energy_db < Abstract_calc_db
    % calculates some measure for connected to energy in sliding
    % windows. Current version takes the energy in
    % sliding windows after detrending: energy=sum(detrended_data.^2)./window;
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params()
    % FILL IN these with your values
    
    % name of calculation performed by database
    name_calc = 'energy';     % best keep it short
    
    % parameters 
    default_params.wndw_len = 128;
    default_params.wndw_overlap = 0;

    % example version
    vers_id = 'v1';
    vers_cmmt = 'energy=sum(detrended_data.^2,2)./window';
    
end


function tbl = calc_seg(data_tbl,params,fig_hndl) 
    
    % assign params
    wndw_len = params.wndw_len;
    wndw_overlap = params.wndw_overlap;
    fs = data_tbl.segment_fs(1);
    pre = data_tbl.segment_pre(1);
    data = data_tbl.segment_data{1};
    nchns = size(data,1);
    
    % mandatory: calculations for whole segment
    [block_ind, ~, t_wndw] = sliding_window(size(data,2), wndw_len, 'wndw_overlap',wndw_overlap,...
    'fs',fs, 'pre',pre);
    nwnds = numel(t_wndw);
    data_wndw = zeros(nchns, nwnds);
    % data = abs(data);
    for w = 1:nwnds
        ddata=detrend(data(:,block_ind(w,1):block_ind(w,2)),'constant')';
        data_wndw(:,w)=sum(ddata.^2,1)./wndw_len;
    end
    
    % create output table
    tbl = table({data_wndw},{t_wndw},'VariableNames',{'energy_ms','t_wndw'});

    % optional: plotting on segment level
    if ~isempty(fig_hndl)
        figure(fig_hndl);
        
        ax1 = subplot(211);
        plot_sz(data_tbl,'str','energy plot','fig_ind',0.5)
        ax2 = subplot(212);
        imagesc(data_wndw,'XData',t_wndw)
        title('Heatmap of energy')
        linkaxes([ax1 ax2],'x');
        colorbar
    end
end

function tbl = calc_chn(data_tbl,params,chn,fig_hndl)
      
    % optional: calculations for specific channel with index content
    warning('No calculations on channel level available')
    tbl = [];

    % optional: plotting on channel level
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        warning('No plotting on channel level available'); 
    end
end


% Don't change anything here
function [tbl,time_series] = plot(data,varargin)
    % for documentation look at Abstract_calc_db.plot
    varargin = [varargin {'subclass'} {mfilename('class')}];
    [tbl,time_series] = Abstract_calc_db.plot(data,varargin{:});
end

end
end