classdef LL_db < Abstract_calc_db
% database that calculates and store the line length of segments.
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params()
    
    % name of calculation performed by database
    name_calc = 'LL';    
    
    % example parameters
    default_params.wndw_len = 256;
    default_params.wndw_overlap = 192;

    % example version
    vers_id = 'v2';
    vers_cmmt = 'changed the sliding windows';
    
end


function [tbl,time_series] = calc_seg(tbl_data,params,fig_hndl) 

    % assign params
    wndw_len = params.wndw_len;
    wndw_overlap = params.wndw_overlap;
    
    % calc
    data = tbl_data.segment_data{1};
    [LL] = univarsw_line_length(data, wndw_len, wndw_overlap);
    [~,~,t_wndw] = sliding_window(size(data,2),wndw_len,'wndw_overlap',wndw_overlap,...
        'pre',tbl_data.segment_pre(1),'fs',tbl_data.segment_fs(1));

    % assign outputs
    tbl = table({LL},{t_wndw},'VariableNames',["LL_ms","t_wndw"]);
    time_series = [];
    
    % plotting
    % Your plotting here (instead of the warnings if available)
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        % plot segment analysis
        ax1 = subplot(211);
        plot_sz(tbl_data,'str','line length analysis plot','fig_ind',0.5)
        ax2 = subplot(212);
        imagesc(LL,'XData',t_wndw)
        title('Heatmap of line length')
        linkaxes([ax1 ax2],'x');
        colorbar

        % way to perfectly align plots with colorbars, maybe try that
        % out at some point:
        % https://uk.mathworks.com/matlabcentral/answers/763711-subplot-aligment-automatic-method
    end
    
end


function tbl = calc_chn(data_tbl,params,chn,fig_hndl)
      
    % optional: calculations for specific channel with index content
    warning('No calculations on channel level available')
    tbl = [];

    % optional: plotting on channel level
    if ~isempty(fig_hndl)
        warning('No plotting on channel level available'); 
    end
end


function [tbl,parset] = plot(data,varargin)
    % for documentation look at Abstract_calc_db.plot
    subclass_str = mfilename('class');
    varargin = [varargin {'subclass'} {subclass_str}];
    [tbl,parset] = Abstract_calc_db.plot(data,varargin{:});
end

end
end