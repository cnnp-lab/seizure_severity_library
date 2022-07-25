classdef BP_db < Abstract_calc_db
% data structure to calculate the bandpower of a ieeg signal
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params(self)
    % FILL IN
    name_calc = 'bp';
    
    default_params.wndw_len = 64;
    default_params.wndw_overlap = 0;
    default_params.bandbounds = [0.5,120];
    % default_params.bandbounds = def_classic_bandbounds;
        % also possible to put multiple bounds into one paramset if desired

    vers_id = 'v1';
    vers_cmmt = 'original';
end


function tbl = calc_seg(data_tbl,params,fig_hndl) 
       
    % calc bandpower
    [band_power] = univarsw_band_power(data_tbl.segment_data{1}, data_tbl.segment_fs(1), ...
        params.wndw_len, params.wndw_overlap, params.bandbounds) ;
    [~,~, t_wndw] = sliding_window(size(data_tbl.segment_data{1},2),params.wndw_len,...
        'wndw_overlap',params.wndw_overlap,'fs',data_tbl.segment_fs(1),'pre',data_tbl.segment_pre(1));
    
    % create output table
    tbl = table({band_power},{t_wndw},'VariableNames',{'bp','t_wndw'});

    % plotting
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        nbands = size(params.bandbounds,1);
        
        % plot the time series
        ax1 = subplot(nbands+1,1,1);
        plot_sz(data_tbl,'str','Bandpower analysis','fig_ind',0.5);
        
        % plot all bands
        ax = nan(nbands,1);
        for b = 1:nbands
            ax(b) = subplot(nbands+1,1,1+b);
            imagesc(band_power(:,:,b),'XData',t_wndw);
            title([num2str(params.bandbounds(b,1)) ' Hz to ' ...
                num2str(params.bandbounds(b,2)) ' Hz'])
            colorbar
        end
        
        % linkaxes
        linkaxes([ax(:); ax1],'x');
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
function tbl = plot(data,varargin)
    % for documentation look at Abstract_calc_db.plot
    varargin = [varargin {'subclass'} {mfilename('class')}];
    tbl = Abstract_calc_db.plot(data,varargin{:});
end

end
end