classdef YourClassName < Abstract_calc_db
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params()
    % FILL IN these with your values
    
    % name of calculation performed by database
    name_calc = 'calc';     % best keep it short
    
    % parameters 
    default_params.wndw_len = 256;
    default_params.wndw_shift = 64;
    default_params.a = 2;

    % example version
    vers_id = 'v1';
    vers_cmmt = 'initial version';
    
    error('Not implemented.')
end


function tbl = calc_seg(data_tbl,params,fig_hndl) 
    
    % mandatory: calculations for whole segment
    error('not implemented')

    % optional: plotting on segment level
    if ~isempty(fig_hndl)
        figure(fig_hndl)
        warning('No plotting on segment level available');
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
    [tbl,time_series] = Abstract_calc_db.plot(data,varargin{:});
end

end
end