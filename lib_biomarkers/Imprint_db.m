classdef Imprint_db < Abstract_calc_db
    
methods (Static = true)

function [default_params,vers_id,vers_cmmt,name_calc] = get_algo_params()
    % FILL IN these with your values
    
    % name of calculation performed by database
    name_calc = 'imprint';     % best keep it short
    
    % parameters 
    default_params.wl = 1;
    default_params.ol = 0;
    default_params.rec_thresh = 0.1;
    default_params.ict_buffer = 10;
    default_params.mad_thresh = 3;
    default_params.use_individual_pre = 1;

    % example version
    vers_id = 'v1';
    vers_cmmt = 'initial version';
   
end


function tbl = calc_seg(data_tbl,params,fig_hndl) 
    
    % mandatory: calculations for whole segment
    nv_pairs = namedargs2cell(params);
    [imprint,rawfeat_allsz,scoredfeat_allsz,pre_ids,ictal_ids] = sz_imprint(data_tbl, nv_pairs{:});

    % create output
    tbl = table(imprint,rawfeat_allsz,scoredfeat_allsz,pre_ids,ictal_ids, ...
        'VariableNames',{'imprint','rawfeat_allsz','scoredfeat_allsz',...
        'pre_ids','ictal_ids'});
    
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
function [tbl,time_series] = plot(data,varargin)
    % for documentation look at Abstract_calc_db.plot
    varargin = [varargin {'subclass'} {mfilename('class')}];
    [tbl,time_series] = Abstract_calc_db.plot(data,varargin{:});
end

end
end