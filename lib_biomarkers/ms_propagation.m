function [tbl_onset,tbl_recruitment,optout_chanrecruit] = ms_propagation(meta_tbl,tbl_imprint,tbl_t,opts)
% calculates propagation metrics based on imprint of seizures.
% 
% inputs:
%   - meta_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - tbl_imprint: output from ms_imprint
%   - tbl_t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below
%outputs are in order of meta_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        meta_tbl
        tbl_imprint
        tbl_t

        %opts.avgtw (1,1) double {mustBeNumeric, mustBePositive} = 5%5 seconds needs to average over to detect persistent recruitment.
        %opts.rthresh (1,1) double {mustBeNumeric, mustBePositive} = 0.6 %average number of channels needed to be recruited in the time set in avgtw. E.g. 1 channel needs to be recruited in 5 sec = 1/5 = 0.2
        opts.fig_ind = -1 % -1: no plotting, 0:new figure, k>0: figure with index k
    end
    
    %fill in optional arguments
    %avgtw=opts.avgtw;    
    %rthresh=opts.rthresh;
    fig_ind=opts.fig_ind;
    
    %% fill in some info from meta_tbl


    if size(meta_tbl,1) ~= size(tbl_imprint,1) || size(meta_tbl,1) ~= numel(tbl_t) || ...
        size(tbl_imprint,1) ~= numel(tbl_t)
        error('All inputs must have same amount of rows.')
    end
    
    nsegs = size(meta_tbl,1);
    
    
    
    %outputs
    tbl_onset = zeros(nsegs,1);
    tbl_recruitment = zeros(nsegs,5);
    optout_chanrecruit = cell(nsegs,1);
    
    
    for s = 1:nsegs
        
        disp(['Calculating propagation in ' meta_tbl.segment_id{s}])
        imprint=tbl_imprint{s};
        tw=tbl_t{s};%timing for each seizure, in actual seconds
        if size(imprint,2)~=size(tw,2)
            error('mismatched dimensions in imprint and time')
        end        
        if numel(imprint)>0
            
            nchinv=sum(imprint,1);%get n channels in each window
            kid=nchinv>1;%remove 0 recruitement time windows
            min(nchinv(kid))
            max(nchinv(kid))
        end
        
        
        
    end
    
    
    
    
    
    
    
    
    
    
end