function [Sz_out, rm_sz_ind] = rm_non_part_chn(Sz_preproc,rm_sz_flag)
% removes non participating channels from Sz_preproc and returns Sz_out
% with size smaller eqausl than Sz_preproc. 
%   rm_sz_ind = indices of seizures that were removed because there were
%   no participating channenls left.

    if nargin < 2
        rm_sz_flag = true;
    end

    [~, pchns_ind] = detect_participating_chan(Sz_preproc);
    nsz = size(Sz_preproc,1);
    Sz_empty = nan(nsz,1);
    rm_sz_ind = nan(size(Sz_preproc,1),1);
    for s = 1:size(Sz_preproc,1)
        Sz_preproc = rm_chns(Sz_preproc,s,pchns_ind{s});
        if rm_sz_flag && isempty(Sz_preproc.segment_data{s})
            rm_sz_ind(s) = s;        
        end
    end
    rm_sz_ind = rmmissing(rm_sz_ind);
    Sz_preproc(rm_sz_ind,:) = [];
    Sz_out = Sz_preproc;
end