function [block_ind, wndw_ind, t_val] = sliding_window(data_len, wndw_len, opts)

% PURPOSE: calculate sliding window indice to facilitate sliding window
% approaches
% INPUT:    data_len = length of data
%           wndw_len = length of the window (in data points)
%           wndw_shift OR wndw_overlap = indicates window shift or overlap,
%               if both specified error will be thrown (default if none
%               inputted: overlap = 0)
%           associate = double from 0 to 1, denotes to which index of the
%               window the calculated value of the window will be assigned
%               where 0 corresponds to left (i.e. beginning of window), 1 to
%               the right and 0.5 to the middle (floored).
% OUTPUT:   block_ind = noBlocks by 2 array with start and end indices of blocks
%           t_val = indices that the blocks are associated to (on sample
%           scale!

arguments
    data_len (1,1) {mustBeInteger}
    wndw_len (1,1) {mustBeInteger}
    
    opts.wndw_shift (1,1) {mustBeInteger} = -1 
    opts.wndw_overlap (1,1) {mustBeInteger} = -1 
    opts.associate double = 0
    opts.fs (1,1) {mustBeInteger} = -1
    opts.pre (1,1) {mustBeInteger} = -1
end

% extract wndw_shift
if opts.wndw_overlap == -1 && opts.wndw_shift == -1
    wndw_shift = wndw_len;
elseif opts.wndw_overlap ~= -1 && opts.wndw_shift ~= -1
    error('Both overlap and shift specify. Specify either or!')
elseif opts.wndw_overlap ~= -1
    wndw_shift = wndw_len - opts.wndw_overlap;
else
    wndw_shift = opts.wndw_shift;
end
if wndw_shift > wndw_len
    warning('Some data will be unused as the window shift is larger than the window length!')
end
if wndw_shift == 0
    error('Window shift needs to be larger than 0!')
end

% calc sliding windows
if data_len < wndw_len
    error('Data length smaller than window length!')
end
nslid_wndw = ceil((data_len - wndw_len + 1) / wndw_shift);
w = (1:nslid_wndw)';
block_ind = [(w-1)*wndw_shift + 1, (w-1)*wndw_shift + wndw_len];

% calc associated points
wndw_ind = (block_ind(:,1) + floor(opts.associate * (block_ind(:,2) - block_ind(:,1))));
wndw_ind = wndw_ind';

% calc time points
if opts.fs >= 0
    t_val = (wndw_ind-1) / opts.fs;
    if opts.pre >= 0
        t_val = t_val - opts.pre;
    end
end

end
