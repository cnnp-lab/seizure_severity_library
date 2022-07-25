function  [band_power] = univarsw_band_power(data, fs, window, overlap, bandbounds)   
%       INPUTS:
%       data: matrix of signals (# channels x # time points)
%       window: window size in terms of number of samples
%       overlap: (optional) number of samples that consecutive windows overlap
%       bandbounds: (optional) bounds for freq bands
% 
%       OUTPUTS: 
%       band_power as a 3D matrix of # channels x # windows x # bands
%
% Sarah Jane Gascoigne/Yujiang Wang 23 July 2021
% Newcastle University School of Computing 
% To do: make psd calculation optional inputs?
%
    arguments % Declare the class and size (where appropriate of the input 
                % variables
        data double {mustBeNumeric,mustBeReal}
        fs (1,1) double {mustBeNumeric,mustBeGreaterThan(fs,0)} 
        window (1,1) double {mustBeNumeric,mustBeGreaterThan(window,0)} 
        overlap (1,1) double {mustBeNumeric,mustBeLessThan(overlap,window),mustBeGreaterThanOrEqual(overlap,0)} = 0 %Default overlap is 0
        bandbounds (:,2) double {mustBeNumeric,mustBeReal} = def_classic_bandbounds()
    end

    % number of time series/signals/channels
    num_channels=size(data,1);    

    % windows for computing measures
    series_length=size(data,2); % length of time series
    
    % number of freq bands requested
    num_bands=size(bandbounds,1);

    % find index of star and stop times for each window
    w_start=1:(window-overlap):(series_length-(window-1));
    w_stop=w_start+(window-1);
    
    num_windows=length(w_start); % number of windows

    % initialize array for holding band powers
    band_power=zeros(num_channels,num_windows,num_bands);

    % compute line length for each window

    for i=1:num_windows%for each window
        % window start & stop ids
        w = w_start(i):w_stop(i); 
        for b=1:num_bands % for each band

            % calculate band power 
            band_power(:,i,b) = bandpower(data(:,w)',fs, bandbounds(b,:));
            
        end

    end




end
