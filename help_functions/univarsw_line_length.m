function  [lineLength] = univarsw_line_length(data, window, overlap)   
%       INPUTS:
%       data: matrix of signals (# channels x # time points)
%       window: window size in terms of number of samples
%       overlap: number of samples that consecutive windows overlap
% 
%       OUTPUTS: 
%       lineLength: line length of each window of each signal (# channels x
%           # windows)
%
% Gabrielle M. Schroeder 12 June 2017
% Sarah Jane Gascoigne 16 July 2021
% Newcastle University School of Computing


    arguments % Declare the class and size (where appropriate of the input 
                % variables
        data double {mustBeNumeric,mustBeReal}
        window (1,1) double {mustBeInteger,mustBeGreaterThan(window,0)} 
        overlap (1,1) double {mustBeInteger,mustBeLessThan(overlap,window),mustBeGreaterThanOrEqual(overlap,0)} = 0 %Default overlap is 0
    end

    % number of time series/signals/channels
    numChannels=size(data,1);    

    % windows for computing measures
    seriesLength=size(data,2); % length of time series

    % find index of star and stop times for each window
    wStart=1:(window-overlap):(seriesLength-(window-1));
    wStop=wStart+(window-1);
    
    numWindows=length(wStart); % number of windows

    % initialize array for holding line length
    lineLength=zeros(numChannels,numWindows);

    % compute line length for each window
    for i=1:numWindows

        % window
        w = wStart(i):wStop(i); 

        % distance between consecutive points
        d = abs(diff(data(:,w),1,2));

        % average distance in window
        lineLength(:,i) = sum(d,2)./size(d,2);

    end


end
