function  [energy] = univarsw_energy(data, window, overlap, demean)   
%       INPUTS:
%       data: matrix of signals (# channels x # time points)
%       window: window size in terms of number of samples
%       overlap: number of samples that consecutive windows overlap
% 
%       OUTPUTS: 
%       energy: line length of each window of each signal (# channels x
%           # windows)
%
% Gabrielle M. Schroeder 12 June 2017
% Yujiang Wang 27th July 2021
% Newcastle University School of Computing


    arguments % Declare the class and size (where appropriate of the input 
                % variables
        data double {mustBeNumeric,mustBeReal}
        window (1,1) double {mustBeInteger,mustBeGreaterThan(window,0)} 
        overlap (1,1) double {mustBeInteger,mustBeLessThan(overlap,window),mustBeGreaterThanOrEqual(overlap,0)} = 0 %Default overlap is 0
        demean (1,1) logical {mustBeNumericOrLogical} = 1 %default is to demean
    end

    % number of time series/signals/channels
    numChannels=size(data,1);    

    % windows for computing measures
    seriesLength=size(data,2); % length of time series

    % find index of start and stop times for each window
    wStart=1:(window-overlap):(seriesLength-(window-1));
    wStop=wStart+(window-1);
    
    numWindows=length(wStart); % number of windows

    % initialize array for holding line length
    energy=zeros(numChannels,numWindows);

    % compute energy for each window
    for i=1:numWindows

        % window start & stop ids
        w = wStart(i):wStop(i); 

        if demean==1
            d=detrend(data(:,w)','constant')';
        else
            d=data(:,w);
        end
    
        % average signal energy
        energy(:,i)=sum(d.^2,2)./window;

    end


end
