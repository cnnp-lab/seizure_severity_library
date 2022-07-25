function [tbl_plhg] = ms_PLHG(data_tbl,opts)
% calculates phase-locked high gamma. Implementation based on Diamond 2021
% code (accessed Feb 2022).
% 
% inputs:
%   - data_tbl: full data table - consider splitting into calcs and measure
%   in future, e.g. we calculate HG already.
%   - t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below
%outputs are in order of meta_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl

%         opts.rec_thresh (1,1) double {mustBeNumeric, mustBePositive} = 0.1
%         opts.ict_buffer (1,1) double {mustBeInteger, mustBePositive} = 10 %in seconds
%         opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
%         opts.hsc       (1,1) double {mustBeNumeric, mustBePositive} = 0.1 %percentage of seizure to be ignored at the start for calculation of max spread
        opts.fig_ind = -1 % -1: no plotting, 0:new figure, k>0: figure with index k
        opts.CVthresh (1,1) double {mustBeNumeric, mustBePositive} = 2.5
    end
    
    %fill in optional arguments
%     rec_thresh=opts.rec_thresh;    
%     mad_thresh=opts.mad_thresh;
%     ict_buffer=opts.ict_buffer;
%     hsc=opts.hsc;
    fig_ind=opts.fig_ind;
    CVthresh=opts.CVthresh;
    
    %% fill in some info from meta_tbl


    nsegs = size(data_tbl,1);
    
    tbl_plhg = cell(nsegs,2);
%     nchr = zeros(nsegs,1);
%     nmaxchr = zeros(nsegs,1);
%     tmaxchr = zeros(nsegs,1);
%     cell_madscores = cell(nsegs,1);
%     cell_t = cell(nsegs,1);
    
    for s = 1:nsegs
        
        sztbl=data_tbl(s,:);
        Fs = sztbl.segment_fs;            % Sampling frequency
        szd=sztbl.duration;
        
        if szd>6

            FO=round(Fs/2); %filter order
            uf=min(150,floor(Fs/2.1));
            tsdata=sztbl.segment_data{1};
            tsBaseline = tsdata(:,1:(sztbl.segment_pre*Fs))';
            ts = tsdata(:,(sztbl.segment_pre*Fs)+[1:szd*Fs])';
            %% Preparing our values  

            sizeWindow = 3*Fs; 
            sizeSkip = round(Fs/2); 


            stepsBuffer = buffer(1:size(ts,1),sizeWindow,sizeWindow - sizeSkip,'nodelay');
            timeVals = stepsBuffer(1,:); 


            myFilter = buildFilter(Fs, [4 30],FO);
            filter_4_30 = filtfilt(myFilter,1,ts);
            phi_4_30 = angle(hilbert(filter_4_30));


            myFilter = buildFilter(Fs, [80 uf],FO);
            filter_80_150 = filtfilt(myFilter,1,ts);
            a_80_150 = abs(hilbert(filter_80_150));
            phi_a_80_150 = angle(hilbert(a_80_150));



            %% Pre seizure baseline 

            filter_80_150_baseline = filtfilt(myFilter,1,tsBaseline);
            a_80_150_baseline = abs(hilbert(filter_80_150_baseline));
            a_80_150_baseline = mean(a_80_150_baseline);

            %% Computation of PLV 
            PLVMaster = zeros(size(stepsBuffer,2),size(ts,2)); 
            plhgMaster = zeros(size(PLVMaster));



            for jj = 1:size(stepsBuffer,2)
                currentTime = stepsBuffer(:,jj);
                currentTime = currentTime(logical(currentTime));    

                phi_4_30_jj = phi_4_30(currentTime,:);
                phi_a_80_150_jj = phi_a_80_150(currentTime,:);

                PLV = abs(mean(exp(1i * (phi_4_30_jj - phi_a_80_150_jj))));

                a_80_150_jj = a_80_150(currentTime,:); 
                a_80_150_norm_jj = mean(a_80_150_jj) ./ a_80_150_baseline;

                PLHG = a_80_150_norm_jj .* PLV;

                PLVMaster(jj,:) = PLV; 
                plhgMaster(jj,:) = PLHG;

            end

            plhgMaster = zscore(plhgMaster);

            %% Significant leads 

            maxVal = max(plhgMaster); 

            [~,mu, sigma] = zscore(maxVal); 
            coeffVar = sigma / mu; 
            sigThresh = mu + coeffVar * CVthresh; 

            sigLeads = maxVal > sigThresh; 

            %% Key times 

            sigTime = nan(1,size(plhgMaster,2));
            for jj = 1:length(sigTime)
                currentInd = find(plhgMaster(:,jj) >= sigThresh,1); 
                if isempty(currentInd); continue; end

                sigTime(jj) = timeVals(currentInd) / Fs; 
            end


            %% write output
            tbl_plhg{s,1}=sigLeads;
            tbl_plhg{s,2}=sigTime;
        end
        
    end
    
    tbl_plhg=cell2table(tbl_plhg);
    tbl_plhg.Properties.VariableNames{1} = 'has_plhg';
    tbl_plhg.Properties.VariableNames{2} = 'when_plhg';
    
    
    %doing plotting here if only one input is given
    if nsegs == 1 && fig_ind > -1
        if fig_ind == 0
            figure()
        else
            figure(fig_ind)
        end
        ax1 = subplot(12,1,1:6);
        plot_sz(data_tbl(s,:),'fig_ind',0.5)%you'd have to pass in the data table here as meta_tbl!
        standard_title(data_tbl(s,:),'str', 'PLHG detection')

        ax2 = subplot(12,1,7:8);
        imagesc(PLVMaster','XData',timeVals./Fs)
        ylabel('PLV Master')
        colorbar
        
        ax3 = subplot(12,1,9:10);
        imagesc(plhgMaster','XData',timeVals./Fs)
        ylabel('PLHG Master')
        colorbar

        ax4 = subplot(12,1,11:12);
        imagesc((plhgMaster>sigThresh)','XData',timeVals./Fs)
        ylabel('PLHG thresholded')


        linkaxes([ax1,ax2, ax3, ax4],'x')
        
        
    end
    
    

    
end