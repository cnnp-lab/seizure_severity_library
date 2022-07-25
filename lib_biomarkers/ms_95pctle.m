function [tbl_top] = ms_95pctle(meta_tbl,calcs,opts)
% calculates the max value across channels of 95th percentile
% 
% inputs:
%   - meta_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - cals: from the calc databases
%   - opts: for more info on parameters look below

    arguments
        meta_tbl
        calcs
        
        opts.pct (1,1) double {mustBeNumeric, mustBePositive} = 0.1%top 10% default
        opts.fig_ind = -1 % -1: no plotting, 0:new figure, k>0: figure with index k
    end
    
    %fill in optional arguments
    pct=opts.pct;    
    fig_ind=opts.fig_ind;
    %% fill in some info from meta_tbl

    
    nsegs = size(meta_tbl,1);
    
    tbl_top = zeros(nsegs,1);

    
    for s = 1:nsegs

        %read out seizure segment & compress all features into one 3D array
        feature=cell2mat(calcs{s,2});
        tw=cell2mat(calcs{s,3});%timing for each seizure, in actual seconds
        
        
        %pull out segment ids
        pre_ids=find(tw<=0);
        ictal_ids=find(  tw>0 & tw<=meta_tbl.duration(s)  );
        
        %------------------------------------------------------------------
        % Calculate 95th percentile in each channel icatlly
        feature_ictal=feature(:,ictal_ids);
        nepoch=size(feature_ictal,2); %number of epochs of seizure activity
        k=ceil(nepoch*pct);
        B=maxk(feature_ictal,k,2); %get the top 10% values in B for each channel
        tbl_top(s)=max(median(B,2)); %take median of top 10% (this gives 95th percentile, then
        % take max across channels
    end
    
    
    
    %doing plotting here if only one input is given
    if nsegs == 1 && fig_ind > -1
        if fig_ind == 0
            figure()
        else
            figure(fig_ind)
        end
        ax1 = subplot(12,1,1:4);
        plot_sz(meta_tbl(s,:),'fig_ind',0.5)%you'd have to pass in the data table here as meta_tbl!
        standard_title(meta_tbl(s,:),'str', 'Median of top 10% values of feature')

        ax2 = subplot(12,1,5:8);
        imagesc(feature_ictal,'XData',tw(ictal_ids))
        ylabel('feature')

        ax3 = subplot(12,1,9:12);
        barh(flipud(median(B,2)))
        ylabel('Top 10% in each channel')

        linkaxes([ax1,ax2],'x')
        
        
    end

    
end
