function [tbl_imprint_out,cell_imprint,cell_t,cell_madscores] = ms_imprint(meta_tbl,val_tbl,t,opts)
% calculates "imprint" and other (initial) recruitment metrics of seizures.
% 
% inputs:
%   - meta_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - val_cell: cell array with data matrix chns-by-time for each segment,
%       cell columns are different features, rows are different segments
%   - t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below
%outputs are in order of meta_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        meta_tbl
        val_tbl
        t

        opts.rec_thresh (1,1) double {mustBeNumeric, mustBePositive} = 0.1
        opts.ict_buffer (1,1) double {mustBeInteger, mustBePositive} = 10 %in seconds
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
        opts.hsc       (1,1) double {mustBeNumeric, mustBePositive} = 0.1 %percentage of seizure to be ignored at the start for calculation of max spread
        opts.fig_ind = -1 % -1: no plotting, 0:new figure, k>0: figure with index k
    end
    
    %fill in optional arguments
    rec_thresh=opts.rec_thresh;    
    mad_thresh=opts.mad_thresh;
    ict_buffer=opts.ict_buffer;
    hsc=opts.hsc;
    fig_ind=opts.fig_ind;
    
    mc=-1/(sqrt(2)*erfcinv(3/2)); % fixed factor for MAD score calculation
    
    %% fill in some info from meta_tbl

%     nch = length(meta_tbl.segment_channel_labels{1});%!!!!!!!
%     fs=meta_tbl.segment_fs(1);
%     nsz = size(meta_tbl,1);%number of seizures
    
   
    if size(meta_tbl,1) ~= size(val_tbl,1) || size(meta_tbl,1) ~= numel(t) || ...
        size(val_tbl,1) ~= numel(t)
        error('All inputs must have same amount of rows.')
    end
    
    nsegs = size(meta_tbl,1);
    
    cell_imprint = cell(nsegs,1);
    nchr = zeros(nsegs,1);
    nmaxchr = zeros(nsegs,1);
    tmaxchr = zeros(nsegs,1);
    cell_madscores = cell(nsegs,1);
    cell_t = cell(nsegs,1);
    
    for s = 1:nsegs
        
        
        
        %read out seizure segment & compress all features into one 3D array
        features=log(cat(3,val_tbl{s,:}));%assuming all features require log transform - for now ok, needs changing maybe at input level
        fs=meta_tbl.segment_fs(s);%sampling freq
        tw=t{s};%timing for each seizure, in actual seconds
        
        
        %pull out segment ids
        pre_ids=find(tw<=-1*ict_buffer);
        ictal_ids=find(  tw>0 & tw<=meta_tbl.duration(s)  );
        
        %------------------------------------------------------------------
        %Calculate MAD score for each feature & ignores preictal spikes by
        %ignoring outliers
        disp(['Calculating feature MAD scores in ' meta_tbl.segment_id{s}])
        pre_features=features(:,pre_ids,:);
%         pre_outl=isoutlier(pre_features,2,'ThresholdFactor',5);
%         pre_features(pre_outl)=NaN;
        feat_pre_m=median(pre_features,2,'omitnan');
        feat_pre_smad=mc*mad(pre_features,1,2);
        mad_score_features = (features-feat_pre_m)./feat_pre_smad;%score ictal to median & scaled mad
        
        % Pull out ictal segment
        ict_features= features(:,ictal_ids,:);
        ict_features_mad = mad_score_features(:,ictal_ids,:);
        
        
        %------------------------------------------------------------------
        %calculate binary imprint measure based on all features
        
        seizure_label = meta_tbl.segment_id{s};
        duration = meta_tbl.duration(s);
        wl=(tw(2)-tw(1));
        recruitment_threshold = ceil(rec_thresh*duration/wl);
        
        % Cap recruitment threshold to be between 2 and 5 seconds
        if recruitment_threshold < 2/wl 
            recruitment_threshold = 2;
        elseif recruitment_threshold > 5/wl
            recruitment_threshold = 5;
        end
        fprintf('Seizure %s recruitment threshold: %d seconds \n',seizure_label,recruitment_threshold)

                
        %if at least 2 features are abnormal in madscore, movum them & get
        %recruited channel at each point.
        SA=max(ict_features_mad,[],3)>mad_thresh;%max abnormality across the features, then check if that max is above threshold
        ms_a=movsum(SA,[0 recruitment_threshold-1],2);%forward looking sum
        ms_b=movsum(SA,[recruitment_threshold-1 0],2);%backward looking sum
        %plot(sum(ms==recruitment_threshold))
%         imprint=(ms==recruitment_threshold);
        imprint=((ms_a+ms_b)>=recruitment_threshold);%combining the backward and forward looking sum has the advantage of not cutting off early activity, or neglecting late activity
        
        
            nchr(s)=sum(sum(imprint,2)>=1);%get number of channels ever in imprint
            
            nt=size(imprint,2);%number of time windows
            sid=max(round(nt*hsc),1);%how many windows to ignore at the start of the seizure, 1 or hsc*nt, whichever is larger

            nchinv=sum(imprint,1);%number of channels in imprint at each time window
        if numel(nchinv(sid:end))>1
            [m,fid]=max(nchinv(sid:end));%find max number of channels in 
            nmaxchr(s)=m;%store the max number of channels
            tmaxchr(s)=tw(ictal_ids(fid+sid-1));%store the timing of when max occurred (accounting for windows we ignored at the start)
        end
        
        % write output
        cell_imprint{s}=imprint;
        cell_madscores{s}=ict_features_mad;
        cell_t{s}=tw(ictal_ids);
    end
    
    tbl_imprint_out=table(nchr, nmaxchr, tmaxchr,'VariableNames',{'num_chan_imprint','num_chan_max_concurrent','time_till_max_concurrent'});
    
    
    %doing plotting here if only one input is given
    if nsegs == 1 && fig_ind > -1
        if fig_ind == 0
            figure()
        else
            figure(fig_ind)
        end
        ax1 = subplot(12,1,1:5);
        plot_sz(meta_tbl(s,:),'fig_ind',0.5)%you'd have to pass in the data table here as meta_tbl!
        standard_title(meta_tbl(s,:),'str', 'Imprint detection')

        ax2 = subplot(12,1,6:7);
        imagesc(SA,'XData',tw(ictal_ids))
        ylabel('max abnormality in any feature')

        ax3 = subplot(12,1,8:9);
        imagesc(ms_a+ms_b,'XData',tw(ictal_ids))
        ylabel('moving average of above')

        ax4 = subplot(12,1,10:12);
        imagesc(imprint,'XData',tw(ictal_ids))
        %plot time of max
        hold on
        plot([tmaxchr(s) tmaxchr(s)],[0 size(imprint,1)],'-k')
        hold off
        ylabel('Imprint final output')


        linkaxes([ax1,ax2,ax3,ax4],'x')
        
        
    end
    
    

    
end