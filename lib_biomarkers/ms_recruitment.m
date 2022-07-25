function [tbl_onset,tbl_recruitment,optout_chanrecruit] = ms_recruitment(meta_tbl,tbl_imprint,tbl_t,opts)
% calculates "imprint" and other (initial) recruitment metrics of seizures.
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

        opts.avgtw (1,1) double {mustBeNumeric, mustBePositive} = 5%5 seconds needs to average over to detect persistent recruitment.
        opts.rthresh (1,1) double {mustBeNumeric, mustBePositive} = 0.6 %average number of channels needed to be recruited in the time set in avgtw. E.g. 1 channel needs to be recruited in 5 sec = 1/5 = 0.2
        opts.fig_ind = -1 % -1: no plotting, 0:new figure, k>0: figure with index k
    end
    
    %fill in optional arguments
    avgtw=opts.avgtw;    
    rthresh=opts.rthresh;
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
        
        disp(['Calculating recruitment in ' meta_tbl.segment_id{s}])
        
        %read out seizure segment & compress all features into one 3D array
        imprint=tbl_imprint{s};
        tw=tbl_t{s};%timing for each seizure, in actual seconds
        if size(imprint,2)~=size(tw,2)
            error('mismatched dimensions in imprint and time')
        end        
        if numel(imprint)>0 && length(tw)>=avgtw
           
        avgw = avgtw/(tw(2)-tw(1));
        
        
        
        S=sum(imprint,1);
        DS=diff(S);
        MDS=movmean(DS,[avgw-1 0]);
        recruiting=[0 MDS>rthresh 0];
        DR=diff(recruiting);
        
        
        
        
        
        if ~isempty(find(DR==1,1))%if any recruitment happens at all
            fid_start=find(DR==1,1);
            fid_end=find(DR==-1,1);
            fid_end=min(fid_end,length(S));
            recruitphase=zeros(size(S));
            recruitphase(fid_start:fid_end)=1;
            recruitphase=logical(recruitphase);
    
            MS=S(recruitphase);
            [m,mid]=max(MS);
            tbl_recruitment(s,1)=m;%max number of channels recruited in this phase of first detected recruitment.
            tbl_recruitment(s,2)=numel(MS)*(tw(2)-tw(1));%number seconds in this phase of first detected recruitment.
            tbl_recruitment(s,3)=tw(fid_start);%time of recruitment start
            tbl_recruitment(s,4)=tw(fid_start+mid-1);%time of max recruitment
            
            optout_chanrecruit{s}=imprint(:,fid_start+mid-1);
            
            fid=find(MS>0,1);
            tbl_recruitment(s,5)=MS(fid);%number of channels at first instance of recruitment
        else
            tbl_recruitment(s,1)=NaN;
            tbl_recruitment(s,2)=NaN;
            tbl_recruitment(s,3)=NaN;
            tbl_recruitment(s,4)=NaN;
            tbl_recruitment(s,5)=NaN;
            fid_start=size(imprint,2);%this is just so we can still calculate an onset number of channels
        end
        
        
        
        
        %now get the onset, which is anything before the first recruitment
            %phase
            MO=imprint(:,1:fid_start-1);
            MOc=sum(MO,2)>=1;
            tbl_onset(s)=sum(MOc); %number of onset channels
        else
             warning('No imprint found for this seizure')
             tbl_recruitment(s,1)=NaN;
            tbl_recruitment(s,2)=NaN;
            tbl_recruitment(s,3)=NaN;
            tbl_recruitment(s,4)=NaN;
            tbl_recruitment(s,5)=NaN;
        end
    end
    
    
    tbl_onset=table(tbl_onset,'VariableNames',{'num_onset_chan'});
    tbl_recruitment=table(tbl_recruitment(:,1),tbl_recruitment(:,2),tbl_recruitment(:,3),tbl_recruitment(:,4),tbl_recruitment(:,5),'VariableNames',{'num_max_recruit_chan','dur_recruit','time_recruit_start','time_max_recruit','num_chan_start_recruit'});
    
    %doing plotting here if only one input is given
    if nsegs == 1 && fig_ind > -1
        if fig_ind == 0
            figure()
        else
            figure(fig_ind)
        end
        ax1 = subplot(12,1,1:5);
        plot_sz(meta_tbl(s,:),'fig_ind',0.5)%you'd have to pass in the data table here as meta_tbl!
        standard_title(meta_tbl(s,:),'str', 'Onset/Recruitment detection')

        ax2 = subplot(12,1,6:8);
        imagesc(imprint,'XData',tw)
        title('imprint in each channel')

        ax3 = subplot(12,1,9);
        plot(tw,S)
        title('S - summed imprint in all channels')

        ax4 = subplot(12,1,10);
        plot(tw(1:end-1)+(tw(2)-tw(1)),DS)
        title('dS - derivative of S')
        

        ax5 = subplot(12,1,11);
        plot(tw(1:end-1)+(tw(2)-tw(1)),MDS)
        title('mdS - averaged dS')

        ax6 = subplot(12,1,12);
        plot(tw,recruiting(1:end-1))
        title('mdS over a certain threshold')

        linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')
        
        
    end
    
    
end