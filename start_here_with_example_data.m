%The overall rational of this library is that several seizure markers
%should be calculated quickly and easily. Even when new seizures are added,
%we should only need to run the analysis on these new seizures. Similarly
%when new markers are added, their calculation should re-use as much of
%what is existing already (e.g. if the new marker is based on line-length,
%you will not need to re-calculate line length again if that was done
%already for another marker). Thus the following code keeps track of
%seizures, and calculations (e.g. line length, or bandpower or similar -
%you can add to it!) in a database, and can use them for subsequent
%severity measure calculations. See the readme for more details.

addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
%% setup database and put in first dataset
% set base folder for database
basefolder = './db_severity_markers/';%will create this folder if it does not exist. If you want to start from scratch, change this path to something that doesn't exist yet!

% segment database connect
seg_db = Seg_db([basefolder '/seg_db']);     % create database 


%add subject ID02's seizures (in total 2)
load('SWEZ_data/ID02_Sz_preproc.mat')
seg_db.add(Sz_preproc);                  % add segments
metadata = seg_db.get_meta('ind',-1);    % get only meta data so far - not the seizure EEG itself
data = seg_db.get_data('ind',-1);        % getting all the data so far including the sz EEG

    
%% calculate basic measures like LL, energy, bp etc in 1sec window, no overlap

sampling_rate=metadata.segment_fs(1);%just using the firts, as all sampling was the same in our data after preproc.

linelength_db = LL_db([basefolder '/LL_db']); %setup folder for all Line Length measures
linelength_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0);
linelength_db.paramset_tbl                       % display all currently tracked paramsets
linelength_db.calc(data,[]);                       % calculate all parametersets for all segments in data


energy_db = Energy_db([basefolder '/Energy_db']);
energy_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0);
energy_db.paramset_tbl                       % display all currently tracked paramsets
energy_db.calc(data,[]);                       % calculate all parametersets for all segments in data


%bands: [1 4; 4 8; 8 13; 13 30; 30 60, 60 100]; 
bandpower_db = BP_db([basefolder '/BP_db']);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[1 4]);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[4 8]);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[8 13]);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[13 30]);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[30 60]);
bandpower_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',0,'bandbounds',[60 100]);
bandpower_db.paramset_tbl                       % display all currently tracked paramsets
bandpower_db.calc(data,[]);                       % calculate all parametersets for all segments in data


ampl_db = Ampl_db([basefolder '/ampl_db']);         % calculated the signal range, i.e. max - min of signal in a window
ampl_db.add_paramset();                      % add default paramset
ampl_db.paramset_tbl                       % display all currently tracked paramsets
ampl_db.calc(data,[]);                       



%% add new patient

load('SWEZ_data/ID03_Sz_preproc.mat')
seg_db.add(Sz_preproc);                  % add segments
metadata = seg_db.get_meta('ind',-1);    % see how meta data updated?
data = seg_db.get_data('ind',-1);        % see how data updated?



%% re-run all the measures - will only rerun on new seizures


linelength_db.calc(data,[]);%only need to re-run this line to re-run the line-length with the previous setting on the new seizures
energy_db.calc(data,[]);%same for the other measures
bandpower_db.calc(data,[]);  
ampl_db.calc(data,[]); 


%% add new setting for line length and re-calculate
%we won't need this, but just for demo
linelength_db.add_paramset('wndw_len',sampling_rate,'wndw_overlap',sampling_rate/2);%this time, we use 1sec windows with 0.5 sec overlap
linelength_db.paramset_tbl                       % display all currently tracked paramsets
linelength_db.calc(data,[]);                       % calculate all parametersets for all segments again



%% get calcs
calcs_ll = linelength_db.get(1);    % get all calculation outputs in variable. my parameset of interest was the 1st in this case, change accordingly if your parameter set was a later one
calcs_energy = energy_db.get(1);
calcs_bp_delta = bandpower_db.get(1);
calcs_bp_theta = bandpower_db.get(2);
calcs_bp_alpha = bandpower_db.get(3);
calcs_bp_beta = bandpower_db.get(4);
calcs_bp_gamma = bandpower_db.get(5);
calcs_bp_hgamma = bandpower_db.get(6);
calcs_ampl = ampl_db.get(1); 




%% calculate all the univariate seizure severity measures, "top markers"

[top_ll]=ms_95pctle(metadata,calcs_ll);
[top_energy]=ms_95pctle(metadata,calcs_energy);
[top_delta]=ms_95pctle(metadata,calcs_bp_delta);
[top_theta]=ms_95pctle(metadata,calcs_bp_theta);
[top_alpha]=ms_95pctle(metadata,calcs_bp_alpha);
[top_beta]=ms_95pctle(metadata,calcs_bp_beta);
[top_low_gamma]=ms_95pctle(metadata,calcs_bp_gamma);
[top_high_gamma]=ms_95pctle(metadata,calcs_bp_hgamma);
%Join together into table
top = table(top_ll,top_energy,top_delta,...
    top_theta,top_alpha,top_beta,top_low_gamma,...
    top_high_gamma);

TBL_top=[metadata.segment_id top];
TBL_top.Properties.VariableNames{1}='segment_id';
%% calculate imprint and recruitment markers

val_tbl=[calcs_ll.LL_ms calcs_energy.energy_ms calcs_bp_delta.bp calcs_bp_theta.bp calcs_bp_alpha.bp calcs_bp_beta.bp calcs_bp_gamma.bp calcs_bp_hgamma.bp];
[imprint_out,cell_imprint,cell_t,cell_madscores] = ms_imprint(metadata,val_tbl,calcs_ll.t_wndw);  % using the same t_wndw for all features as using same window length or overlap


% [onset,recruitment,chanr] = ms_recruitment(metadata,cell_imprint,cell_t);

TBL_onset_recr=[metadata.segment_id imprint_out];
TBL_onset_recr.Properties.VariableNames{1}='segment_id';
%% calculate post-ictal suppression markers - using ampl_db

tbl_suppr = ms_suppr(metadata,calcs_ampl.ampl_ms,calcs_ampl.t_wndw);

%% merge all severity markers into one & write into xlsx file
%remove channel label from this spreadsheet in order to write it to
%csv/excel
metadata_wo_chlabel=metadata;
for s=1:size(metadata,1)
    metadata_wo_chlabel.segment_channel_labels{s}=size(data.segment_data{s},1);
end
metadata_wo_chlabel.Properties.VariableNames{'segment_channel_labels'}='num_total_chan';


export_tbl = tbl_join(metadata_wo_chlabel,TBL_top);  %checked before that these are in order
export_tbl = tbl_join(export_tbl,TBL_onset_recr);  
export_tbl = tbl_join(export_tbl,tbl_suppr);  


writetable(export_tbl,'SeverityTable.xlsx')

%% visualise some seizures

%visualise second seizure in database in terms of imprint
s=2;
ms_imprint(data(s,:),val_tbl(s,:),calcs_ll.t_wndw(s),'fig_ind',1);


%visualise fifth seizure in database in terms of post ictal suppression
s=5;
ttbl=ms_suppr(data(s,:),calcs_ampl.ampl_ms(s),calcs_ampl.t_wndw(s),'fig_ind',2);

%% Next-time re-use

%Once you have calculated a database you can load it again easily.
%let's clear all the variables and re-connect to the database:
clear all
close all

basefolder = './db_severity_markers/';%will create this folder if it does not exist. If you want to start from scratch, change this path to something that doesn't exist yet!
seg_db = Seg_db([basefolder '/seg_db']);     % create database 
metadata = seg_db.get_meta('ind',-1);    % get only meta data so far - not the seizure EEG itself
data = seg_db.get_data('ind',-1);        % getting all the data so far including the sz EEG

%to connect to the LL database, simply do:
linelength_db = LL_db([basefolder '/LL_db']); %setup folder for all Line Length measures
%and now you can get the calculation outputs in variable again:
calcs_ll = linelength_db.get(1);  
%now you can do with it what you want, like create new severity measures!
