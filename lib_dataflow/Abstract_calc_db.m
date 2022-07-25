classdef (Abstract) Abstract_calc_db < handle
    properties
        name_calc           % short string that indicates the calc performed
        files_path           % path to directory at which info file and calculations
                            %   will be saved
        paramset_tbl        % table with parameters in columns. Each row 
                            %   corresponds to a set of parameters
        curr_ind_parset     % ind of latest paramset added to db
        seg_ids             % string array with all segment_ids in the db
        calc_status_mat     % array with dimensions n_seg by n_paramsets. 
                            %   Contains 0 if not calculated, 1 if
                            %   calculated and 2 if in RAM.
        RAM_cell            % cell array with dimensions of calc_status_mat.
                            %   Is either empty or contains calc specified
                            %   by column and row index.
        vers_lookup         % cell array with first column vers_id and second
                            %   column vers_comment of all tracked
                            %   versions.
    end
    
    % METHODS TO BE IMPLEMENTED IN CHILD CLASSES!
    methods (Abstract = true, Static = true)
        
        [default_params,vers_id, vers_comment,name_calc] = get_algo_params()
            % returns the default parameters and the current version of 
            % the main calculation function calc_seg as id with 
            % eventually a comment for further description. Also returns a
            % short name of the calculation performed without underscores.
            % These need to be set by the user when changes
            % are made to the function in order to correctly track data in
            % the data base! In the database, only version_id used (i.e.
            % changing version_comment has no influence)
            
        tbl = calc_seg(tbl_data,params,fig_hndl)
            % This is the heart of the database. This function performs the
            % fundamental calculation of a segment the database keeps track of.
            % It receieves the data of one segment as input and should
            % produce a table with one row with the different results in
            % the columns as output.
            % input: 
            %   - tbl_data: tbl with all data from one sz
            %   - params: struct with parameters for algo specified 
            %   - fig_hndl: figure handle for figure to plot in. If empty,
            %       no plotting (in code just call figure(fig_hndl) to
            %       select the figure if non empty).
            % ouput:
            %   - tbl: table with one row that contains the different outputs of
            %       the calculation in columns (can also contain time series in
            %       cells)
            
        tbl = calc_chn(tbl_data,params,chn,fig_hndl)
            % optional function; same as calc_seg, but analyses a specific
            % channel. Has additional input chn that is the index of the
            % channel to be calculated/plotted
    end
    
    % methods implemented in abstract class
    methods  
        function obj = Abstract_calc_db(files_path)
            % basic constructor 
            
            if isempty(files_path)
                return
            end
            
            % parse files_path
            files_path(ismember(files_path,'\')) = '/';
            if ~strcmp(files_path(end),'/')
                files_path = [files_path '/'];
                warning('"/" added to filepath!')
            end
            
            % check if struct param name
            default_params = obj.get_algo_params();
            if ismember('struct',fieldnames(default_params))
                error('Parameter name "struct" illegal.')
            end
            
            % check if database already exists
            obj.files_path = files_path;
            [~,~,~,name_calc] = obj.get_algo_params();
            info_file_path = [files_path name_calc '_info.mat'];
            if isfile(info_file_path)
                % load existing database
                d = load(info_file_path);
                vars = fieldnames(d);
                for f = 1:numel(vars)
                    obj.(vars{f}) = d.(vars{f});
                end
                disp('Connected to database.')
            else
                % create folder if not existent
                if ~exist(files_path,'dir')
                    mkdir(files_path)
                    disp('Specified directory created')
                end
                
                % create new database
                obj.name_calc = name_calc;
                nparams = numel(fieldnames(default_params));
                obj.paramset_tbl = cell2table(cell(0,1+nparams+3),'VariableNames',...
                    [{'paramset_ind'},fieldnames(default_params)',{'vers_id','vers_comment','comment'}]);
                obj.calc_status_mat = zeros(0);
                obj.curr_ind_parset = 0;
                obj.seg_ids = [];
                obj.update_props();
                obj.vers_lookup = cell(0,2);
                disp('Created new database at specified location.')
            end
            obj.files_path = files_path;
            obj.RAM_cell = cell(size(obj.calc_status_mat));     
        end

        % add 
        function params = add_paramset(self,varargin)
            % Takes params as struct and adds them to the database. NO
            % calculation performed in this step. Specify either by struct
            % OR name value parameter pairs, if both specified, struct will be
            % prefered. Also look at documentation from get for more info.
            % input:
            %   - positional arg (optional): paramset as struct. Not specified params 
            %       will be set to default in the struct
            %   - [parameter_names] (nv): give parameter names as name value
            %       pairs. If both parset_struct and name value pairs
            %       given, error.
            %   - comment (nv): optional comment to calculation as chars. A
            %       string will not work!
            % outputs:
            %   - params: struct with all parameters and comments that was
            %       added to the datbase. Can be used for later reference.
            
            % parse inputs
            p = inputParser();
            p.addParameter('comment','')
            [params,p] = self.parse_parset_namevalue(p,varargin);
            comment = char(p.Results.comment);
            
            % check if paramset already exists        
            row_ind = self.get_parset_row_ind(params);
            if row_ind > 0
                warning('Parameter set already exists in database. Nothing added.')
                params = self.get_paramset(self.paramset_tbl.paramset_ind(row_ind));
            else
                
                % add paramset to paramset_tbl
                self.curr_ind_parset = self.curr_ind_parset + 1;
                params.paramset_ind = self.curr_ind_parset;
                self.add_curr_vers();
                vers_comment = self.lookup_vers_comment(params.vers_id);
                params.vers_comment = vers_comment;
                params.comment = comment;
                self.paramset_tbl = [self.paramset_tbl; ...
                        struct2table(params,'AsArray',true)];


                % add paramset to calc_status_mat
                self.calc_status_mat = [self.calc_status_mat ...
                    zeros(numel(self.seg_ids),1)];

                % add paramset to RAM_cell
                self.RAM_cell = [self.RAM_cell cell(numel(self.seg_ids),1)];

                % save new properties
                self.update_props()
                
                % user feedback
                disp(['Added parameter set with index ' num2str(params.paramset_ind)])
            end
        end
              
        function add_segs(self,segs)
            % adds the specified segments to database. NO calculations performed. 
            % input:
            %   - segs: segments baseline
            
            % parse inputs
            segs = select_segs(segs,'ind',-1);
            
            % add segs
            if ~isempty(segs)
                % determine nonexistent segments
                if isempty(self.seg_ids)
                    existent = zeros(size(segs));
                else
                    existent = ismember(segs,self.seg_ids);
                end
                segs = segs(~existent);

                % add them to the database
                self.seg_ids = [self.seg_ids; segs];
                self.calc_status_mat = [self.calc_status_mat; ...
                    zeros(numel(segs), size(self.calc_status_mat,2))];
                self.RAM_cell = [self.RAM_cell; ...
                    cell(numel(segs), size(self.RAM_cell,2))];

                % update
                self.update_props()
            end
            disp([num2str(numel(segs)) ' segment(s) added. '...
                num2str(sum(existent)) ' segment(s) already existed.']);
        end               
 
        
        % delete
        function del_paramset(self,parsets,opts)
            % deletes parameter set specified by either its unique
            % paramset_ind or a struct that contains the parameters. Not
            % existing parameters will be filled with default values.
            % input:
            %   - parsets: paramset_ind as number or struct with
            %       params or cell array of the first two. Not specified params 
            %       will be set to default when given struct or cell with
            %       structures.
            %   - user_warning: If 1, show extra user warning that requires user
            %       to confirm deletion. Default 1.
            
            arguments
                self
                parsets
                opts.user_warning = 1
                
            end
            user_warning = opts.user_warning;
            
            % determine parset_row_ind
            parset_row_ind = parse_parsets(self, parsets, 'warning');
            if numel(parset_row_ind) > 1
                error('Deletion of multiple paramsets at once not supported.')
            end
            
            
            % user warning
            if parset_row_ind > 0
                if user_warning
                    disp('Specified parameter set(s): ')
                    disp(self.paramset_tbl(parset_row_ind,:))
                    del_flag = input(['Are you sure you want to delete all calculations'...
                        ' of the specified parameter set(s)? ["yes"]: ']);
                else
                    del_flag = 'yes';
                end

                % remove paramset from database
                if strcmp(del_flag,'yes')
                    % delete data from folder structure
                    folderpath = [self.files_path self.name_calc '_calc' ...
                        num2str(self.paramset_tbl.paramset_ind(parset_row_ind))];
                    if exist(folderpath,'dir')
                        rmdir(folderpath,'s')
                    end

                    % delete from info
                    self.paramset_tbl(parset_row_ind,:) = [];
                    self.calc_status_mat(:,parset_row_ind) = [];
                    self.RAM_cell(:,parset_row_ind) = [];

                    % save new properties
                    self.update_props()

                    % feedback
                    disp('Deleted specified data.')
                else
                    disp('No data deleted.')
                end
            else
                disp('No data deleted.')
            end
                
        end

        function del_segs(self,segs,opts)
            % deletes specified segments.
            % input:
            %   - segs: segments to delete
            %   - user_warning (nv): If 1, then extra warning user needs to
            %   confirm
            
            arguments
                self
                segs
                opts.user_warning = 1
            end
            
            % parse segs
            [segs,n_req] = self.parse_segs(segs,'warning');
            n_post = numel(segs);
            
            % find indices
            delete_ind = find(ismember(self.seg_ids,segs));
            
            % user warning
            if opts.user_warning
                delete_flag = input(['Are you sure you want to delete ' num2str(numel(segs)) ' segments? ["yes"]: ']);
            else
                delete_flag = "yes";
            end
            
            % actual deletion
            if strcmp(delete_flag,"yes")
                % delete calc from folder structure
                for s = 1:numel(delete_ind)
                    for p = 1:size(self.paramset_tbl,1)
                        if self.calc_status_mat(delete_ind(s),p) > 0
                            self.calc_status_mat(delete_ind(s),p) = ...
                                1- self.del_calc_data(self.seg_ids{delete_ind(s)},...
                                    self.paramset_tbl.paramset_ind(p));
                        end
                    end
                end

                % delete deleted segments from tables
                ind_deleted = sum(self.calc_status_mat(delete_ind,:),2) == 0;
                self.calc_status_mat(ind_deleted,:) = [];
                self.RAM_cell(ind_deleted,:) = [];
                self.seg_ids(ind_deleted) = [];
                ndel = sum(ind_deleted);
                
                % save 
                if sum(ind_deleted) > 0
                    self.update_props();
                end
                
            else
                ndel = 0;
            end
            
            % feedback
            disp([num2str(n_req-n_post) ' segment(s) not in database. '...
                num2str(ndel) ' segment(s) deleted.']);
            
%             disp(['segment(s): ' num2str(ndel) '/' num2str(nreq) ' (deleted/requested)'])
            
%             disp([num2str(numel(delete_ind)) ' segment(s) deleted. '])
%                 num2str(numel(segs)-numel(delete_ind)) ' segment(s) not existent.'])
        
        end
        
        function del_calcs(self,segs,parsets,opts)
            % deletes calculations specified by any combination of the given 
            % segments and parametersets. If exactly one input not specified,
            % all segments or all parameters sets will be selected. If no
            % args, error thrown. However, the specified entries will 
            % still exist in database but without calculation. 
            % To fully delete paramsets or segments, call
            % del_seg or del_paramset.
            % inputs:
            %   - segs: cell/string array containing ids or table with 
            %       column segment_id. If -1, all segments selected.
            %   - parsets: paramset_ind as number or struct with
            %       params or cell array of the first two. Not specified params 
            %       will be set to default.
            
            arguments
                self
                segs
                parsets
                opts.user_warning = 1
            end
            
            % input checking
            if nargin == 1
                error('Specify at least one of the two arguments!')
            end
            
            % parse segs
            if segs == -1
                segs = self.seg_ids;
            else
                segs = self.parse_segs(segs,'warning');
            end
            seg_ind = find(ismember(self.seg_ids,segs));
            
            % parse paramsets
            if nargin == 2 || isempty(parsets)
                parset_row_ind = 1:size(self.paramset_tbl,1);
            else
                parset_row_ind = parse_parsets(self, parsets);
            end
            parset_ind = self.paramset_tbl.paramset_ind(parset_row_ind);
        
            % user warning
            nreq = numel(segs) * numel(parset_ind);
            if opts.user_warning
                delete_flag = input(['Are you sure you want to delete ' num2str(nreq) ' calculations? ["yes"]: ']);
            else
                delete_flag = "yes";
            end
            
            % actual deletion
            if strcmp(delete_flag,"yes")
                % delete calculations
                for s = 1:numel(segs)
                    for p = 1:numel(parset_ind)
                        if self.calc_status_mat(seg_ind(s),parset_ind(p)) > 0
                            self.calc_status_mat(seg_ind(s),parset_ind(p)) = ...
                                1-self.del_calc_data(self.seg_ids{seg_ind(s)},parset_ind(p));
                            self.RAM_cell(seg_ind(s),parset_row_ind(p)) = cell(1);
                        end
                    end
                end
                ndel = sum(sum(self.calc_status_mat(seg_ind,parset_ind)==0));

                % update
                self.update_props();
            else
                ndel = 0;
            end
                
            % user feedback
            disp(['Segment-parameterset combination(s): ' num2str(nreq) '/' ...
                num2str(ndel) '/' num2str(nreq-ndel) '(requested/deleted/not delete)']) 
        end
        
        
        % main functionality
        function [tbl_out,parset_out] = calc(self,data_inp,paramsets,opts)
%             calculates values for all parameter sets and segments
%             specified. The segments are specified by giving their data as
%             input. By default, values are not recalculated if already
%             existent in database. If parallelisation used, please pass a
%             seg_db for the segments for optimal performance
%             inputs:
%               - data_inp: table with column 'segment_id' and
%               'segment_data' or instance of Seg_db
%               - paramsets: cellarray where each cell specifies a
%                   parameterset, i.e. each cell contains either number or
%                   parameter struct. If empty, all paramsets are selected.
%                   Sets not existing in the database will be attempted to
%                   add.
%               - id/ind/pat/row (op,nv): specifies subset of segments from data_tbl 
%                   that should be calculated. For more info about inputs look 
%                   into doc of select_segs. If empty, all segments from data_inp selected.
%               - recalc (op,nv): boolean. If true, then already calculated
%                   values are recalculated.
%               - parfor (op,nv): boolean. If true, calculation parallelized.
%                   However, good choice for batchsize required for actual speed! 
%                   advantage! Requires Parallel Computing Toolbox!
%               - batchsize (op,nv): integer. Denotes number of calculations
%                   that are passed to a parpool worker in the parfor loop.
%                   To take place, still needs opts.parfor flag true.
%               - only_seg_in_db (op,nv): Applies after selecting segments via
%                   indexing. If true, algo will not calculate
%                   segments not already contained in the database. If false
%                   (default), the specified segments are added to the database
%             outputs:
%               - tbl_out: if only one parameterset specified, algo will
%                   return the result table for this paramset
            
            
            arguments
                self
                data_inp
                paramsets
                opts.recalc = false
                
                % specify subset of segments
                opts.id
                opts.ind
                opts.pat
                opts.row
                
                % parallel computing
                opts.parfor = false
                opts.batch_size = -1
                
                % merging behaviour
                opts.only_seg_in_db = 0
                
            end
            disp('------calc------')
            
            % ----
            % parse inputs
            % ----
            % * parse data_inp
            % ** parse databse
            if nargin < 3
                error('No data inputted.')
            end
            is_database = isa(data_inp,'Seg_db');
            if is_database
                seg_db = data_inp;
                data_inp = seg_db.get_meta('ind',-1);
            else
                seg_db = [];
            end
            
            % ** subsetting input
            if ~isfield(opts,'id') && ~isfield(opts,'ind') ...
                    && ~isfield(opts,'pat') && ~isfield(opts,'row')
                % no subestting -> select all
                ids = select_segs(data_inp,'ind',-1);
            else
                % subsetting
                ids = select_segs(data_inp,'struct',opts);
            end
            n_segs_inp = numel(ids);
            
            % ** compare to ids in internal database and enventually sort out
            if opts.only_seg_in_db
                inters = intersect(self.seg_ids,ids);
                if numel(inters) < numel(ids)
                    disp(['Removed ' num2str(numel(ids)-numel(inters)) ...
                        ' segments from request since they were not in the ' self.name_calc ' database.'])
                end
                ids = inters;
            else
                self.add_segs(ids);
            end
            [~,~,seg_ind_db] = intersect(ids,self.seg_ids,'stable');
            [~,~,seg_ind_input] = intersect(ids,data_inp.segment_id,'stable');
            n_segs_sorted = numel(seg_ind_input);
            
%             ff = table(); % debugging
%             ff.ids = ids;
%             ff.db_ind = self.seg_ids(seg_ind_db);
%             ff.input_ind = data_inp.segment_id(seg_ind_input);
            
            
            % * parse paramsets
            if isempty(paramsets)
                [~,vers_id,~] = self.get_algo_params();
                parset_rows_db = find(strcmp(self.paramset_tbl.vers_id,vers_id));
                disp(['No parameter sets specified. Selected the ' num2str(numel(parset_rows_db)) ' out of '...
                    num2str(size(self.paramset_tbl,1)) ' parametersets with the current algo version.']);
            else
                if ~iscell(paramsets)
                    paramsets = {paramsets};
                end
                parset_rows_db = parse_parsets(self, paramsets,'error');
                not_existent = find(parset_rows_db == -1);
                for i = 1:numel(not_existent)
                    inp = paramsets{not_existent(i)};
                    if isstruct(inp)
                        self.add_paramset(inp)
                        disp(['Added paramset at index ' num2str(not_existent(i))]);
                    else
                        disp(['Couldnt add paramset at index ' num2str(not_existent(i))]);
                    end
                end
            end
            
            % * parse parfor params
            if opts.batch_size > 0 && ~opts.parfor
                warning('No parfor loop executed. Set parfor flag to true.')
            elseif opts.batch_size < 0 && opts.parfor
                opts.batch_size = 1;
            end
            
            
            % ----
            % create tasklist (idea: tasklist ist cell array with 3 columns; 
            % 1st: seg_ind in db, 2nd: seg_ind in data_inp, 3rd: array of 
            % ind of paramsets that need to be calculated
            % ----
            % * reduce to paramsets and segs specified
            calc_status_tocalc = self.calc_status_mat(seg_ind_db,parset_rows_db);
            
            % * create tasklist
            nseg_tocalc = numel(seg_ind_db);
            task_tbl = cell(nseg_tocalc, 3);
            if opts.recalc == true
                for s = 1:nseg_tocalc
                    % add every combination
                    task_tbl{s,1} = seg_ind_db(s);
                    task_tbl{s,2} = seg_ind_input(s);
                    task_tbl{s,3} = parset_rows_db;
                end
            else
                % check calc_status list
                to_remove = false(nseg_tocalc,1);
                for s = 1:nseg_tocalc
                    % find not calculated combinations
                    parsets_to_calc = find(calc_status_tocalc(s,:) == 0);
                    if isempty(parsets_to_calc)
                        to_remove(s) = true;
                    else
                        task_tbl{s,1} = seg_ind_db(s);
                        task_tbl{s,2} = seg_ind_input(s);
                        task_tbl{s,3} = parset_rows_db(parsets_to_calc);
                    end
                end
                task_tbl(to_remove,:) = [];
            end
            
            
            % ----
            % calculating and saving
            % ----
            % * create calc folder if not existent
            parset_ind = self.paramset_tbl.paramset_ind(parset_rows_db);
            for p = 1:numel(parset_ind)
                filepath = [self.files_path self.name_calc '_calc' num2str(parset_ind(p))];
                if ~exist(filepath,'dir')
                    mkdir(filepath);
                end
            end
            
            % * main calc loop
            ntasks = size(task_tbl,1);
            disp(['Start calculation and saving of segments']) 
            tic
            if ~opts.parfor
                
                % loop over tasks
                for s = 1:ntasks
                    disp([num2str(s) '/' num2str(ntasks)])
                    
                    if is_database
                        data = seg_db.get_data('row',task_tbl{s,2});
                    else
                        data = data_inp(task_tbl{s,2},:);
                    end
                    self.exec_calc_task(data,task_tbl(s,:));
%                     self.calc_status_mat() = 1;
                end
                
            else
%                 gcp('nocreate')  % get current parpool
                % create batches
                if ntasks == 0
                    block_ind = [];
                else
                    ind = 1:opts.batch_size:ntasks;
                    block_ind = [ind' [ind(2:end)-1 ntasks]'];
                end
                nbatches = size(block_ind,1);
                task_batches = cell(nbatches,1);
                for b = 1:nbatches  % prebuild batches for less overhead
                    task_batches{b} = task_tbl(block_ind(b,1):block_ind(b,2),:);
                end
                
                % loop over tasks
                parfor b = 1:nbatches
                    
                    % get inputs
                    tasks = task_batches{b};
                    ntasks_in_batch = size(tasks,1);
                    
                    % calculate
                    for t = 1:ntasks_in_batch
                        if is_database
                            data = seg_db.get_data('row',tasks{t,2});
                        else
                            data = data_inp(tasks{t,2},:);
                        end
                        self.exec_calc_task(data,tasks(t,:));
                            % bc of parfor, not possible to change
                            % calc_status_mat
                    end
                end
                
                % change calc_status_mat
                self.calc_status_mat(seg_ind_db,parset_rows_db) = 1;
            end
            disp('Finished calculation and saving of segments.') 
            toc
            
            % update properties
            self.update_props();
            
            % user feedback
            if opts.recalc
                n_calc = numel(calc_status_tocalc);  
                n_retrieved = 0;
            else
                n_retrieved = sum(sum(calc_status_tocalc > 0));
                n_calc = sum(sum(calc_status_tocalc == 0));                
            end
            n_sorted_out = (n_segs_inp - n_segs_sorted) * numel(parset_rows_db);
            disp(['seg-paramset combinations: ' num2str(n_calc) '/' num2str(n_retrieved) '/' num2str(n_sorted_out)...
                ' (calculated/already in database/sorted out)'])          
            disp('----------------')

            % generate output
            if numel(parset_rows_db) == 1
                [tbl_out,parset_out] = self.get(...
                    self.paramset_tbl.paramset_ind(parset_rows_db),...
                    'segs',ids(seg_ind_input));
            end
        end
             
        function [tbl_out,params_full] = get(self,varargin)
            % retrieve data specified by segs and paramset from the
            % database. Specify paramset by either a positional first arg OR name
            % value pairs of parameters. 
            % example calls:
            %   get(parset_struct)
            %   get('wndw_len',30)
            % inputs:
            %   - positional (optional): paramset_ind as number or struct with
            %       params. Not specified params will be set to default in the
            %       struct
            %   - [parameter_names] (nv): give parameter names as name value
            %       pairs. If both parset_struct and name value pairs
            %       given, error.
            %   - segs (nv): cell/string array containing ids or table with 
            %       column segment_id. If segs empty/not given, all segments
            %       returned
            % outputs:
            %   - data_tbl: table with one row per segment inputted,
            %   	olumns contain calculation data
            %   - params_full: struct with complete paramset and comments
            
            % parse inputs
            % * extract inputs with parser
            p = inputParser();
            p.addParameter('segs',self.seg_ids);
            p.addParameter('ignore_missing',1)
            [paramset,p] = self.parse_parset_namevalue(p,varargin);
            ignore_missing = p.Results.ignore_missing;
            segs = p.Results.segs;
            
            % * parse segments
            segs = select_segs(segs,'ind',-1);
            if isempty(segs)
                tbl_out = table();
                params_full = [];
                disp('Specified empty set of segments.')
                return
            end
            if ignore_missing
                [segs,n_seg_req] = self.parse_segs(segs,'warning');
            else
                [segs,n_seg_req] = self.parse_segs(segs,'error');
            end
            
            % * get row indices
            parset_row_ind = self.parse_parsets(paramset,'warning');
            
            
            % main part
            if parset_row_ind > 0
                parset_ind = self.paramset_tbl.paramset_ind(parset_row_ind);
                params_full = table2struct(self.paramset_tbl(parset_row_ind,:));

                % check if all requested data calculated
                [seg_not_calc,ind] = intersect(segs, self.seg_ids(self.calc_status_mat(:,parset_row_ind) == 0),'stable');
                if ~isempty(seg_not_calc)
                    str = ['The following segments have not been calculated: ' self.cell2str(seg_not_calc)];
                    if ignore_missing
                        segs(ind) = [];
                    else
                        error(str)
                    end
                else
                    str = '';
                end

                % load data
                nseg = numel(segs);
                data_tbl = cell(1,nseg);
                for s = 1:nseg
                    seg_ind = find(strcmp(self.seg_ids,segs{s}),1);
                    if self.calc_status_mat(seg_ind,parset_row_ind) == 1
                        data_tbl{s} = self.load_calc(parset_ind,segs{s});
                    elseif self.calc_status_mat(seg_ind,parset_row_ind) == 2
                        data_tbl{s} = self.RAM_cell(seg_ind,parset_row_ind);
                    else
                        error('This should not happen.')
                    end
                end
                data_tbl = vertcat(data_tbl{:});

                % merge with ids
                id_tbl = table(segs,'VariableNames',"segment_id");
                tbl_out = [id_tbl data_tbl];

                % user feedback
                disp([num2str(numel(segs)) ' out of ' num2str(n_seg_req) ' segments retrieved.'])
                warning(str);
            else
                tbl_out = [];
                disp('No data could be retrieved')
            end
             
        end
        
        
        % other utility
        function add_defaultparam(self,param_name,default_value)
            % checks if parameter given by param_name and default_value in
            % get_algo_params but not tracked in table. If so, adds it with
            % its default values to all parametets in the table.
            
            % basically force the user to enter the changes twice 
            default_params = self.get_algo_params();
            if ~isfield(default_params,param_name)
                error('Parameter not found in get_algo_params.')
            end
            if default_params.(param_name) ~= default_value
                error('Value in default parameters and function input are not the same.')
            end
            if ismember(param_name, self.paramset_tbl.Properties.VariableNames)
                error('Parameter already exists in table.')
            end
            
            % add to table
            npar = numel(fieldnames(default_params));
            tbl_row = table(default_value,'VariableNames',string(param_name));
            tbl = table();
            for r = 1:size(self.paramset_tbl,1) % maybe goes easier but this loop works for any value input
                tbl = [tbl; tbl_row];
            end
            self.paramset_tbl = [self.paramset_tbl(:,1:npar) tbl self.paramset_tbl(:,npar+1:end)];
        end
        
        function update_props(self)
            % saves important properties in the info file
            info_file_path = [self.files_path self.name_calc '_info.mat'];
            
            % properties to be saved
            props = {
                'name_calc',...
                'paramset_tbl',...
                'curr_ind_parset',...
                'seg_ids',...
                'calc_status_mat',...
                'vers_lookup'};
            
            % create save struct obj
            for p = 1:numel(props)
                obj.(props{p}) = self.(props{p});
            end
            
            % save
            save(info_file_path,'-struct','obj')
        end

        function params = get_paramset(self,parset_ind)
            % returns complete parameterset as struct with index parset_ind
            
            row_ind = self.parse_parsets(parset_ind);
            if row_ind < 1
                error('Paramset with this row index doesnt exist.')
            end
            params = table2struct(self.paramset_tbl(row_ind,:));
        end
        
    end
    
    methods (Static = true)
        
         function [tbl,parset] = plot(data,varargin)
            % Basically a nice wrapper to call calc_seg with one segment.
            % Takes one table with one row in data specifying one segment
            % and takes one parameter set and calls calc_seg. Main usecases
            % should be viewing the analysis plots or just calculate one
            % segment without the database infrastructure. In the latter
            % case, plotting can be disabled by fig_ind = -1.
            % inputs:
            %   - data: table with one row that specifies all needed data
            %       in columns
            %   - struct (nv): specify a parameterset as struct
            %   - [param_names] (nv,opt): specify parameters by name value
            %       pairs. If parset given, parset will be preferred.
            %   - fig_ind (nv,opt): -1 is no plot, 0 = new figure, 
            %       >0 is the matlab index of the figure to plot in.
            %       Default new figure.
            %   - chn (nv,opt): determines the content to be calculated,
            %       plotted and returned. If 0, whole segment. If chn>0 then
            %       only channel with index chn. Default 0
            % outputs:
            %   - tbl: output tbl as in calc_seg or calc_chn
            %   - parset: structure of full parameterset used in
            %       calculation
            
            % extract subclass string by hand
            ind = find(strcmp(varargin, 'subclass'),1);
            subclass_str = varargin{ind+1};
            
            % simulate static function call by creating empty object of
                % subclass that called this function
            obj = eval([subclass_str '([])']);
            
            % input parser
            p = inputParser();
            p.addParameter('chn',-1);
            p.addParameter('fig_ind',0);
            p.addParameter('subclass_str','Abstract_calc_db');
            [parset,p] = obj.parse_parset_namevalue(p,varargin);
            chn = p.Results.chn;
            fig_ind = p.Results.fig_ind;
            
            % asserting inputs
            if ~istable(data)
                error('First input argument has to be table.')
            end
            if size(data,1) ~= 1
                error('Only exactly one row in data input allowed.')
            end
            
            % open correct plot (or none)
            if fig_ind == 0
                fig_hndl = figure();
            elseif fig_ind > 0
                fig_hndl = figure(fig_ind);
            else
                fig_hndl = [];
            end
            
            % call calc_seg function from the current subclass
            if chn > 0
                 if 1 > chn || chn > size(data.segment_data{1},1)
                    error('Channel out of bound for input data!')
                end
                tbl = obj.calc_chn(data,parset,chn,fig_hndl);
            else
                tbl = obj.calc_seg(data,parset,fig_hndl);
            end
         end
        
    end
    
    methods (Access = 'protected')
       
        function params_out = parse_param_struct(self,params)
            % takes struct, extracts specfied params and adds non
            % specified params with default values. Returns struct with
            % exactly one field for each parameter and one field for
            % current algo version id.
            % inputs:
            %   - params: any struct
            
            default_params = self.get_algo_params();
            [~,vers_id,~] = self.get_algo_params();
            default_params.vers_id = vers_id;
            default = fieldnames(default_params);
            inp = fieldnames(params);
            params_out = params;
            
            % Remove fields from input with no default params
            ind_keep = ismember(inp,default);
            rm = inp(~ind_keep);
            if ~isempty(rm)
                params_out = rmfield(params_out,rm);
                warning(['Removed fields from input: ' self.cell2str(rm)])
            end
            
            % fill not specified fields with defaults
            ind_default = find(~ismember(default,inp));
            for i = 1:numel(ind_default)
                p = ind_default(i);
                params_out.(default{p}) = default_params.(default{p});
            end
            if ~isempty(ind_default)
%                 warning(['Added fields to input: ' self.cell2str(default(ind_default))])
            end
            
            % reorder
            params_out = orderfields(params_out,default_params);
        end
        
        
        function ind_out = get_parset_row_ind(self,paramset)
            % checks if paramset given as argument existent in databse.
            % Gives specific feedback if only versions mismatch.
            % input:
            %   - paramset: paramset_ind as number or struct with
            %   params. Not specified params will be set to default in the
            %   struct
            % output:
            %   - ind_out: If -1, not contained in db. If 0, parameters 
            %   contained but different algo version. If >0, row_ind of parset.
            
            % parse number and struct version into row_ind
            if isnumeric(paramset)
                parset_ind = paramset;
                [~,parset_row_ind_alleq] = ismember(parset_ind,self.paramset_tbl.paramset_ind);
                parset_params_eq = 0;
            else
                params = self.parse_param_struct(paramset);
                params = struct2table(params,'AsArray',true);
                

%                 [parset_params_eq,~] = ismember(params(1,1:end-1),self.paramset_tbl(:,2:end-3));
%                 [~,parset_row_ind_alleq] = ismember(params,self.paramset_tbl(:,2:end-2));
                parset_params_eq = self.ismember_cell(params(1,1:end-1),self.paramset_tbl(:,2:end-3));
                parset_row_ind_alleq = self.ismember_cell(params,self.paramset_tbl(:,2:end-2));
            end
            
            ind_out = (parset_row_ind_alleq == 0 & parset_params_eq == 0) *(-1) ...
                + parset_row_ind_alleq;
        end
          
        function parset_row_ind = parse_parsets(self, parsets, warnerror)
            % Parses any paramsetinput that is not name value pair
            % basically loop over get_parset_row_ind that throws errors
            % if conflicting versions (0) or not found in database (-1)
            
            % parse error warning
            if nargin < 3 || strcmp(warnerror,'error') == 1
                is_error = 1;
            else
                is_error = 0;
            end       

            
            % check if array then transform to cell
            if isnumeric(parsets)
                cell_parsets = cell(size(parsets));
                for i = 1:numel(parsets)
                    cell_parsets{i} = parsets(i);
                end
                parsets = cell_parsets;
            elseif ~iscell(parsets)
                % check if user didnt to put cell for single paramset input
                parsets = {parsets};
            end
            

            % check if parsets exist in database
            parset_row_ind = nan(numel(parsets),1);
            for p = 1:numel(parsets)
                parset_row_ind(p) = get_parset_row_ind(self,parsets{p});
                if parset_row_ind(p) == -1
                    str = ['The inputted paramset at index '...
                        num2str(p) ' does not exist in database'];
                    if is_error
                        error(str)
                    else 
                        warning(str)
                    end
                elseif parset_row_ind(p) == 0
                    req_vers = parsets{p}.vers_id;
                    [~,curr_vers,~] = self.get_algo_params();
                    str = ['Conflicting algo versions at cell index ' num2str(p) ...
                        '(requested version: ' req_vers ', current version: ' curr_vers];
                    if is_error
                        error(str)
                    else 
                        warning(str)
                    end
                end
            end
        end

        function [params,p] = parse_parset_namevalue(self,p,args)
            % parses struct and args into a paramstruct. If both struct and
            % args nonempty, struct prefered. If none given, default
            % arguments taken. p is parser with other arguments from function.
            
            % parse inputs
            p.StructExpand = 0;
            [default_params,vers_id] = self.get_algo_params();
            default_params.vers_id = vers_id;
            par_names = fieldnames(default_params);
            p.addOptional('struct',[]);
            for par = 1:numel(par_names)
                p.addParameter(par_names{par},default_params.(par_names{par}));
            end
%             p.addParameter('struct',[]);
            parse(p,args{:});
            
            % extract struct or params
            if ~isempty(p.Results.struct) 
                % struct (really ugly, as gets row_ind, translates it to
                % param struct just to get translated back to row_ind in
                % get)
                struct = p.Results.struct;
                if isnumeric(struct)
                    row_ind = self.parse_parsets(struct,'warning');
                    params = table2struct(self.paramset_tbl(row_ind,2:end-2));
                else
                    params = self.parse_param_struct(struct);
                end
            else
                % name value pair or default
                % extract the paarams from other potential arguments
                clear params;
                for par = 1:numel(par_names)
                    params.(par_names{par}) = p.Results.(par_names{par});
                end
                params.vers_id = p.Results.vers_id;
            end
            
        end
             
        function [segs_out,n_segs_in] = parse_segs(self,segs,str)
            % parses segs input, checks if all segs exist in database and
            % eventually throws warning/error
            
            % parse input into string array
            segs = select_segs(segs,'ind',-1);
            n_segs_in = numel(segs);
            
            % check which are existent
            in_db = ismember(segs,self.seg_ids);
            seg_not_in_db = segs(~in_db);
            if ~isempty(seg_not_in_db)
                if strcmp(str,'warning')
                    warning(['The following segments are not in database and have been removed ' ...
                    'from the request: ' self.cell2str(seg_not_in_db)]);
                    segs(~in_db) = [];
                else
                    error(['The following segments are not in database: ' self.cell2str(seg_not_in_db)]);
                end
            end
            segs_out = segs;
        end
        
        
        function save_calc(self,parset_ind,seg_id,data_tbl)
            % saves the given data in the current data structure
            % inputs: 
            %   - parset_ind: index of parameter set 
            %   - seg_id: cell/string array with seg_ids
            %   - data_tbl: calculation results to be saved as table 
            
            filename = [self.files_path self.name_calc '_calc' num2str(parset_ind) ...
                '/' seg_id '.mat'];
            tbl = data_tbl;
            save(filename,'tbl')
        end
        
        function tbl_out = load_calc(self,parset_ind,seg_id)
            % loads the specified data from data structure
            % inputs: 
            %   - parset_ind: index of parameter set 
            %   - seg_ids: seg_id
            % outputs:
            %   - tbl_out: tbl with calculation results in columns, one
            %   column for seg_ids and the values for the specified
            %   segments in the rows
            
            % main loop
            filename = [self.files_path self.name_calc '_calc' num2str(parset_ind) ...
                '/' seg_id '.mat'];
            tbl_out = load(filename).tbl;
        end
        
        function ret = del_calc_data(self,seg_id,parset_ind)
            % deletes the specified segment parameter set combination
            % from the database
            % inputs:
            %   - parset_ind: index of parameter set 
            %   - seg_id: cell/string array with seg_ids
            % outputs:
            %   - ret: returns 1 if successful and 0 if not.
            
            filename = [self.files_path self.name_calc '_calc' num2str(parset_ind) ...
                '/' seg_id '.mat'];
            try
                delete(filename);
                ret = 1;
            catch
                warning(['Couldnt delete seg ' seg_id ' at parset ' num2str(parset_ind)])
                ret = 0;
            end
            
        end
        
        function exec_calc_task(self,data,task_tbl)
            % executes a task, i.e. calculate multiple parametersets for
            % one segment and save the results
            
            % parse inputs
            paramsets_row_ind = task_tbl{1,3};
            seg_ind = task_tbl{1,1};
            
            nparsets = numel(paramsets_row_ind);
            for p = 1:nparsets
                
                % calc
                parset = self.paramset_tbl(paramsets_row_ind(p),1:end-3);
                tbl = self.calc_seg(data,table2struct(parset(1,2:end)),[]);
                
                % save
                self.save_calc(parset.paramset_ind,data.segment_id{1},tbl);
                
                % change calc_tbl status
                self.calc_status_mat(seg_ind,paramsets_row_ind) = 1;
            end
        end
        
        
        function add_curr_vers(self)
            % adds the current version of the algorithm to the lookup tbl.
            % If already existent, comment renewed.
            
            [~,vers_id,vers_cmm] = self.get_algo_params();
            vers_tbl = self.vers_lookup;
            if isempty(vers_tbl)
                row = 0;
            else
                [~,row] = ismember(vers_id,vers_tbl(:,1));
            end
            if row == 0
                self.vers_lookup = [[{vers_id},{vers_cmm}]; self.vers_lookup];
            else
                self.vers_lookup{row,2} = vers_cmm;
            end
        
        end
        
        function cmm = lookup_vers_comment(self,vers_id)
            % returns the version comment if version recorded before
            
            [~,row] = ismember(vers_id,self.vers_lookup(:,1));
            if row > 0 
                cmm = self.vers_lookup{row,2};
            else
                cmm = '';
            end
        end
        
        function row_ind = ismember_cell(self,row,tbl)
            % returns the row index of row in table if existent. If not
            % existent return 0. Compares cells with isequal
            
            % identify cell columns
            if isempty(tbl)
                row_ind = 0;
                return
            end
            ncols = size(tbl,2);
%             is_cell = zeros(ncols,1);
%             for c = 1:ncols
%                 is_cell(c) = iscell(tbl{1,c});
%             end
%             
%             % compare with removing cell columns
%             tbl_ind = ismember(tbl(:,~is_cell),row(~is_cell));

            % compare just cell columns separately
%             ncellcols = sum(is_cell);
            nrows = size(tbl,1);
            vals = zeros(nrows,ncols);
            
            % loop over all rows
            for c = 1:ncols
                for r = 1:nrows
                    e1 = row{1,c};
                    e2 = tbl{r,c};
                    if isnumeric(e1) && iscell(e2)
                        e1 = {e1};
                    elseif isnumeric(e2) && iscell(e1)
                        e2 = {d2};
                    end
                    vals(r,c) = isequal(e1,e2);
                end
            end
            
            % find row where all equal
            row_ind = find(all(vals,2));
            if isempty(row_ind)
                row_ind = 0;
            end
            
            
            
            
            
            % OLD VERSION
            
%             for c = 1:ncols
%                 for r = 1:nrows
%                     e1 = row{1,c};
%                     e2 = tbl{r,c};
%                     vals(r,c) = isequal(e1,e2);
%                     % account for weird matlab behaviour where 1-dim arrays
%                     % are not saved as cells in tables where higher
%                     % dimensional arrays are
%                     if ~iscell(e1) && iscell(e2)
%                         vals(r,c) = isequal({e1},e2);
%                     elseif iscell(e1) && ~iscell(e2)
%                         vals(r,c) = isequal(e1,{e2});
%                     else
%                         vals(r,c) = isequal(e1,e2);
%                     end
%                 end
%             end
%             row_ind = find(all(vals,2));
%             if isempty(row_ind)
%                 row_ind = 0;
%             end
%             
%             % merge both comparisons
            
        end
                
        function str = cell2str(self,c)
            % creates string from cell array
            str = '';
            for i = 1:numel(c)
                str = [str c{i} ',' ];
            end
            str(end) = '';
        end
        

    end
   
end