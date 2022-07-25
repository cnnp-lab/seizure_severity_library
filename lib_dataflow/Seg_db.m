classdef Seg_db < handle
    % defines a classical database for the seizure data and metadata. As of
    % now, only functions for minimum functionality are defined. Allows to
    % easily handle seizure data from multiple seizures and patients.
    
    % potential improvements:
    % - substitute through actual SQL database
    % - current structure saves the names of ext_data columns in each
    % seizure field -> unnecessary storage. Maybe other possibility
    % (tall tables)?
    % - handle missing values for heavy data (as of now, error will be
    % thrown)
    
    properties
        files_path      % path at which data heavy columns are stored (with '/' in end!)
        Seg_meta         % table with metadata about szs stored
        ext_data        % names of columns to be stored externally as cell/string array. 
                        % Should not be changed after creation.
        ind_curr        % denotes db_ind of last seizure added
    end
    
    methods 
        function obj = Seg_db(files_path,opts)
            % constructor
            arguments
                files_path
                opts.ext_data
            end
            
            % files_path
            files_path(ismember(files_path,'\')) = '/';
            if ~strcmp(files_path(end),'/')
                files_path = [files_path '/'];
                warning('"/" added to filepath!')
            end
            obj.files_path = files_path;
            

            % see if database already exists
            try
                d = load([files_path 'Sz_db_meta.mat']);
                obj.Seg_meta = d.Seg_meta;
                obj.ind_curr = d.ind_curr;
                obj.ext_data = d.ext_data;
                disp('Connected to database.')
            catch
                % create folder if not existent
                if ~exist(files_path,'dir')
                    mkdir(files_path)
                    disp('Specified directory created')
                end
                
                % create database
                Seg_meta = table();
                ind_curr = 0;
                if nargin > 1 && isfield(opts,'ext_data')
                    ext_data = opts.ext_data;
                else
                    ext_data = {'segment_data'};
                end
                save([files_path 'Sz_db_meta.mat'], 'Seg_meta','ind_curr','ext_data');
                
                obj.Seg_meta = Seg_meta;
                obj.ind_curr = 0;
                obj.ext_data = ext_data;
                disp('Created new database at specified location.')
            end
            

        end
        
        function add(self, Sz_inp)
            % adds the data of Sz_inp (table) into the database. Columns are not
            % checked for completeness in both directions, i.e.
            % non-existing columns in database are created and filled with
            % nans for other seizures and non-existing columns in Sz_inp
            % are created and filled with nans. 
            
            % check if input unique
            if sum(strcmp(Sz_inp.Properties.VariableNames,'segment_id')) == 0
                error('Sz_inp has no column "segment_id"!')
            end
            if sum(isempty(Sz_inp.segment_id)) > 0
                error('Sz_inp has segments without segment_id!')
            end
            if numel(unique(Sz_inp.segment_id)) ~= size(Sz_inp,1)
                error('Sz_inp has not unique set of segment_ids!')
            end
            
            % remove szs from Sz_inp already existing in Seg_meta
            if self.ind_curr > 0
                nsz_inp = size(Sz_inp,1);
                del_ind = zeros(nsz_inp,1);
                for s = 1:nsz_inp
                    if sum(strcmp(self.Seg_meta.segment_id,Sz_inp.segment_id{s}))>0
                        del_ind(s) = 1;
                    end 
                end
                if sum(del_ind) > 0
                    warning(['The following ' num2str(sum(del_ind)) ' segment(s) already exist in the database and were not added: '...
                        self.cell2str(Sz_inp.segment_id(logical(del_ind)))])
                end
                Sz_inp(logical(del_ind),:) = [];
            end
                
            % extract meta data
            ind_tbl = table((self.ind_curr+1:self.ind_curr+size(Sz_inp,1))', 'VariableNames',{'db_ind'});
            data_ind = [];
            for i = 1:numel(self.ext_data)
                tmp = find(strcmp(Sz_inp.Properties.VariableNames,self.ext_data{i}));
                if isempty(tmp)
                    error(['Data column "' self.ext_data{i} '" not found in input data.'])
                end
                data_ind = [data_ind tmp];
            end
            ind_meta = setdiff(1:size(Sz_inp,2),data_ind);
            Sz_prep = [ind_tbl Sz_inp(:,ind_meta)];
            
            % add metadata to db
            if isempty(self.Seg_meta)
                self.Seg_meta = Sz_prep;
            else
                % merge Seg_meta and Sz_prep by creating missing columns
                % * find missing columns
                add_to_meta = setdiff(Sz_prep.Properties.VariableNames, self.Seg_meta.Properties.VariableNames);
                add_to_prep = setdiff(self.Seg_meta.Properties.VariableNames, Sz_prep.Properties.VariableNames);
                
                % * create missing columns
                nprep = size(Sz_prep,1);
                for  c_name = add_to_meta
    
                    % choose correct datatype to fill missing values with
                    data_type = Sz_prep.(c_name{1})(1);
                    if iscategorical(data_type)
                        self.Seg_meta.(c_name{1}) = categorical(nan(size(self.Seg_meta,1),1));
                    elseif iscell(data_type)
                        self.Seg_meta.(c_name{1}) = cell(size(self.Seg_meta,1),1);
                    else
                        data_type = standardizeMissing(data_type,data_type);
                        self.Seg_meta.(c_name{1}) = repmat(data_type,size(self.Seg_meta,1),1);
                    end
                    warning(['Column "' num2str(c_name{1}) '" added to db!'])
                end
                for  c_name = add_to_prep
                    
                    % choose correct datatype to fill missing values with
                    data_type = self.Seg_meta.(c_name{1})(1);
                    if iscategorical(data_type)
                        Sz_prep.(c_name{1}) = categorical(nan(size(Sz_prep,1),1));
                    elseif iscell(data_type)
                        Sz_prep.(c_name{1}) = cell(size(Sz_prep,1),1);
                    else
                        data_type = standardizeMissing(data_type,data_type);
                        Sz_prep.(c_name{1}) = repmat(data_type,nprep,1);
                    end
                    warning(['Column "' num2str(c_name{1}) '" of new seizures filled with nans!'])
                end
                
                % * concat
                self.Seg_meta = vertcat(self.Seg_meta, Sz_prep);
            end
            self.ind_curr = self.ind_curr + size(Sz_inp,1);
            
            % save updated meta data
            self.update_props();
            
            % save heavy data
            for s = 1:size(Sz_inp,1)    % probably possible without loop
                
                % create patient folder
                filename = [self.files_path Sz_inp.patient_id{s} '/'];
                if ~exist(filename,'dir')
                    mkdir(filename)
                end
                
                % save data
                tbl = Sz_inp(s,data_ind);
                save([filename  Sz_inp.segment_id{s} '.mat'],'tbl','-v7.3');
            end
            
            % display feedback for user
            disp(['Added ' num2str(size(Sz_inp,1)) ' segment(s)'])
               
        end
        
        
        function seg_out = get_data(self, opts)
            % returns table with data from segments specified in
            % ind_struct.
            arguments
                self

                % indexing
                opts.id
                opts.ind
                opts.pat
                opts.row
            end
            
            [~,ind] = select_segs(self.Seg_meta,'struct',opts);
            seg_out = self.Seg_meta(ind,:);
            
            % get data
            data_tbl = cell(1,numel(self.ext_data));
            for s = 1:numel(ind)
%                     filename = [self.files_path self.Seg_meta.patient_id{ind(s)} '.mat'];
%                     data_tbl{s} = load(filename,self.Seg_meta.segment_id{ind(s)});
                filename = [self.files_path self.Seg_meta.patient_id{ind(s)} '/' ...
                    self.Seg_meta.segment_id{ind(s)} '.mat'];
                data_tbl{s} = load(filename).tbl;
            end
            data_tbl = vertcat(data_tbl{:});
            seg_out = [seg_out, data_tbl];
            
            
        end
        
        function seg_out = get_meta(self,opts)
            % returns table with data from segments specified in
            % ind_struct. No external data retrieved from database. Much
            % faster than get_data.
            arguments
                self

                % indexing
                opts.id
                opts.ind
                opts.pat
                opts.row
            end
            
            [~,ind] = select_segs(self.Seg_meta,'struct',opts);
            seg_out = self.Seg_meta(ind,:);
        end
        
        function del(self,opts)
            % deletes segments specified by standard specification

            arguments
                self
                
                opts.id
                opts.ind
                opts.pat
                opts.row
            end
            
            
            % parse indices
            [~,ind] = select_segs(self.Seg_meta, 'struct',opts);
            
            % delete data
            for s = 1:numel(ind)
                folder_loc = [self.files_path self.Seg_meta.patient_id{ind(s)}];
                filename = [folder_loc '/' self.Seg_meta.segment_id{ind(s)} '.mat'];
                delete(filename);
                
                % delete folder location
                if numel(dir(folder_loc)) == 2
                    rmdir(folder_loc)
                end
            end
                     
            % delete seg from meta
            self.Seg_meta(ind,:) = [];
            
            % save props
            self.update_props();
            disp(['Deleted ' num2str(numel(ind)) ' segment(s).']);
            
        end
        
        function refactor(self)
            % refactors the db_ind to the row indices of the current Seg_meta 
            % table. This will delete gaps in the db_ind and thus 
            % WILL change the db_ind of some seizures!
            
            inp = input('Refactoring will change db_ind. Are you sure you want to continue? ["y"/"n"] ');
            if strcmp(inp,'y')
                if size(self.Seg_meta,1) > 0
                    self.Seg_meta.db_ind = (1:size(self.Seg_meta,1))';
                    self.ind_curr = size(self.Seg_meta,1);
                end
                filepath = [self.files_path 'Sz_db_meta.mat'];
                d = load(filepath);
                d.Seg_meta = self.Seg_meta;
                d.ind_curr = self.ind_curr;
                save(filepath,'-struct','d')
                warning('Db_ind refactored!')
            end
        end
        
        function update_props(self)
            Seg_meta = self.Seg_meta;
            ind_curr = self.ind_curr;
            ext_data = self.ext_data;
            save([self.files_path 'Sz_db_meta.mat'],'Seg_meta','ind_curr','ext_data')
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