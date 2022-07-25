function [seg_ids,row_ind] = select_segs(seg_inp, opts)
% takes multiple segments in various data representations and allows selecting
% specific segments by segment_id, patient_id, database_ind or row index 
% (last one for convenience such that everything in single function).
% Returns the row index and the segment_ids of selected segments.
% Focus on intuitive indexing, not performance. 
% Inputs: 
%   - seg_inp: one of the following: chars,string,string/cell array,table,Sz_db
%   - opts.id: single string or string array or cell array specifiying
%   seizures by their segment_id. If opts.id == -1, then all seizures
%   specified
%   - opts.ind: single number or array specifiying seizures by their
%   db_ind
%   - opts.pat: single string or string array or cell array containint
%   patient_ids that specify all seizures of those patients
%   - opts.struct: struct that might contain the fields id,ind,pat handled
%   as above. If both struct and other optionals given, opts.id and
%   opts.struct.id will be concatenated together.
% Outputs:
%   - row_ind: array with indices of specified segments in seg_inp
%   - seg_ids: cell array with segment ids of specified segments


arguments
    % lookup table/array/database
    seg_inp
    
    % name value pairs
    opts.id
    opts.ind
    opts.pat
    opts.row
    
    % struct with fields id,ind,pat
    opts.struct
end
if isempty(seg_inp)
    seg_ids = [];
    row_ind = [];
    return
end

% parse segment input to table
% * parse chars to string
if ischar(seg_inp)
    seg_inp = {seg_inp};
end
    
% * cell and string arrays
if iscell(seg_inp) || isstring(seg_inp)
    % ceck for one dim
    if numel(seg_inp) ~= max(size(seg_inp))
        error('String/cell array must only have one dimension!')
    end
    seg_inp = table(reshape(seg_inp,[],1),'VariableNames',{'segment_id'});
end

% * sz_db
if isa(seg_inp,'Seg_db')
    seg_inp = seg_inp.get_meta('ind',-1); 
end


% parse subsetting inputs into id,ind,pat
id = [];
ind = [];
pat = [];
row = [];
if isfield(opts,'id')
    if ischar(opts.id)
        opts.id = {opts.id};
    end
    id = opts.id(:);
end
if isfield(opts,'ind')
    ind = opts.ind(:);
end
if isfield(opts,'pat')
    if ischar(opts.pat)
        opts.pat = {opts.pat};
    end
    pat = opts.pat(:);  
end
if isfield(opts,'row')
    row = opts.row(:);
end
if isfield(opts,'struct')   % could probably be done easier
    if isfield(opts.struct,'id')
        if ischar(opts.struct.id)
            opts.struct.id = {opts.struct.id};
        end
        id = [id; opts.struct.id(:)];
    end
    if isfield(opts.struct,'ind')
        ind = [ind; opts.struct.ind(:)];
    end
    if isfield(opts.struct,'pat')
        if ischar(opts.struct.pat)
            opts.struct.pat = {opts.struct.pat};
        end
        pat = [pat; opts.struct.pat(:)];
    end
    if isfield(opts.struct,'row')
        row = [row; opts.struct.row(:)];
    end
end


% get row ind
if ind == -1
    row_ind = 1:size(seg_inp,1);
    seg_ids = seg_inp.segment_id(row_ind);
else
    
    % get row ind by id
    row_ind = [];
    if ~isempty(id)
        if ~ismember('segment_id',seg_inp.Properties.VariableNames)
            error('Segment ids specified but no column "segment_id" in input table.')
        end
        [idx,~] = ismember(seg_inp.segment_id,id);
        [~,check] = ismember(id,seg_inp.segment_id);
        if sum(check == 0) > 0 
            not_found = id(check==0);
            warning(['The following seg_ids were not found: ' cell2str(not_found)]);
        end
        row_ind = [row_ind; find(idx)];
    end

    % get row ind by ind
    if ~isempty(ind)
        if ~ismember('db_ind',seg_inp.Properties.VariableNames)
            error('Database indices specified but no column "db_ind" in input table.')
        end
        [idx,~] = ismember(seg_inp.db_ind,ind);
        [~,check] = ismember(ind,seg_inp.db_ind);
        if sum(check == 0) > 0 
            not_found = ind(check==0);
            warning(['The following db_inds were not found: ' numarray2str(not_found)]);
        end
        row_ind = [row_ind; find(idx)]; 
    end

    % get row ind by pat
    if ~isempty(pat)
        if ~ismember('patient_id',seg_inp.Properties.VariableNames)
            error('Patients specified but no column "patient_id" in input table.')
        end
        [idx,~] = ismember(seg_inp.patient_id,pat);
        [~,check] = ismember(pat,seg_inp.patient_id);
        if sum(check == 0) > 0 
            not_found = pat(check==0);
            warning(['The following patients were not found: ' cell2str(not_found)]);
        end
        row_ind = [row_ind; find(idx)];
    end
    
    % get row ind by ind
    if ~isempty(row)
        check =  1 <= row & row <= size(seg_inp,1);
        if sum(check == 0) > 0 
            not_found = row(check==0);
            warning(['The following rows were not found: ' numarray2str(not_found)]);
        end
        row_ind = [row_ind; row(check)]; 
    end

    % merge all together
    row_ind = unique(row_ind);
    seg_ids = seg_inp.segment_id(row_ind);
end


 % local function for nice warning/error messages
    function str = cell2str(c)
       str = '';
       for i = 1:numel(c)
           str = [str c{i} ',' ];
       end
       str(end) = '';
    end
    function str = numarray2str(a)
       str = '';
       for i = 1:numel(a)
           str = [str num2str(a(i)) ',' ];
       end
       str(end) = '';
    end

end