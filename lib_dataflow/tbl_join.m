function tbl_out = tbl_join(tbl1,tbl2)
    % joins together the 2 tables tbl1,tbl2.
    % inputs:
    %   - tbl1,tbl2: tables to join with at least column segment_id
    
    arguments
        tbl1
        tbl2
    end
    
    % actual innerjoin
    names1 = tbl1.Properties.VariableNames;
    names2 = tbl2.Properties.VariableNames;
    db_ind = ismember('db_ind',names1) && ismember('db_ind',names2);
    patient_id = ismember('patient_id',names1) && ismember('patient_id',names2);
    if db_ind && patient_id
        tbl_out = innerjoin(tbl1,tbl2,'keys',{'segment_id','db_ind','patient_id'});
    elseif db_ind && ~patient_id
        tbl_out = innerjoin(tbl1,tbl2,'keys',{'segment_id','db_ind'});
    elseif ~db_ind && patient_id
        tbl_out = innerjoin(tbl1,tbl2,'keys',{'segment_id','patient_id'});
    else
        tbl_out = innerjoin(tbl1,tbl2,'keys',{'segment_id'});
    end
    
    
    % user feedback
    inters = intersect(names1,names2);
    inters = setdiff(inters,{'segment_id','db_ind','patient_id'}); % remove common keys
    if numel(inters) > 0
        warning(['Following columns conained in both tables: ' cell2str(inters)])
    end
        
    % internal functions
    function str = cell2str(c)
        % creates string from cell array
        str = '';
        for i = 1:numel(c)
            str = [str c{i} ',' ];
        end
        str(end) = '';
    end
        
end