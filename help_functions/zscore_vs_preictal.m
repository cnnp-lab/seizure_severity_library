function zval = zscore_vs_preictal(data,pre_ind,opts)
% data has maximum 2 dimensions, pre_ind denotes index of pre, along dimension dim,
% default is 2 (i.e. along the columns). If column vector with no
% dimension inputted, then algo will zscore along rows

    arguments
        data
        pre_ind
        
        opts.is_mad = 0
        opts.thrhld_outlier = 0.025
        opts.dim = 2
    end
    
    % set data to always have operation dimension in rows
    isColVec = size(data,2) == 1;        
    if isColVec || opts.dim == 1
        data = data';
    end
    
    
    
    
    data_pre = data(:,1:pre_ind);
    
    ind_inlier = ~isoutlier(data_pre,'percentiles',[opts.thrhld_outlier*100 (1-opts.thrhld_outlier)*100],2);
    tmp = nan(size(data_pre));
    tmp(ind_inlier) = data_pre(ind_inlier);
    data_pre = tmp;
    
    
    if opts.is_mad
        mean_norm = median(data_pre,2,'omitnan');
        std_norm = mad(data_pre,[],2);  %'omitnan' already in method
    else
        mean_norm = mean(data_pre,2,'omitnan');
        std_norm = std(data_pre,[],2,'omitnan');
    end
        
    zval = ( data - mean_norm ) ./ std_norm;
    if numel(find(std_norm == Inf)) > 0
        warning(['indizes ' num2str(find(std_norm == Inf)) ' have 0 std!'])
    end
    
    
    
    
    
    
    if isColVec || opts.dim == 1
        zval = zval';
    end
    
end