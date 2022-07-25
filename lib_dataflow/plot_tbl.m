function plot_tbl(tbl,rows,cols,opts)
    % function that visualises the interaction between specified columns 
    % in tbl with subplot tables. Depending on the datatype, either
    % correlation plot, boxplot, or table plotted in subplot.
    % Function is basically a wrapper for the gramm package 
    % https://github.com/piermorel/gramm, official citation:
    %   Morel, (2018). Gramm: grammar of graphics plotting in Matlab. 
    %   Journal of Open Source Software, 3(23), 568, 
    %   https://doi.org/10.21105/joss.00568
    % For installation of the package see github page (basically clone
    % repository and add the cloned folder 'gramm' (NOT @gramm) to the
    % matlab path.
    
    % input: columns are specified by one of the following: column index,
    % column name as char or string, string array, numeric array, or cell array of
    % combinations of the first two. If either cols or rows empty, they
    % will be filled with input from the other.
    %   - tbl: table with data to visualize
    %   - rows: data that should be plotted in rows 
    %   - cols: data that should be plotted in cols 
    %   - opts.color: color all plots by column interpreted as categorical data
    %   - opts.adv_corrplt: If true, takes all inputs from rows
    %       and plots nice dense correlation plot. Requires Econometrics
    %       toolbox! Doesnt work with opts.color argument. If not given,
    %       program will decide whats best.
    %   - opts.adv_corrplt_names: adv_corr_plot will shorten column names to
    %       5 chars. Give optional new names of all columns in rows here.
    %   - opts.fig_ind: 0 = new figure, 0.5 = into current figure, >=1 into
    %       figure with that figure index.
    
    % potential improvements: 
    %   - add optional to highlight specific obersvations with triangles

    
    arguments
        tbl
        rows
        cols
        
        opts.color
        opts.catcont = 'boxplot'   % jitter, violin
        opts.fig_ind = 0
        opts.adv_corrplt = 0   % default: let program decide
        opts.adv_corrplt_names
    end
    
    % find data specified
    var_names = tbl.Properties.VariableNames;
    row_ind = find_col_ind(var_names,rows);
    col_ind = find_col_ind(var_names,cols);
    
    % fill either rows or cols empty
    if isempty(row_ind) && isempty(col_ind)
        error('Specify either col or row ind!')
    elseif isempty(row_ind)
        row_ind = col_ind;
    elseif isempty(col_ind)
        col_ind = row_ind;
    end
    
    % parse color
    if isfield(opts,'color')
        ind = find_col_ind(var_names,opts.color);
        color = tbl{:,ind};
    else
        color = [];
    end

    % automatically check if adv_corrplt better
    if ~isfield(opts,'adv_corrplt')
        
        % check that no colour given
        auto_choice = 1;
        if ~isempty(color)
            auto_choice = 0;
        end
        
        % check if identical elements
        if numel(row_ind) ~= numel(col_ind)
            auto_choice = 0;
        elseif sum(sort(row_ind)==sort(col_ind)) ~= numel(unique([row_ind col_ind]))
            auto_choice = 0;
        end
            
        % check if all data numeric
        ind = row_ind;
        for i = ind
            if ~isnumeric(tbl{:,i})
                auto_choice = 0;
            end
        end
        
        % overwrite
        opts.adv_corrplt = auto_choice;
        if auto_choice
            disp('Using advanced correlation plot. Consider using argument "adv_corrplt_names" for better readability.')
        end
    end
    
    
    % open figure
    if opts.fig_ind >= 1
        figure(opts.fig_ind)
        clf(opts.fig_ind)
    elseif opts.fig_ind == 0
        figure()
    end
    
    % plotting
    if opts.adv_corrplt
        
        % check if all numeric
        ind = row_ind;
        for i = ind
            if ~isnumeric(tbl{:,i})
                error('All inputs have to be numeric!')
            end
        end
        
        % plot better correlation plots
        if isfield(opts,'adv_corrplt_names')
            [~, ~, h] = corrplot(tbl(:,ind),'varNames',opts.adv_corrplt_names);
        else
            [~, ~, h] = corrplot(tbl(:,ind));
        end
        
        % make nice axes (copied from stackoverflow 1:1)
        lineHandles = h(strcmp(get(h, 'type'), 'line'));       %get handles for scatter plots only
        % Loop through each scatter plot
        for i = 1:numel(lineHandles)
            x = lineHandles(i).XData;                         %x data 
            y = lineHandles(i).YData;                         %y data
            xlim(lineHandles(i).Parent, [min(x), max(x)]);    % set x limit to range of x data
            ylim(lineHandles(i).Parent, [min(y), max(y)]);    % set y limit to range of y data

            % To convince yourself that the axis scales are still the same within rows/cols,
            % include these two lines of code that will display tick marks.
            %lineHandles(i).Parent.Position(3:4) = lineHandles(i).Parent.Position(3:4) * .8; 
            %set(lineHandles(i).Parent, 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
        end
        % now take care of the x axis limits of the histogram plots
        histHandles = h(strcmp(get(h, 'type'), 'histogram'));     %handles to all hist plots
        % loop through hist plots
        for j = 1:numel(histHandles)
            x = histHandles(j).BinEdges;                         %bin edges
            xlim(histHandles(j).Parent, [min(x), max(x)]);       %set x limits
        end
        
    else
    
        % plot normal plots
        ncols = numel(col_ind);
        nrows = numel(row_ind);
        clear g
        for c = 1:ncols
            for r = 1:nrows

                % find out what data to plot
                x = tbl{:,col_ind(c)};
                y = tbl{:,row_ind(r)};

                % remove nans
                rm_ind = false(size(x));
                if isnumeric(x)
                    rm_ind = rm_ind | isnan(x);
                end
                if isnumeric(y)
                    rm_ind = rm_ind | isnan(y);
                end
                if sum(rm_ind) > 0 
                    x(rm_ind) = [];
                    y(rm_ind) = [];
                    warning('Removed nans.')
                end

                % invalid input handling
                if ~isnumeric(x) && ~iscategorical(x)
                    error(['Invalid data type at ' var_names{col_ind(c)}]);
                end
                if ~isnumeric(y) && ~iscategorical(y)
                    error(['Invalid data type at ' var_names{row_ind(r)}]);
                end

                % main plot function
                if isnumeric(x) && isnumeric(y)
                    % cont vs cont -> correlation plot

                    % plot points
                    p = gramm('x',x,'y',y,'color',color);
                    p.geom_point();

                    % plot correlation line as linear regression
                    beta = polyfit(x,y,1);
                    corr_fit = @(x) beta(2) + beta(1)*x;
                    p.geom_funline("fun",corr_fit);

                    % add pearson correlation to plot
                    [pearson_r,pearson_p] = corr(x,y);

                    % set title
                    p.set_title(['Pearson r=' num2str(pearson_r)]);


                elseif iscategorical(x) && iscategorical(y)
                    % cat vs cat -> table plot
                    warning(['Plotting categorical vs categorical values not supported. ' ...
                        'Look at crosstable in console (undefined not shown).'])
                    disp(['crosstable: ' var_names{row_ind(r)} ' vs ' var_names{col_ind(c)}])
                    cross_tbl = crosstab(y,x);
                    cross_tbl = array2table(cross_tbl);
                    cross_tbl.Properties.VariableNames = categories(x);
                    cross_tbl.Properties.RowNames = categories(y);
                    disp(cross_tbl);
                    p = gramm();
                else 
                    % cont vs cat -> boxplot/volin plot

                    % flip
                    flip_data = iscategorical(y);
                    if flip_data
                        tmp = x;
                        x = y;
                        y = tmp;
                    end

                    % plot with categorical on x axis
                    p = gramm('x',x,'y',y,'color',color);
                    if strcmp(opts.catcont,'boxplot')
                        p.stat_boxplot();
                    elseif strcmp(opts.catcont,'jitter')
                        p.geom_jitter();
                    elseif strcmp(opts.catcont,'violin')
                        p.stat_violin();
                    else
                        error('Invalid catcont option')
                    end


                    % flip
                    if flip_data
                        p.coord_flip();
                    end
                end

                % set labels
                if c == 1
                    ystr = var_names{row_ind(r)};
                else
                    ystr = '';
                end
                if r == nrows
                    xstr = var_names{col_ind(c)};
                else
                    xstr = '';
                end
                p.set_names('x',xstr,'y',ystr);
                
                % assign to final plot grid
                g(r,c) = p;

            end
        end


        % draw final plot
        g.draw();
    end
    
    
    %---internal functions----
    function col_ind = find_col_ind(col_names,cell_specs)
        % parses columns specified by cell_specs into the indices of the
        % corresponding names in col_names
        
        % input: columns are specified by one of the following: column index,
        % column name as char or string, string array, numeric array, or cell array of
        % combinations of the first two.
        
        % if single name in chars
        if ischar(cell_specs)
            cell_specs = {cell_specs};
        end
        
        % if numeric array
        if isnumeric(cell_specs)
            invalid = find(cell_specs < 1 | numel(col_names) < cell_specs,1);
            if ~isempty(invalid)
                error(['Invalid input: ' num2str(cell_specs(invalid))]);
            end
            col_ind = cell_specs;
            return;
        end
        
        % else: string or cell array
        n = numel(cell_specs);
        col_ind = nan(n,1);
        for i = 1:n
            val = cell_specs{i};
            if isnumeric(val)
                if val < 1 || numel(col_names) < val
                    error(['Invalid input: ' string(val)]);
                end
                col_ind(i) = val;
            elseif isstring(val) || ischar(val)
                pos = find(strcmp(col_names,val),1);
                if isempty(pos)
                    error(['Invalid input: ' num2str(val)]);
                end
                col_ind(i) = pos;
            else
                error(['Invalid input: ' string(val)])
            end
        end
    end

end