function vals = gridsearch(fct_hndl,par_cell,opts)

    % function that performs a gridsearch for the fct_hndl function for
    % arbitrary number of parameters and ranges. If plot handle given and
    % less or equal than 2 parameters, plots gridmatrix of the plots
    % generate by plot handle. If value handle, given amd less or equal
    % than 2 parameters, plots 2D or 3D plot of results from value handle
    % against parameters.
    
    % input:
    %   - fct_hndl: takes exactly p parameters as input. Only first ouput
    %       used!
    %   - par_cell: cell array with exactly p cells, where each cell
    %       contains parameter range. 
    %   - opts.plt_hndl: function that takes the output of fct_hndl and
    %       generates a plot into the current figure. The plot function
    %       needs to be compatible with matlab subplots (e.g. gramm doesnt 
    %       work here!). If p<=2, then algo
    %       plots these analysis plots for the parameter grid.
    %   - opts.val_hndl: function that takes the output of fct_hndl and
    %       converts it to a single numeric value. If p<=2, then algo plots 
    %       a 3D or 2D plot of results against parameters.
    %   - opts.fig_ind: 0 is new figure, value larger zero is that figure.
    %       If both plt_hndl val_hndl given, plt_hndl is plotted into
    %       fig_ind and val_hndl plot plotted into new figure.
    %   - opts.param_names: parameter names as cell or string array for
    %       plotting
    % output:
    %   - vals: cell of dimension p, that contains the
    %       first outputs of fct_hndl for all input parameter combinations.
    %       It holds: vals{a,b,c,...} =
    %       fct_hndl(parcell{1}(a),parcell{2}(b),...)

    arguments
        fct_hndl
        par_cell
        
        opts.param_names
        opts.plt_hndl
        opts.val_hndl
        opts.fig_ind = 0
    end
    
    % parse par_cell
    if ~iscell(par_cell)
        warning('parameter input interpreted as one parameter')
        par_cell = {par_cell};
    end
    npar = numel(par_cell);
    ninp = zeros(npar,1);
    for p = 1:npar
        ninp(p) = numel(par_cell{p});
    end
    
    % parse fig_ind
    val_flag = isfield(opts,'val_hndl') && npar <= 2;
    plt_flag = isfield(opts,'plt_hndl') && npar <= 2;
    if plt_flag
        if opts.fig_ind == 0
            figure()
        else
            figure(opts.fig_ind)
        end
    end
    
    % parse param_names
    if isfield(opts,'param_names')
        param_names = opts.param_names;
    else
        param_names = cell(npar,1);
        for p = 1:npar
            param_names{p} = ['par' num2str(p)];
        end
    end
    
    
    % calculate all parameter combinations
    comb_inp = cell(npar,1);
    for c = 1:npar
        comb_inp{c} = 1:ninp(c);
    end
    comb = allcomb(comb_inp{:});
        
    
    % main loop
    vals = cell(ninp');
    ncomb = size(comb,1);
    for c = 1:ncomb
        % specify arguments
        ind = comb(c,:);
        ind_cell = cell(1,npar);
        args = cell(npar,1); % as cell
        for a = 1:npar
            args{a} = par_cell{a}(ind(a));
            ind_cell{a} = ind(a);
        end
        
        % calc
        res = fct_hndl(args{:});
        
        % save calc
        vals{ind_cell{:}} = res;
        
        % plot the plot_hndl
        if plt_flag
            subplot(ninp(1),ninp(2),(ind(1)-1)*ninp(2)+ind(2))
            opts.plt_hndl(res);
            title([param_names{1} '=' num2str(args{1}) ', ' ...
                param_names{2} '=' num2str(args{2})],'Interpreter','None');
        end
    end
    
    
    % plot fct_hndl
    if val_flag
        
        % open plot
        if opts.fig_ind == 0
            figure();
        else
            figure(opts.fig_ind+1);
        end
        

        % actual plotting
        if npar == 1
            
             % extract params
            data_mat = nan(ninp(1),1);
            for p1 = 1:ninp(1)
                data_mat(p1) = opts.val_hndl(vals{p1});
            end
            
            % actual plotting
            plot(par_cell{1},data_mat)
            xlabel(param_names{1})
            ylabel('val_hndl value','Interpreter','None')
            
        else  % npar == 2
            
            % extract params
            data_mat = nan(ninp(1),ninp(2));
            for p1 = 1:ninp(1)
                for p2 = 1:ninp(2)
                    data_mat(p1,p2) = opts.val_hndl(vals{p1,p2});
                end
            end

            % actual plotting
%             imagesc(data_mat, 'YData',par_cell{1},'XData',par_cell{2})
            imagesc(data_mat)
            set(gca,'XTick',1:numel(par_cell{2}))
            set(gca,'YTick',1:numel(par_cell{1}))
            set(gca,'XTickLabels',round(par_cell{2},2))
            set(gca,'YTickLabels',par_cell{1})
            
            xlabel(param_names{2},'Interpreter','None')
            ylabel(param_names{1},'Interpreter','None')
            colorbar
        end
        
        title('Gridsearch')
    end
    
end