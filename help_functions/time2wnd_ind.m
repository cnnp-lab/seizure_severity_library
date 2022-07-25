function ind = time2wnd_ind(t_wndw,t,opts) 
% transforms time to window index for given window time points t_wndw, i.e
% returns largest index less or lessequals than t

arguments
    
   t_wndw
   t
   
   opts.op = 'lesseq'
end


    if strcmp(opts.op,'lesseq')
        % if else since maybe faster because for loop not initialized
        if length(t) == 1
            ind = sum(t_wndw <= t);
        else
            ind = zeros(length(t));
            for i = 1:length(t)
               ind(i) = sum(t_wndw <= t(i));
            end
        end
        
    elseif strcmp(opts.op,'less')
        
        if length(t) == 1
            ind = sum(t_wndw < t);
        else
            ind = zeros(length(t));
            for i = 1:length(t)
               ind(i) = sum(t_wndw < t(i));
            end
        end
    else
        error("Optional must be 'less' or 'lesseq'")
    end
        
    
    
    % idea: min(find(...
end


