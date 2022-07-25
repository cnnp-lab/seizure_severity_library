function standard_title(Sz_preproc,opts)
% generates standardized title that contains all the important information 

    arguments
        Sz_preproc
        
        opts.freq = -1
        opts.str = ' '
    end
    s = 1;

    if opts.freq == [-1 -1]
        freqstr = '; ';
    else 
        freqstr = [', ' num2str(opts.freq(1)) '-' num2str(opts.freq(2)) ' Hz; '];
    end
    
    title(['s' Sz_preproc.patient_id{s} ', ' Sz_preproc.segment_id{s} ...
    ', ' convertStringsToChars(string(Sz_preproc.ilae_sz_type(s))) ...
    freqstr opts.str],...
    'Interpreter','none');

end