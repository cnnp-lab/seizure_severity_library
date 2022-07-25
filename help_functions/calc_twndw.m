function t_wndw = calc_twndw(Sz,s,wndw_len, wndw_shift)
% calculates the window timepoints for the seizure s in Sz, 
    
    fs = Sz.segment_fs(s);
    pre = Sz.segment_pre(s);
    dur = Sz.duration(s);
    post = Sz.segment_post(s);
    
    wndw_shift_s = wndw_shift / fs;

    t_wndw = -pre:wndw_shift_s:(dur+post-(wndw_len-1)/fs);

end