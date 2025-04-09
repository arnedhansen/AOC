function sendtrigger(trg,port,SITE,stayup)
if SITE == 'A' % A = ANT Neuro
    ppdev_mex('Write', port, trg);
    if ~stayup
        WaitSecs(0.010);
        ppdev_mex('Write', 1, 0);
    end
end


