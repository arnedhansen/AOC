%% compute fooof for repetitions
tic
clear ff allff tmpdat
% iterate every 50 ms and compute fooof

VAR = ftdata.trialinfo(:, 1); % REPETITIONS

% get all combinations facilitate ft_analysis within parfor
var = unique(VAR);
c1 = 1;
all_combs = [];
for r1 = 1 : length(var)
    numStim = length(find(ismember(VAR, var(r1))));
    if numStim >= 5
        all_combs(c1, 1) = var(r1);
        all_combs(c1, 2) = acc(r1);
        all_combs(c1, 3) = lr(r1);
        all_combs(c1, 4) = numStim;
        c1 = c1 + 1;
    end
end

if not(isempty(all_combs))
    
    tot = round(TOI/STEP) + 1; 
    
    for t = 1 : tot 

        % select data
        cfg = [];
        cfg.latency = startlat+0.05*t;
        tmpdat = ft_selectdata(cfg,ftdata);

        clear tmpfreq tmptrialinfo tmpfooof
        % loop over trials
        parfor r = 1 : size(all_combs, 1)

            cfg = [];
            cfg.method     = 'mtmfft';
            cfg.taper      = 'hanning';
            cfg.foilim     = FOI;
            cfg.pad        = 5;
            cfg.output     = 'fooof'; 
            cfg.trials     = find(ismember(VAR, all_combs(r, 1)));

            % fooof
            tmpfreq{r} = ft_freqanalysis_methlab(cfg, tmpdat);
            tmpfooof{r, 1}  =  tmpfreq{r}.fooofparams; % save fooof params
            
            % test
            % mean(cell2mat({tmpfreq{r}.fooofparams.r_squared}))

        end

        % compute avg over trials
        cfg = [];
        cfg.keepindividual = 'yes';
        ff = ft_freqgrandaverage(cfg, tmpfreq{:});
        ff.fooofparams = tmpfooof;
        % ff.trialinfo = tmptrialinfo;
        ff.cfg=[];
        allff{t}=ff;
    end
    toc

    % concat data over time points
    tfr_ff = tmpfreq{1};
    fooofed_power = [];
    ff_foof = {};
    for t = 1 : length(allff)
        fooofed_power(:, :, :, t) = allff{t}.powspctrm;
        ff_foof{t} = allff{t}.fooofparams;
    end

    % extract aperiodic signal
    tmp_aperiodic = [];
    power_spectrum = [];
    for t = 1 : length(ff_foof)

        tdata = ff_foof{t};

        for trl = 1 : length(tdata)

            repdata = tdata{trl};

            tmpaperdiodic = {repdata.aperiodic_params};
            tmperror = {repdata.error};
            tmpr_sq = {repdata.r_squared};
            tmp_pwr_spec = {repdata.power_spectrum};
            elec_data = [];
            datafit = [];
            pwr_spec = [];
            for e = 1 : size(tfr_ff.label, 2)
                elec_data(1, e, :) = tmpaperdiodic{e};
                datafit(1, e, 1) = tmperror{e};
                datafit(1, e, 2) = tmpr_sq{e};
                pwr_spec(e, :) = tmp_pwr_spec{e};
            end

            tmp_aperiodic(trl, :, 1, t) = elec_data(1, :, 1); % intercept
            tmp_aperiodic(trl, :, 2, t) = elec_data(1, :, 2); % slope
            tmp_aperiodic(trl, :, 3, t) = datafit(1, :, 1); % error
            tmp_aperiodic(trl, :, 4, t) = datafit(1, :, 2); % r squared
            power_spectrum(trl, :, :, t) = pwr_spec;
        end
    end


    tfr_ff.dimord = 'rpt_chan_freq_time';
    tfr_ff.powspctrm = fooofed_power;
    tfr_ff.power_spectrum = power_spectrum;
    tfr_ff.trialinfo = all_combs; % the same for all time points
    tfr_ff.old_trialinfo = ftdata.trialinfo;
    tfr_ff.fooofparams = tmp_aperiodic;
    tfr_ff.time = TIME; 


end




%% tests


%%