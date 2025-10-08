clear all
close all

% Subject IDs
load('/Volumes/Homestore/OCC/arne/subjects.mat');

base_dir = '/Volumes/Homestore/OCC/arne/merged';

% Tasks and conditions
tasks = {'Sternberg', 'Nback'};
task_conditions = struct(...
    'Sternberg', [22, 24, 26], ...
    'Nback', [21, 22, 23]);

%% Loop over subjects
for s = 91%1:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    fprintf('\n--- Processing subject %s ---\n', subj);
    
    if ~exist(subj_dir, 'dir')
        fprintf('  ✘ Directory not found: %s\n', subj_dir);
        continue;
    end
    cd(subj_dir);
    
    for t = 1:length(tasks)
        task = tasks{t};
        conds = task_conditions.(task);
        
        datafile = fullfile(subj_dir, [subj '_' task '_all.mat']);
        if ~exist(datafile, 'file')
            fprintf('  > Skipping %s (no data found)\n', task);
            continue;
        end
        
        load(datafile);  % loads dataSternberg or dataNback
        data = eval(['data' task]);
        
        for c = 1:length(conds)
            cond = conds(c);
            
            % Trial selection
            cfg = [];
            cfg.trials = find(data.trialinfo == cond);
            if isempty(cfg.trials)
                fprintf('  > No trials for condition %d in task %s\n', cond, task);
                continue;
            end
            dataload = ft_selectdata(cfg, data);
            
            % Remove eye-tracking channels (assume last 6)
            eye_labels = dataload.label(end-5:end);  % channels 130–135
            cfg_rm = [];
            cfg_rm.channel = setdiff(dataload.label, eye_labels);
            dataload = ft_selectdata(cfg_rm, dataload);
            
            % Average reference
            cfg_ref = [];
            cfg_ref.reref = 'yes';
            cfg_ref.refchannel = 'all';
            dataload = ft_preprocessing(cfg_ref, dataload);
            
            %% Time-Frequency Analysis
            cfg = [];
            cfg.output     = 'pow';
            cfg.method     = 'mtmconvol';
            cfg.taper      = 'hanning';
            cfg.tapsmofrq  = 5;
            cfg.foi        = 1:1:40;
            cfg.t_ftimwin  = ones(length(cfg.foi), 1) * 0.5;
            cfg.toi        = -1.5:0.05:3;
            cfg.keeptrials = 'no';
            tfr = ft_freqanalysis(cfg, dataload);
            
            %% FOOOF (commented for now)
            clear fspctrm
            settings = struct();
            settings.peak_width_limits = [2 12];
            for t = 1:length(tfr.time)
                cfg = [];
                cfg.latency = tfr.time(t);
                tmp = ft_selectdata(cfg, tfr);
                
                for chan = 1:length(tmp.label)
                    freqs = tmp.freq';
                    psd = tmp.powspctrm(chan,:)';
                    fooof_results = fooof(freqs, psd, [min(freqs), max(freqs)], settings, true);
                    powspctrmff(chan,:) = fooof_results.fooofed_spectrum - fooof_results.ap_fit;
                end
                fspctrm(:,:,t) = powspctrmff;
            end
            
            tfr_fooof = tfr;
            tfr_fooof.powspctrm = fspctrm;
            save(fullfile(subj_dir, sprintf('%s_%s_cond%d_fooof.mat', subj, task, cond)), 'tfr_fooof');
         
            %% Save only the normal TFR (no FOOOF)
%             savefile = fullfile(subj_dir, sprintf('%s_%s_cond%d_tfr.mat', subj, task, cond));
%             save(savefile, 'tfr');
%             fprintf('  ✔ Saved: %s\n', savefile);
        end
    end
end

%% Example Plot
cfg = [];
cfg.figure = 'gcf';
cfg.linecolor = 'brk';
cfg.ylim = [3 40];

% Optional: pick a loaded file to visualize
% load('/Volumes/TOURO/arne/merged/301/301_Sternberg_cond52_tfr.mat');
% figure; ft_multiplotTFR(cfg, tfr);
