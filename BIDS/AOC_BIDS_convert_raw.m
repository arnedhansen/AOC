function AOC_BIDS_convert_raw(cfg)
%AOC_BIDS_CONVERT_RAW  Convert AOC task-segmented raw files to BIDS layout.
%
% Recommended source: task-segmented .mat files in each subject folder
% (e.g. 301_AOC_Nback_block1_task_EEG.mat), NOT the continuous .cnt in archive/.
%
% Rationale:
%   - archive/*.cnt is the continuous ANT recording and would require users to
%     re-run the cutting pipeline to obtain task segments.
%   - archive/*.edf are EyeLink exports, not EEG.
%   - Task .mat EEG files are already segmented, include channel locations, and
%     match the analysis pipeline entry point after AOC_DataCuttingFunction.
%
% Excluded by design: training, resting, photodiode.
%
% Usage:
%   startup
%   setup('AOC')
%   addpath('/Users/Arne/Documents/GitHub/AOC/BIDS')
%   AOC_BIDS_convert_raw()
%   AOC_BIDS_convert_raw(struct('subjects', {'301'}))
%
% Outputs per subject:
%   sub-XXXX/eeg/sub-XXXX_task-<task>_run-XX_eeg.set (+ .fdt if needed)
%   sub-XXXX/eeg/sub-XXXX_task-<task>_run-XX_eeg.json
%   sub-XXXX/eeg/sub-XXXX_task-<task>_run-XX_channels.tsv
%   sub-XXXX/eeg/sub-XXXX_task-<task>_run-XX_events.tsv
%   sub-XXXX/eyetrack/sub-XXXX_task-<task>_run-XX_eyetrack.mat
%   sub-XXXX/eyetrack/sub-XXXX_task-<task>_run-XX_eyetrack.json
%   sub-XXXX/beh/sub-XXXX_task-<task>_run-XX_beh.mat
%   sub-XXXX/beh/sub-XXXX_task-<task>_run-XX_beh.json

if nargin < 1
    cfg = struct();
end

defaults = struct( ...
    'raw_root', '/Volumes/g_psyplafor_methlab_data$/OCC/AOC/data', ...
    'bids_root', '/Volumes/g_psyplafor_methlab_data$/OCC/AOC_BIDS', ...
    'subjects', {{}}, ...
    'overwrite', false, ...
    'use_eeglab', true);

cfg = merge_struct(defaults, cfg);

if isempty(cfg.subjects)
    d = dir(fullfile(cfg.raw_root, '[0-9]*'));
    d = d([d.isdir]);
    cfg.subjects = sort({d.name});
end

if cfg.use_eeglab
    if exist('eeglab', 'file') ~= 2
        error('EEGLAB not found on path. Run startup/setup or set use_eeglab to false.');
    end
    evalc('eeglab nogui');
end

fprintf('Converting %d subjects from %s\n', numel(cfg.subjects), cfg.raw_root);

for i = 1:numel(cfg.subjects)
    sub_id = cfg.subjects{i};
    fprintf('Subject %s (%d/%d)\n', sub_id, i, numel(cfg.subjects));
    convert_subject(sub_id, cfg);
end

fprintf('Raw conversion complete.\n');
end

function convert_subject(sub_id, cfg)
sub_dir = fullfile(cfg.raw_root, sub_id);
bids_sub = sprintf('sub-%04d', str2double(sub_id));

eeg_dir = fullfile(cfg.bids_root, bids_sub, 'eeg');
et_dir = fullfile(cfg.bids_root, bids_sub, 'eyetrack');
beh_dir = fullfile(cfg.bids_root, bids_sub, 'beh');
mkdir(eeg_dir);
mkdir(et_dir);
mkdir(beh_dir);

patterns = { ...
    struct('regex', '^(?<id>\d+)_AOC_Nback_block(?<run>\d+)_task_EEG\.mat$', 'task', 'nback'), ...
    struct('regex', '^(?<id>\d+)_AOC_Sternberg_block(?<run>\d+)_task_EEG\.mat$', 'task', 'sternberg'), ...
    struct('regex', '^(?<id>\d+)_AOC_Nback_block(?<run>\d+)_task_ET\.mat$', 'task', 'nback', 'modality', 'eyetrack'), ...
    struct('regex', '^(?<id>\d+)_AOC_Sternberg_block(?<run>\d+)_task_ET\.mat$', 'task', 'sternberg', 'modality', 'eyetrack'), ...
    struct('regex', '^(?<id>\d+)_AOC_NBack_block(?<run>\d+)_(?<load>\d+back)_task\.mat$', 'task', 'nback', 'modality', 'beh'), ...
    struct('regex', '^(?<id>\d+)_AOC_Sternberg_block(?<run>\d+)_task\.mat$', 'task', 'sternberg', 'modality', 'beh') ...
};

files = dir(sub_dir);
for f = 1:numel(files)
    if files(f).isdir
        continue;
    end
    name = files(f).name;
    if should_skip_file(name)
        continue;
    end

    matched = match_pattern(name, patterns);
    if isempty(matched)
        continue;
    end

    src = fullfile(sub_dir, name);
    run_num = str2double(matched.run);
    if isnan(run_num)
        tok = regexp(name, 'block(\d+)', 'tokens', 'once');
        run_num = str2double(tok{1});
    end
    run_label = sprintf('%02d', run_num);
    base = sprintf('%s_task-%s_run-%s', bids_sub, matched.task, run_label);

    switch matched.modality
        case 'eyetrack'
            dst = fullfile(et_dir, [base, '_eyetrack.mat']);
            if ~exist(dst, 'file') || cfg.overwrite
                copyfile(src, dst);
            end
            write_eyetrack_json(fullfile(et_dir, [base, '_eyetrack.json']), matched.task, run_num, src);

        case 'beh'
            dst = fullfile(beh_dir, [base, '_beh.mat']);
            if ~exist(dst, 'file') || cfg.overwrite
                copyfile(src, dst);
            end
            meta = struct('TaskName', matched.task, 'Run', run_num, 'SourceFile', name);
            if isfield(matched, 'load')
                meta.MemoryLoad = matched.load;
            end
            write_json(fullfile(beh_dir, [base, '_beh.json']), meta);

        otherwise
            convert_eeg_mat(src, eeg_dir, base, matched.task, run_num, cfg);
    end
end
end

function convert_eeg_mat(src, eeg_dir, base, task, run_num, cfg)
set_path = fullfile(eeg_dir, [base, '_eeg.set']);
json_path = fullfile(eeg_dir, [base, '_eeg.json']);
chan_path = fullfile(eeg_dir, [base, '_channels.tsv']);
evt_path = fullfile(eeg_dir, [base, '_events.tsv']);

if exist(set_path, 'file') && ~cfg.overwrite
    return;
end

S = load(src);
if isfield(S, 'EEG')
    EEG = S.EEG;
else
    error('No EEG struct in %s', src);
end

if cfg.use_eeglab
    pop_saveset(EEG, 'filename', [base, '_eeg.set'], 'filepath', eeg_dir);
else
    EEG = EEG;
    save(set_path, 'EEG', '-v7.3');
end

write_eeg_json(json_path, EEG, task, run_num, src);
write_channels_tsv(chan_path, EEG);
write_events_tsv(evt_path, EEG);
end

function tf = should_skip_file(name)
skip_tokens = {'Training', 'training', 'block0', 'Resting', 'Photodiode', 'Photodiode', 'Res.edf', 'NTr', 'STr'};
tf = any(contains(name, skip_tokens));
end

function matched = match_pattern(name, patterns)
matched = [];
for i = 1:numel(patterns)
    tok = regexp(name, patterns{i}.regex, 'names');
    if ~isempty(tok)
        matched = patterns{i};
        matched.run = tok.run;
        if isfield(tok, 'load')
            matched.load = tok.load;
        end
        return;
    end
end
end

function write_eeg_json(json_path, EEG, task, run_num, src)
meta = struct();
meta.TaskName = task;
meta.Run = run_num;
meta.RecordingType = 'continuous';
meta.SamplingFrequency = EEG.srate;
meta.PowerLineFrequency = 50;
meta.EEGReference = 'CPz';
meta.EEGGround = 'AFz';
meta.EEGChannelCount = EEG.nbchan;
meta.SourceFormat = 'MATLAB/EEGLAB struct from AOC_DataCuttingFunction';
meta.SourceFile = src;
if isfield(EEG, 'times') && ~isempty(EEG.times)
    meta.Duration = (numel(EEG.times) / EEG.srate);
end
write_json(json_path, meta);
end

function write_eyetrack_json(json_path, task, run_num, src)
meta = struct();
meta.TaskName = task;
meta.Run = run_num;
meta.RecordingType = 'eyetracking';
meta.SamplingFrequency = 1000;
meta.SourceFormat = 'MATLAB/EyeLink parseeyelink output';
meta.SourceFile = src;
write_json(json_path, meta);
end

function write_channels_tsv(chan_path, EEG)
fid = fopen(chan_path, 'w');
fprintf(fid, 'name\ttype\tunits\n');
for i = 1:EEG.nbchan
    label = '';
    if i <= numel(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels')
        label = EEG.chanlocs(i).labels;
    end
    if isempty(label)
        label = sprintf('E%d', i);
    end
    ch_type = 'EEG';
    if strcmpi(label, 'CPz')
        ch_type = 'REF';
    end
    fprintf(fid, '%s\t%s\tuV\n', label, ch_type);
end
fclose(fid);
end

function write_events_tsv(evt_path, EEG)
fid = fopen(evt_path, 'w');
fprintf(fid, 'onset\tduration\ttrial_type\n');
if isfield(EEG, 'event') && ~isempty(EEG.event)
    for i = 1:numel(EEG.event)
        onset = NaN;
        if isfield(EEG.event, 'latency')
            lat = EEG.event(i).latency;
            if isfield(EEG, 'srate') && ~isempty(EEG.srate)
                onset = (lat - 1) / EEG.srate;
            else
                onset = lat;
            end
        end
        etype = '';
        if isfield(EEG.event, 'type')
            etype = EEG.event(i).type;
        end
        fprintf(fid, '%.6f\t0\t%s\n', onset, etype);
    end
end
fclose(fid);
end

function write_json(json_path, s)
txt = jsonencode(s);
fid = fopen(json_path, 'w');
fprintf(fid, '%s\n', txt);
fclose(fid);
end

function out = merge_struct(defaults, overrides)
out = defaults;
if isempty(overrides)
    return;
end
names = fieldnames(overrides);
for i = 1:numel(names)
    out.(names{i}) = overrides.(names{i});
end
end
