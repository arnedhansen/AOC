addpath('/Users/tpopov/Documents/matlabtools/eeglab2020_0');
eeglab; close all hidden
clear all; close all;

% Subject IDs
subjects = {'301','302','304','306','309','310','312','313','314','315', ...
    '319','320','321','325','327','328','329','330','331','335', ...
    '336','339','341','342','343','344','345','346','347','348', ...
    '349','352','355','356','358','359','361','362','364','365', ...
    '368','372','373','374','375','377','378','379','386','388', ...
    '389','390','391','394','397','398','401','403','404','406', ...
    '412','413','416','419'};

base_dir = '/Volumes/TOURO/arne/merged';
%%
for s = 57:length(subjects)
    subj = subjects{s};
    subj_dir = fullfile(base_dir, subj);
    fprintf('\n--- Processing subject %s in %s ---\n', subj, subj_dir);
    
    if ~exist(subj_dir, 'dir')
        warning('Subject directory %s not found. Skipping...', subj_dir);
        
    end
    
    cd(subj_dir);
    
    %% ----- Sternberg Task ----- %%
    block_data = struct('load2', [], 'load4', [], 'load6', []);
    for b = 1:10
        filename = sprintf('%s_EEG_ET_Sternberg_block%d_merged.mat', subj, b);
        if exist(filename, 'file')
            tmp = load(filename); EEG = tmp.EEG;
            
            if b <= 2
                trigger = '52';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                block_data.load2{end+1} = data;
            elseif b <= 4
                trigger = '54';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                block_data.load4{end+1} = data;
            elseif b <= 6
                trigger = '56';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                block_data.load6{end+1} = data;
            end
        end
    end
    
    %     if ~isempty(block_data.load2)
    %         data2 = ft_appenddata([], block_data.load2{:});
    %         save(fullfile(subj_dir, [subj '_Sternberg_load2.mat']), 'data2');
    %     end
    %     if ~isempty(block_data.load4)
    %         data4 = ft_appenddata([], block_data.load4{:});
    %         save(fullfile(subj_dir, [subj '_Sternberg_load4.mat']), 'data4');
    %     end
    %     if ~isempty(block_data.load6)
    %         data6 = ft_appenddata([], block_data.load6{:});
    %         save(fullfile(subj_dir, [subj '_Sternberg_load6.mat']), 'data6');
    %     end
    % Append all Sternberg data regardless of load
    allSternberg = [block_data.load2, block_data.load4, block_data.load6];
    if ~isempty(allSternberg)
        dataSternberg = ft_appenddata([], allSternberg{:});
        save(fullfile(subj_dir, [subj '_Sternberg_all.mat']), 'dataSternberg');
    end
    %% ----- N-back Task ----- %%
    nback_data = struct('oneback', [], 'twoback', [], 'threeback', []);
    for b = 1:10
        filename = sprintf('%s_EEG_ET_Nback_block%d_merged.mat', subj, b);
        if exist(filename, 'file')
            tmp = load(filename); EEG = tmp.EEG;
            eventtypes = unique({EEG.event.type});
            
            if any(strcmp(eventtypes, '21'))
                trigger = '21';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                nback_data.oneback{end+1} = data;
            elseif any(strcmp(eventtypes, '22'))
                trigger = '22';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                nback_data.twoback{end+1} = data;
            elseif any(strcmp(eventtypes, '23'))
                trigger = '23';
                EEGload = pop_epoch(EEG, {trigger}, [-2 4]);
                data = eeglab2fieldtrip(EEGload, 'raw');
                data.trialinfo = repmat(str2double(trigger), length(data.trial), 1);
                nback_data.threeback{end+1} = data;
            else
                warning('No known N-back triggers (21/22/23) in %s', filename);
            end
        end
    end
    
    %     if ~isempty(nback_data.oneback)
    %         data1b = ft_appenddata([], nback_data.oneback{:});
    %         save(fullfile(subj_dir, [subj '_Nback_oneback.mat']), 'data1b');
    %     end
    %     if ~isempty(nback_data.twoback)
    %         data2b = ft_appenddata([], nback_data.twoback{:});
    %         save(fullfile(subj_dir, [subj '_Nback_twoback.mat']), 'data2b');
    %     end
    %     if ~isempty(nback_data.threeback)
    %         data3b = ft_appenddata([], nback_data.threeback{:});
    %         save(fullfile(subj_dir, [subj '_Nback_threeback.mat']), 'data3b');
    %     end
    allNback = [nback_data.oneback, nback_data.twoback, nback_data.threeback];
    if ~isempty(allNback)
      
        dataNback = ft_appenddata([], allNback{:});
        save(fullfile(subj_dir, [subj '_Nback_all.mat']), 'dataNback');
    end
    
end
