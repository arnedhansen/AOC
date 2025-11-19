%% AOC Split of Nback EEG Data
% Split Nback participants by condition (1-back, 2-back, 3-back) order.
% Check if gaze metrics differ between these groups.

%% Load data
clc
clear
close all
data = readtable('/Volumes/methlab/Students/Arne/AOC/data/features/merged_data_sternberg.csv');
data.ReactionTime = data.ReactionTime .* 1000;
data.Condition = data.Condition .* 0.5;
data.GazeStd = (data.GazeStdX + data.GazeStdY) ./ 2;
variables = {'Accuracy', 'ReactionTime', 'GazeDeviation', 'GazeStd', 'MSRate', 'Fixations', 'Saccades', 'AlphaPower', 'IAF'};
save_names = {'acc', 'rt', 'gazedev', 'ms', 'blink', 'fix', 'sacc', 'pow', 'iaf'};
colors = color_def('AOC');
subjects = unique(data.ID);

%% Find Condition Order
path = '/Volumes/methlab/Students/Arne/AOC/data/merged/';
for subj = 1:length(subjects)
    subjID = subjects(subj);
    datapath = strcat(path, num2str(subjID));
    cd(datapath)
    files = dir(['/Volumes/methlab_data/OCC/AOC/data/', num2str(subjID), filesep, num2str(subjID), '_AOC_Nback_block*_task.mat']);
    for i = 1:numel(files)
        fname = files(i).name;
        tok = regexp( fname, ...
            '_NBack_block(\d+)_(\d)back_task\.mat$', ...
            'tokens', 'once' );
        blockIdx = str2double( tok{1} );
        backLev  = str2double( tok{2} );
        condOrder(blockIdx) = backLev;
    end
    disp(condOrder)

    % Save in data
    data.CondOrder(data.ID == subjects(subj)) = [{condOrder}; {condOrder}; {condOrder}];
end
disp('Finding  Condition Orders DONE')

%%
percent_struct = struct();
for s = 1:numel(subjects)
    sid = subjects(s);
    T = data(data.ID==sid, :);
    percent_struct(s).ID = sid;
    
    % get starting condition (first element of the CondOrder vector)
    ord = T.CondOrder{1};
    percent_struct(s).StartCond = ord(1);
    
    if height(T)==3
        for v = 1:numel(variables)
            low_val  = T{T.Condition==1, variables{v}};
            high_val = T{T.Condition==3, variables{v}};
            if ~isempty(low_val) && ~isempty(high_val)
                pchg = (high_val - low_val) / low_val * 100;
            else
                pchg = NaN;
            end
            percent_struct(s).(save_names{v}) = pchg;
        end
    else
        % missing data → NaNs
        for v = 1:numel(variables)
            percent_struct(s).(save_names{v}) = NaN;
        end
    end
end

%% PLOT %-CHANGE BARPLOTS, COLOURED BY STARTING CONDITION
figure('Position',[100 100 2000 1200],'Color','w');
conds = unique([percent_struct.StartCond]);
cols = lines(numel(conds));  % distinct colours

for v = 1:numel(variables)
    subplot(3,3,v); hold on;
    values = [percent_struct.(save_names{v})];
    starts = [percent_struct.StartCond];
    
    for s = 1:numel(subjects)
        cidx = find(conds==starts(s));
        bar(s, values(s), 'FaceColor', cols(cidx,:), 'EdgeColor','none');
    end
    
    xlim([0.5 numel(subjects)+0.5]);
    abw = max(abs(values),[],'omitnan');
    ylim([-abw*1.25 abw*1.25]);
    if abw>100, ylim([-100 100]); end
    
    xticks(1:numel(subjects));
    xticklabels(subjects);
    xlabel('Subject');
    ylabel('% Change');
    title(variables{v}, 'FontSize',16);
    hold off;
end

% add legend for starting conditions
%legend(arrayfun(@(c) sprintf('Start Cond %d',c),conds,'uni',false),'Location','bestoutside');
legend()

saveas(gcf, fullfile('/Volumes/methlab/Students/Arne/AOC/figures/tests',...
    'AOC_split_startcond.png'));

%% GLMMs: TEST EFFECT OF STARTING CONDITION ON METRIC CHANGE
% We’ll loop through each metric and fit:
%   pct_change ~ StartCond + (1|ID)
clc
results = struct();

for v = 8%1:numel(variables)
    sn = save_names{v};
    
    % build a table for modelling
    tbl = table();
    tbl.ID        = [percent_struct.ID]';
    tbl.StartCond = categorical([percent_struct.StartCond]');
    tbl.PctChange = [percent_struct.(sn)]';
    
    % remove NaNs
    tbl = tbl(~isnan(tbl.PctChange),:);
    
    % formula
    glme = fitglme(tbl, ...
        'PctChange ~ StartCond + (1|ID)', ...
        'Distribution','Normal', 'Link','identity');
    
    results.(sn) = glme;
    
    % display summary
    fprintf('\n===== GLMM for %s =====\n', variables{v});
    disp(glme);
end
