%% AOC setup
% Resolves project paths and optional toolbox initialization for OSF use.
% Usage:
%   [subjects, paths, colors, headmodel] = setup('AOC');
%   [subjects, paths, colors, headmodel] = setup('AOC', 0);

function [subjects, paths, colors, headmodel] = setup(projectName, initToolboxes)
if nargin < 1 || isempty(projectName)
    projectName = 'AOC';
end
if nargin < 2
    initToolboxes = 1;
end

if ~strcmpi(projectName, 'AOC')
    error('setup:UnsupportedProject', 'Only projectName ''AOC'' is supported.');
end

rootFromEnv = getenv('AOC_OSF_ROOT');
if ~isempty(rootFromEnv) && isfolder(rootFromEnv)
    rootDir = rootFromEnv;
else
    rootDir = fileparts(fileparts(mfilename('fullpath')));
end

paths = struct();
paths.root = rootDir;
paths.functions = fullfile(rootDir, 'functions');
paths.toolboxes = fullfile(rootDir, 'toolboxes');
paths.data = fullfile(rootDir, 'data');
paths.raw = fullfile(paths.data, 'raw');
paths.raw_occ = fullfile(paths.raw, 'occ');
paths.features = fullfile(paths.data, 'features');
paths.figures = fullfile(paths.data, 'figures');
paths.stats = fullfile(rootDir, 'stats');
paths.vp_table = fullfile(paths.raw_occ, 'AOC_VP_List.xlsx');
paths.seb_path = '';
paths.headmodel_path = fullfile(paths.toolboxes, 'fieldtrip', 'template', 'headmodel', 'standard_bem.mat');

if isfolder(paths.features)
    dirs = dir(paths.features);
    folders = dirs([dirs.isdir] & ~ismember({dirs.name}, {'.', '..'}));
    subjects = {folders.name};
else
    subjects = {};
end

colors = struct();
colors.nback = [0.20 0.45 0.80; 0.30 0.65 0.30; 0.85 0.35 0.25];
colors.sternberg = [0.20 0.45 0.80; 0.55 0.55 0.55; 0.85 0.35 0.25];

headmodel = [];
if exist(paths.headmodel_path, 'file') == 2
    S = load(paths.headmodel_path);
    if isfield(S, 'vol')
        headmodel = S.vol;
    end
end

if initToolboxes
    startup;
end
end
