%% AOC OSF startup
% Adds local functions and available toolboxes to the MATLAB path.

function startup()
rootDir = fileparts(fileparts(mfilename('fullpath')));
functionsDir = fullfile(rootDir, 'functions');
toolboxesDir = fullfile(rootDir, 'toolboxes');

if isfolder(functionsDir)
    addpath(genpath(functionsDir));
end

if isfolder(toolboxesDir)
    addpath(genpath(toolboxesDir));
end
end
