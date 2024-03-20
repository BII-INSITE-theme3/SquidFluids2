%% IGNORE THIS FILE - IT IS JUST USED FOR PROTOTYPING CODE

%%

% Define file paths

% Define directory path and pattern (e.g., all .m files)
dirPath = 'outputs/';
pattern = 'parameters_*.m';

% Get list of matching files
fileList = Glob(fullfile(dirPath, pattern));

% Display filenames
disp(length(fileList))
