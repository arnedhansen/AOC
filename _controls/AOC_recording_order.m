% Read the Excel file. (Optionally, use 'VariableNamingRule','preserve' to keep original headers)
T = readtable('/Volumes/methlab_vp/AOC/AOC_VPs.xlsx');

% Combine the 'Datum' (date) and 'Uhrzeit' (time as fraction of day) columns.
% Adding the numeric fraction to a datetime will correctly offset the time.
T.DateTime = T.Datum + T.Uhrzeit;

% Sort the table by the combined datetime
[~, sortOrder] = sort(T.DateTime);

% Display the order of measurements (row indices) in the command window.
disp('Measurement order (row indices):');
disp(sortOrder);