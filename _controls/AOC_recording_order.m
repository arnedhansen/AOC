% Read the Excel file. (Optionally, use 'VariableNamingRule','preserve' to keep original headers)
T = readtable('/Volumes/methlab_vp/AOC/AOC_VPs.xlsx');

% Combine the 'Datum' (date) and 'Uhrzeit' (time as fraction of day) columns.
% Adding the numeric fraction to a datetime will correctly offset the time.
T.DateTime = T.Datum + T.Zeit;

% Sort the table by DateTime
[~, sortOrder] = sort(T.DateTime);

% Sort the entire table
T_sorted = T(sortOrder, :);

% Extract the sorted IDs
sorted_IDs = T_sorted.ID;
disp(sorted_IDs);

% Save
T_sorted = T_sorted(:, ["x_", "ID", "Vorname", "Nachname", "DateTime"]);
writetable(T_sorted, '/Volumes/methlab/Students/Arne/AOC/data/controls/dataAcquisitionOrder/AOC_sorted_table.xlsx');
