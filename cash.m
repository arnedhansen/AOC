%% Cash calculation for OCC_NBack.m and OCC_Sternberg.m

%% Create dialog box to ask for subjectID

defAns = {'999'};
while true
    prompt = {'Subject Number'};
    box = inputdlg(prompt, 'Enter Subject Information', 1,defAns);
    subjectID = char(box(1));
    if length(subjectID) == 3               % Ensure response made in subject ID
        break
    end
end

% Convert subject data to numeric
subjectID = str2num(subjectID);

%% Load data file and extract cash amount
try
load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block6_1back_task.mat']);
load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block6_2back_task.mat']);
load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_NBack_block6_3back_task.mat']);
end
cashNback = saves.amountCHFextraTotal;

load([DATA_PATH, '/', num2str(subjectID), '/', num2str(subjectID), '_OCC_Sternberg_block8_task.mat']);
cashSternberg = saves.amountCHFextraTotal;

cashTotal = cashNback + cashSternberg;
cashTotal = round(cashTotal, 2);

%% Display Cash for Participant

msgbox(['Participant OCC', num2str(subjectID), ' has earned CHF ', num2str(cashTotal), ' in total.'])