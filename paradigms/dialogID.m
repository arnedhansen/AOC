% Create dialog box to record subject data
defAns = {'999','XX','99','X', 'X', 'X', '99'};
while true
    prompt = {'Subject Number', 'Subject Initials', 'Subject Age', 'Subject Gender', 'EEG Cap Size', 'Ocular Dominance', 'Time Hair Washing'};
    box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);
    subject.ID = str2double(box{1}); % Convert directly to a number
    subject.initials = char(box{2});
    subject.age = str2double(box{3}); % Convert directly to a number
    subject.gender = char(box{4});
    subject.capsize = char(box{5});
    subject.ocDom = char(box{6});
    subject.hairwash = str2double(box{7}); % Convert directly to a number
    % Ensure that Subject ID and initials are valid
    if ~isnan(subject.ID) && length(box{1}) == 3 && length(subject.initials) == 2
        break
    else
        errordlg('Invalid Subject ID or Initials. Please try again.');
    end
end

% Assign subjectID for further use
subjectID = subject.ID;