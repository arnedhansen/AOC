% Define letterSequence depending on block iteration
function letterSequence = createLetterSequence(TRAINING, COND)
    if TRAINING == 1
        if COND == 1
            letterSequence = 'METTHLLAABBUZHH';
        elseif COND == 2
            letterSequence = 'SKLKNUNNTRTSPSS';
        elseif COND == 3
            letterSequence = 'SKLSKLNUNNTRTSP';
        end
    else
        letterSequence = generateLetterSequence(COND);
    end
end

function letterSequence = generateLetterSequence(COND)
    % Create consonant sequence
    consonants = char(setdiff('A':'Z', 'AEIOUW')); % Exclude vowels and 'W'
    consonants102 = [consonants consonants consonants consonants(1:end-2)];
    randConsonants = consonants102(randperm(length(consonants102)));

    % Ensure around 33% matching pairs based on COND
    letterSequence = createMatchingPairs(randConsonants, COND);

    % Check pseudorandom match probability and letter grouping
    while ~checkPRMP(letterSequence, COND) || ~checkLetterGrouping(letterSequence)
        randConsonants = consonants102(randperm(length(consonants102)));
        letterSequence = createMatchingPairs(randConsonants, COND);
    end
end

function letterSequence = createMatchingPairs(randConsonants, COND)
    letterSequence = randConsonants;
    for pairs = COND+1:length(letterSequence)
        if rand() <= 0.33
            letterSequence(pairs) = letterSequence(pairs-COND);
        end
        if sum(diff(double(letterSequence)) == 0) / length(letterSequence) * 100 >= 33
            break;
        end
    end
end

function isValid = checkPRMP(letterSequence, COND)
    matchCount = 0;
    for idx = COND+1:length(letterSequence)
        if letterSequence(idx) == letterSequence(idx-COND)
            matchCount = matchCount + 1;
        end
    end
    pseudoRandomMatchProbability = matchCount / length(letterSequence) * 100;
    isValid = pseudoRandomMatchProbability >= 32 && pseudoRandomMatchProbability <= 34;
    disp(['Check for pseudorandom match probability: ' num2str(pseudoRandomMatchProbability) ' %.']);
end

function isValid = checkLetterGrouping(letterSequence)
    grouping = 0;
    for i = 4:length(letterSequence)-3
        if all(letterSequence(i) == letterSequence(i-1:i+2)) || all(letterSequence(i) == letterSequence(i-1:-1:i-3)) || all(letterSequence(i) == letterSequence(i+1:i+3))
            grouping = grouping + 1;
        end
    end
    isValid = grouping == 0;
    if ~isValid
        disp(['Grouping of letters in sequence: ' num2str(grouping) '. Creating new sequence...']);
    else
        disp('No grouping of letters > 3. Continuing...');
    end
end
