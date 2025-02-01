%% Randomization of Subject IDs for AOC
allnums = 301:420;
taken = [331, 368, 397, 348, 388, 320, 373, 404, 310, 352, 419, 325, 355, 366, 312, 336, 412, 369, 391, 390, 306, 386, 401];

% Remove the taken values
remaining = setdiff(allnums, taken);

% Randomize the order
randomized_list = remaining(randperm(length(remaining)));
size(unique(randomized_list), 2)

% Display the randomized list
disp(randomized_list');