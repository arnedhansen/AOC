%% Function to remove time dimension in FieldTrip data structure
function data_no_time = remove_time_dimension(data)
    data_no_time = data;
    % only remove if time is present or dimord mentions 'time'
    hasTimeField   = isfield(data_no_time, 'time');
    hasTimeInDimord = isfield(data_no_time, 'dimord') && contains(data_no_time.dimord, 'time');

    if hasTimeField || hasTimeInDimord
        if hasTimeField
            data_no_time = rmfield(data_no_time, 'time');
        end
        if isfield(data_no_time, 'dimord')
            data_no_time.dimord = strrep(data_no_time.dimord, '_time', '');
        end
    end

    % normalise to chan_freq if applicable
    if isfield(data_no_time, 'powspctrm')
        sz = size(data_no_time.powspctrm);
        if numel(sz) == 2
            data_no_time.dimord = 'chan_freq';
        end
    end
end
