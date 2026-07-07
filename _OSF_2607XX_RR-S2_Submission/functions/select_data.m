%% Function to select data in FieldTrip structures
function selected_data = select_data(latency, frequency, data)
    cfg              = [];
    cfg.latency      = latency;
    cfg.frequency    = frequency;
    cfg.avgovertime  = 'yes';
    cfg.avgoverfreq  = 'no';
    cfg.avgoverchan  = 'no';
    cfg.nanmean      = 'yes';

    selected_data = ft_selectdata(cfg, data);

    % If time artifacts remain, prune safely
    if isfield(selected_data, 'time') || (isfield(selected_data, 'dimord') && contains(selected_data.dimord,'time'))
        selected_data = remove_time_dimension(selected_data);
    end
end
