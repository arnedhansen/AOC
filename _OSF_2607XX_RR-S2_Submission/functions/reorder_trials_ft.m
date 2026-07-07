function data = reorder_trials_ft(data, order)
% Reorder FieldTrip raw-like data by an index vector "order".
% Handles trial, time, trialinfo, sampleinfo, and cfg.trl if present.

    % basic containers
    if isfield(data, 'trial'),     data.trial     = data.trial(order);        end
    if isfield(data, 'time'),      data.time      = data.time(order);         end
    if isfield(data, 'trialinfo'), data.trialinfo = data.trialinfo(order,:);  end
    if isfield(data, 'sampleinfo'),data.sampleinfo= data.sampleinfo(order,:); end

    % cfg.trl appears in many pipelines; keep it aligned if present
    if isfield(data, 'cfg') && isstruct(data.cfg) && isfield(data.cfg, 'trl')
        try
            data.cfg.trl = data.cfg.trl(order, :);
        catch
            % fail-soft if shape differs
        end
    end
end
