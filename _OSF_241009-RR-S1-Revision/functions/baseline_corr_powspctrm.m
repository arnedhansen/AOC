%% Baseline correction for power spectra for TFR

function freqGA_data = baseline_corr_powspctrm(freqGA_data, baseline_period)
    % Find baseline period in data from ft_freqgrandaverage
    baseline_idx = find(freqGA_data.time >= baseline_period(1) & freqGA_data.time <= baseline_period(2));
    baseline_power = mean(freqGA_data.powspctrm(:,:,baseline_idx), 3, 'omitnan');
    freqGA_data.powspctrm = ((freqGA_data.powspctrm - baseline_power) ./ baseline_power) * 100;
end