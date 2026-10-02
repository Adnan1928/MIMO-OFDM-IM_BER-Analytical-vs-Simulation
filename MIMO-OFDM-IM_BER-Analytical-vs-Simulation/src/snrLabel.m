function s = snrLabel(mode)
% SNRLABEL  x-axis label for the selected SNR definition.
    if strcmpi(mode, 'EbN0'), s = 'E_b/N_0 (dB)'; else, s = 'E_s/N_0 (dB)'; end
end
