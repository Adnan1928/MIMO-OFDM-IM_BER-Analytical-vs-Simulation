function N0 = noisePower(sys, snrdB)
% NOISEPOWER  Frequency-domain noise variance N0 for the chosen SNR axis.
    snr = 10.^(snrdB/10);
    if strcmpi(sys.snrMode, 'EbN0')
        Eb = sys.Eg * (sys.NF + sys.CP) / (sys.NF * sys.p);   % energy / bit
        N0 = Eb ./ snr;
    else
        N0 = (sys.Eg / sys.N) ./ snr;                         % energy / subcarrier
    end
end
