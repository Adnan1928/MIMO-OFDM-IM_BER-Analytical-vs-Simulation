function ber = simulateBER(sys, snrdB, MC)
% SIMULATEBER  Monte-Carlo BER of MIMO-OFDM-IM with joint ML detection.
    N = sys.N; G = sys.G; NF = sys.NF; Nt = sys.Nt; Nr = sys.Nr; Q = sys.Q;
    N0v  = noisePower(sys, snrdB);
    tapA = sqrt(sys.pdp(:)/2);                       % tap std per real dim
    ber  = zeros(1, numel(snrdB));

    for s = 1:numel(snrdB)
        N0 = N0v(s); nErr = 0; nBit = 0; f = 0;
        while (nErr < MC.minErr || f < MC.minFrames) && f < MC.maxFrames
            f = f + 1;

            % --- Source + IM/symbol mapping: uniform random bit words are
            %     the same as uniform random codeword indices (G x Nt)
            idx = randi(Q, G, Nt);
            X   = reshape(sys.C(:, idx), N, G, Nt);         % N x G x Nt

            % --- Interleaving: subcarrier (n-1)*G + g  <-  X(n,g,t)
            Xf  = reshape(permute(X, [2 1 3]), NF, 1, Nt);  % NF x 1 x Nt

            % --- Rayleigh multipath channel (L taps) -> frequency response
            h   = (randn(sys.L, Nr, Nt) + 1j*randn(sys.L, Nr, Nt)) .* tapA;
            H   = fft(h, NF, 1);                            % NF x Nr x Nt

            % --- Received signal on every Rx antenna + AWGN (complex!)
            W   = sqrt(N0/2) * (randn(NF, Nr) + 1j*randn(NF, Nr));
            Y   = sum(H .* Xf, 3) + W;                      % NF x Nr

            % --- De-interleaving + grouping
            Yg  = permute(reshape(Y, G, N, Nr),     [2 1 3]);    % N x G x Nr
            Hg  = permute(reshape(H, G, N, Nr, Nt), [2 1 3 4]);  % N x G x Nr x Nt

            % --- Joint ML detection and bit-error counting
            idxHat = mlDetect(Yg, Hg, sys, MC.maxElems);         % G x Nt
            nErr   = nErr + sum(sys.D(sub2ind([Q Q], idx(:), idxHat(:))));
            nBit   = nBit + G * Nt * sys.p;
        end
        ber(s) = nErr / nBit;
        fprintf('  SNR = %5.1f dB | frames = %6d | errors = %6d | BER = %.3e\n', ...
                snrdB(s), f, nErr, ber(s));
        if nErr == 0, ber(s:end) = NaN; break; end     % below resolution
    end
end
