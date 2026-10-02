function ber = analyticalBER(sys, snrdB, AN)
% ANALYTICALBER  Union bound on the BER of joint-ML MIMO-OFDM-IM.
%   Efficient evaluation:
%   1) All single-antenna pairs (x, xhat) are grouped into classes of equal
%      difference vector e (up to a common phase, which does not change the
%      PEP). Class u has c(u) pairs and total Hamming weight s(u).
%   2) Joint errors are Nt-tuples of classes; tuple weight is
%      w = sum_t s(u_t) * prod_{v~=t} c(u_v)  (all-zero tuple has w = 0).
%   3) det(I + q*Mx) = sum_k e_k(Mx) q^k : the characteristic coefficients
%      e_k are computed once per tuple (Newton identities), then the PEP is
%      evaluated for all SNR values at once.
    N = sys.N; Q = sys.Q; Nt = sys.Nt; Nr = sys.Nr;

    % ---- 1) single-antenna difference classes
    [I, J] = ndgrid(1:Q, 1:Q);
    E  = sys.C(:, I(:)) - sys.C(:, J(:));               % N x Q^2
    dH = sys.D(:).';                                    % Hamming distances
    nz = abs(E) > 1e-9;
    [hasNz, first] = max(nz, [], 1);                    % first non-zero entry
    ref = E(sub2ind(size(E), first, 1:size(E,2)));
    ph  = ones(1, size(E,2));
    ph(hasNz) = conj(ref(hasNz)) ./ abs(ref(hasNz));    % remove common phase
    E   = E .* ph;
    key = round([real(E); imag(E)].' * 1e8) / 1e8;
    [~, iu, cls] = unique(key, 'rows');
    Eu  = E(:, iu);                                     % N x U representatives
    U   = numel(iu);
    cnt = accumarray(cls(:), 1).';                      % c(u)
    sw  = accumarray(cls(:), dH(:)).';                  % s(u)

    % ---- SNR dependent constants
    N0  = noisePower(sys, snrdB);
    pw  = (0:N).';
    Pq1 = (1 ./ (4*N0)).^pw;                            % (N+1) x nSNR
    Pq2 = (1 ./ (3*N0)).^pw;

    % ---- 2)+3) sum over Nt-tuples of classes, in memory-friendly chunks
    nTot    = U^Nt;
    sampled = nTot > AN.maxTerms;
    nEval   = min(nTot, AN.maxTerms);
    chunk   = max(1, floor(4e6 / N^3));
    acc     = zeros(1, numel(snrdB));
    fprintf('  union bound: %d difference classes, %g tuples%s\n', U, nTot, ...
            repmat(' (sampled)', 1, sampled));

    for i0 = 1:chunk:nEval
        nb = min(chunk, nEval - i0 + 1);
        if sampled, lin = randi(nTot, nb, 1); else, lin = (i0:i0+nb-1).'; end
        uT = mod(floor((lin - 1) ./ U.^(0:Nt-1)), U) + 1;   % nb x Nt classes

        % weight of every tuple
        cT = reshape(cnt(uT), size(uT));                % class sizes
        sT = reshape(sw(uT),  size(uT));                % class bit weights
        w  = zeros(nb, 1);
        for t = 1:Nt
            w = w + sT(:, t) .* prod(cT(:, [1:t-1, t+1:Nt]), 2);
        end
        keep = w > 0;                                   % drop zero-weight terms
        if ~any(keep), continue; end
        uT = uT(keep, :); w = w(keep); nb = numel(w);

        % Mx = Kc .* sum_t e_t e_t^H        (N x N x nb)
        Gm = zeros(N, N, nb);
        for t = 1:Nt
            e  = reshape(Eu(:, uT(:, t)), N, 1, nb);
            Gm = Gm + e .* conj(reshape(e, 1, N, nb));
        end
        Mx = sys.Kc .* Gm;

        % characteristic coefficients via Newton identities
        ek = charCoeffs(Mx);                            % (N+1) x nb

        % PEP for every tuple and SNR, weighted sum
        PEP = (1/12) * (ek.' * Pq1).^(-Nr) + (1/4) * (ek.' * Pq2).^(-Nr);
        acc = acc + w.' * PEP;
    end
    if sampled, acc = acc * nTot / nEval; end
    ber = acc / (Nt * sys.p * Q^Nt);
end
