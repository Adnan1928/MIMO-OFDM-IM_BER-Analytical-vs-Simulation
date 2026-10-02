function sys = buildSystem(S)
% BUILDSYSTEM  Derive all fixed quantities of one configuration:
%   codebook of one group (N x Q), bit labels, Hamming-distance table,
%   group energy and the channel correlation matrix inside one group.
    sys   = S;
    sys.N = S.NF / S.G;
    assert(mod(sys.N,1) == 0, 'NF must be a multiple of G.');
    assert(S.K >= 1 && S.K <= sys.N, 'K must satisfy 1 <= K <= N.');
    assert(S.CP >= S.L-1, 'CP must be at least L-1 (no ISI).');
    assert(abs(sum(S.pdp)-1) < 1e-9 && numel(S.pdp) == S.L, 'pdp must have L entries summing to 1.');
    N = sys.N; K = S.K; M = S.M;
    mb = log2(M);
    assert(mod(mb,1) == 0, 'M must be a power of 2.');

    % ---- index-modulation patterns: first 2^b of the C(N,K) combinations
    sys.b   = floor(log2(nchoosek(N, K)));
    allPat  = nchoosek(1:N, K);
    sys.pat = allPat(1:2^sys.b, :);                 % (2^b) x K

    % ---- constellations (Gray labelled, unit average energy)
    A = grayConstellation(M, S.modType);            % primary (active)
    if strcmpi(S.modType,'PSK'), rot = exp(1j*pi/M); else, rot = exp(1j*pi/4); end
    B = A * rot;                                    % secondary (DM mode)

    % ---- bits per group and number of codewords
    isDM   = strcmpi(S.scheme, 'DM-OFDM-IM');
    nSym   = K + isDM*(N-K);                        % symbols carrying bits
    sys.p  = sys.b + nSym*mb;
    sys.Q  = 2^sys.p;
    Q      = sys.Q;

    % ---- codebook: column q+1 is the group signal for bit word q (MSB first)
    q      = 0:Q-1;
    patIdx = floor(q / 2^(sys.p - sys.b)) + 1;      % index bits -> pattern
    symInt = mod(q, 2^(sys.p - sys.b));             % symbol bits as integer
    Cb     = zeros(N, Q);
    act    = sys.pat(patIdx, :).';                  % K x Q active positions
    for k = 1:K                                     % active subcarriers
        lab = mod(floor(symInt / M^(nSym-k)), M);
        Cb(sub2ind([N Q], act(k,:), q+1)) = A(lab+1);
    end
    if isDM                                         % idle subcarriers
        for qq = 1:Q
            idle = setdiff(1:N, act(:,qq));
            for k = 1:N-K
                lab = mod(floor(symInt(qq) / M^(nSym-K-k)), M);
                Cb(idle(k), qq) = B(lab+1);
            end
        end
    end
    sys.C  = Cb;
    sys.Eg = mean(sum(abs(Cb).^2, 1));              % avg energy / group / antenna

    % ---- Hamming distance between all bit words (Q x Q)
    D = zeros(Q);
    for k = 1:sys.p
        bk = bitget(q, k);
        D  = D + double(xor(bk.', bk));
    end
    sys.D = D;

    % ---- channel correlation of the N subcarriers of one group
    %      Kc(n1,n2) = sum_l pdp(l) exp(-j2*pi*(l-1)*(n1-n2)*G/NF)
    dn     = ((0:N-1).' - (0:N-1)) * S.G;
    sys.Kc = reshape(exp(-1j*2*pi*dn(:)*(0:S.L-1)/S.NF) * S.pdp(:), N, N);

    if Q^S.Nt > 1e6
        warning('ML search over %d hypotheses per group: simulation will be slow.', Q^S.Nt);
    end
end
