function run_tests()
% RUN_TESTS  Unit and sanity tests (MATLAB and GNU Octave).
%   From the repository root:   cd tests; run_tests
%   Throws an error (non-zero exit in CI) if any test fails.
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'src'));
    rng(7);

    tests = {@test_constellations, @test_codebook, @test_dm_codebook, ...
             @test_char_coeffs, @test_ml_noiseless, @test_ber_agreement};
    nFail = 0;
    for k = 1:numel(tests)
        name = func2str(tests{k});
        try
            tests{k}();
            fprintf('[PASS] %s\n', name);
        catch err
            nFail = nFail + 1;
            fprintf('[FAIL] %s : %s\n', name, err.message);
        end
    end
    if nFail > 0
        error('%d of %d tests failed.', nFail, numel(tests));
    end
    fprintf('All %d tests passed.\n', numel(tests));
end

% -------------------------------------------------------------------------
function P = baseParams()
% Default parameters shared by the tests.
    P.NF = 64; P.G = 16; P.CP = 16; P.L = 10; P.pdp = ones(1,10)/10;
    P.M = 4; P.modType = 'PSK'; P.scheme = 'OFDM-IM'; P.snrMode = 'EbN0';
    P.Nt = 1; P.Nr = 1; P.K = 2;
end

function test_constellations()
% Unit energy and Gray labelling (neighbouring PSK points differ by 1 bit).
    for M = [2 4 8 16]
        A = grayConstellation(M, 'PSK');
        assert(abs(mean(abs(A).^2) - 1) < 1e-12, 'PSK energy');
        [~, ord] = sort(mod(angle(A), 2*pi));          % points by angle
        lab = ord - 1;
        d = sum(dec2bin(bitxor(lab, circshift(lab, 1))) == '1', 2);
        assert(all(d == 1), 'PSK Gray labelling');
    end
    for M = [4 16 64]
        A = grayConstellation(M, 'QAM');
        assert(abs(mean(abs(A).^2) - 1) < 1e-12, 'QAM energy');
        assert(numel(unique(round(A*1e9))) == M, 'QAM points distinct');
    end
end

function test_codebook()
% OFDM-IM codebook: size, K active subcarriers, distinct codewords, energy.
    for K = 1:4
        P = baseParams(); P.K = K;
        s = buildSystem(P);
        assert(isequal(size(s.C), [s.N, s.Q]), 'codebook size');
        assert(s.p == floor(log2(nchoosek(s.N, K))) + K*log2(s.M), 'bits/group');
        assert(all(sum(abs(s.C) > 1e-12, 1) == K), 'K active subcarriers');
        assert(size(unique(round([real(s.C); imag(s.C)].' * 1e9), 'rows'), 1) == s.Q, 'distinct codewords');
        assert(abs(s.Eg - K) < 1e-12, 'group energy');
        assert(all(diag(s.D) == 0) && isequal(s.D, s.D.'), 'Hamming table');
    end
end

function test_dm_codebook()
% DM-OFDM-IM codebook: all subcarriers used, extra bits on idle ones.
    P = baseParams(); P.scheme = 'DM-OFDM-IM';
    s = buildSystem(P);
    assert(s.p == 2 + 4*2, 'DM bits/group');
    assert(all(all(abs(s.C) > 1e-12)), 'all subcarriers active');
    assert(size(unique(round([real(s.C); imag(s.C)].' * 1e9), 'rows'), 1) == s.Q, 'distinct DM codewords');
end

function test_char_coeffs()
% det(I + q*A) from the characteristic coefficients equals the direct det.
    N = 4; P = 20;
    A = randn(N, N, P) + 1j*randn(N, N, P);
    for k = 1:P, A(:,:,k) = A(:,:,k) * A(:,:,k)'; end
    ek = charCoeffs(A);
    q  = 0.37;
    d1 = ek.' * (q.^(0:N)).';
    d2 = zeros(P, 1);
    for k = 1:P, d2(k) = real(det(eye(N) + q*A(:,:,k))); end
    assert(max(abs(d1 - d2) ./ d2) < 1e-9, 'determinant mismatch');
end

function test_ml_noiseless()
% Without noise the joint ML detector must recover every codeword.
    P = baseParams(); P.Nt = 2; P.Nr = 2;
    s = buildSystem(P);
    G = s.G; N = s.N;
    idx = randi(s.Q, G, s.Nt);
    X   = reshape(s.C(:, idx), N, G, s.Nt);
    Hg  = (randn(N, G, s.Nr, s.Nt) + 1j*randn(N, G, s.Nr, s.Nt)) / sqrt(2);
    Yg  = zeros(N, G, s.Nr);
    for r = 1:s.Nr
        for t = 1:s.Nt
            Yg(:,:,r) = Yg(:,:,r) + Hg(:,:,r,t) .* X(:,:,t);
        end
    end
    idxHat = mlDetect(Yg, Hg, s, 4e6);
    assert(isequal(idxHat, idx), 'noiseless ML detection failed');
end

function test_ber_agreement()
% Union bound and simulation agree at medium SNR (1x2, 10 dB).
    P = baseParams(); P.Nr = 2;
    s  = buildSystem(P);
    MC = struct('minErr', 200, 'minFrames', 100, 'maxFrames', 5000, 'maxElems', 4e6);
    AN = struct('maxTerms', 5e7);
    ba = analyticalBER(s, 10, AN);
    bs = simulateBER(s, 10, MC);
    ratio = ba / bs;
    assert(ratio > 0.8 && ratio < 3, sprintf('bound/simulation ratio %.2f out of range', ratio));
end
