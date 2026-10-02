%% ========================================================================
%  MIMO-OFDM-IM : Simulated vs. Analytical BER  (generic Nt, Nr, N, K, M)
%  ------------------------------------------------------------------------
%  System model (per OFDM frame, V-BLAST style spatial multiplexing):
%    * Every transmit antenna sends an independent OFDM-IM frame of NF
%      subcarriers, split into G groups of N = NF/G subcarriers.
%    * In every group, b = floor(log2(nchoosek(N,K))) index bits select K
%      active subcarriers and K*log2(M) bits are carried by M-ary symbols.
%      (Optional 'DM-OFDM-IM': the N-K idle subcarriers carry symbols of a
%      second, rotated constellation -> N*log2(M) symbol bits per group.)
%    * Block interleaving: subcarrier n of group g goes to OFDM subcarrier
%      (n-1)*G + g, so the subcarriers of one group are G bins apart.
%    * Frequency-selective Rayleigh channel with L taps (power profile pdp)
%      for each of the Nr x Nt links; CP >= L-1, so OFDM turns it into a
%      per-subcarrier flat channel  y = sum_t H_t .* x_t + w.
%    * Joint maximum-likelihood (ML) detection per group over all Nt
%      antennas (exhaustive search over Q^Nt hypotheses, Q = 2^p).
%
%  Analytical BER: union bound with ML detection, exact channel correlation
%  inside a group, and the Chiani approximation of the Q-function:
%     Pb <= 1/(Nt*p*Q^Nt) * sum_X sum_Xhat P(X->Xhat) * e(X,Xhat)
%     P(X->Xhat) ~ 1/12*det(I+q1*Mx)^-Nr + 1/4*det(I+q2*Mx)^-Nr
%     Mx = Kc .* (E*E'),  E = [x_1-xh_1, ..., x_Nt-xh_Nt],  q1=1/(4N0), q2=1/(3N0)
%  The bound is tight at medium/high SNR (it may exceed the simulation at
%  low SNR, which is normal for union bounds).
%
%  Requirements: MATLAB R2016b+ or GNU Octave 6+. No toolbox needed.
%  Helper functions live in ./src (added to the path below).
%
%  Repository : MIMO-OFDM-IM_BER (Analytical Vs Simulation)
%  Author     : Adnan Tariq      License: MIT (see LICENSE)
%  ========================================================================
clear; clc; close all;

%% ------------------------------------------------------------------------
%  0) Paths: make ./src visible and create the output folder
%  ------------------------------------------------------------------------
rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir), rootDir = pwd; end
addpath(fullfile(rootDir, 'src'));
outDir  = fullfile(rootDir, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ------------------------------------------------------------------------
%  1) Common system parameters
%  ------------------------------------------------------------------------
P.NF      = 64;            % total number of subcarriers
P.G       = 16;            % number of groups  -> N = NF/G subcarriers/group
P.CP      = 16;            % cyclic prefix length (must be >= L-1)
P.L       = 10;            % number of channel taps (multipath)
P.pdp     = ones(1,P.L)/P.L;   % power delay profile (sums to 1)
P.M       = 4;             % modulation order (power of 2)
P.modType = 'PSK';         % 'PSK' or 'QAM' (square QAM)
P.scheme  = 'OFDM-IM';     % 'OFDM-IM' (idle = 0) or 'DM-OFDM-IM'
P.snrMode = 'EbN0';        % 'EbN0' : x-axis is Eb/N0 (CP energy included)
                           % 'EsN0' : x-axis is Es/N0 per subcarrier
EbN0dB    = 0:2:20;        % SNR points [dB]

%% ------------------------------------------------------------------------
%  2) Configurations to compare  (one simulated + one analytical curve each)
%     Add / remove entries freely: any Nt, Nr and K (1 <= K <= N) work.
%  ------------------------------------------------------------------------
cfg = struct('Nt', {1, 1, 2}, ...
             'Nr', {1, 2, 2}, ...
             'K',  {2, 2, 2});

%% ------------------------------------------------------------------------
%  3) Monte-Carlo stopping rules (stop when enough errors are collected)
%  ------------------------------------------------------------------------
MC.minErr    = 300;        % target number of bit errors per SNR point
MC.minFrames = 100;        % always simulate at least this many frames
MC.maxFrames = 2e4;        % hard limit of frames per SNR point
MC.maxElems  = 4e6;        % memory limit for the vectorised ML search
AN.maxTerms  = 5e7;        % above this, the union sum is sampled (unbiased)
rng(1);                    % reproducible results

%% ------------------------------------------------------------------------
%  4) Run simulation and analysis for every configuration
%  ------------------------------------------------------------------------
nC     = numel(cfg);
berSim = nan(nC, numel(EbN0dB));
berAna = nan(nC, numel(EbN0dB));
lbl    = cell(1, nC);

for c = 1:nC
    S      = P;                      % copy common parameters
    S.Nt   = cfg(c).Nt;
    S.Nr   = cfg(c).Nr;
    S.K    = cfg(c).K;
    sys    = buildSystem(S);         % codebook, bit tables, correlations
    lbl{c} = sprintf('%d\\times%d, N=%d, K=%d, %d-%s', ...
                     sys.Nt, sys.Nr, sys.N, sys.K, sys.M, sys.modType);

    fprintf('\n=== Config %d/%d : Nt=%d Nr=%d N=%d K=%d M=%d  (p=%d bits/group, %d ML hypotheses)\n', ...
            c, nC, sys.Nt, sys.Nr, sys.N, sys.K, sys.M, sys.p, sys.Q^sys.Nt);

    tic; berAna(c,:) = analyticalBER(sys, EbN0dB, AN);
    fprintf('  analytical bound done in %.1f s\n', toc);
    tic; berSim(c,:) = simulateBER(sys, EbN0dB, MC);
    fprintf('  simulation done in %.1f s\n', toc);
end

%% ------------------------------------------------------------------------
%  5) Plot (font 12, line width 2) and save results
%  ------------------------------------------------------------------------
fig = plotBER(EbN0dB, berAna, berSim, lbl, P);

save(fullfile(outDir, 'BER_results.mat'), 'EbN0dB', 'berSim', 'berAna', 'cfg', 'P');
saveas(fig, fullfile(outDir, 'BER_comparison.png'));
fprintf('\nResults saved in %s\n', outDir);
