[README.md](https://github.com/user-attachments/files/32954702/README.md)
# MIMO-OFDM-IM_BER-Analytical-vs-Simulation
MATLAB/Octave BER analysis of MIMO-OFDM with index modulation (OFDM-IM): Monte-Carlo simulation with joint ML detection vs. analytical union bound. Configurable Tx/Rx antennas, subcarriers per group, active subcarriers and modulation order.
# MIMO-OFDM-IM BER — Analytical vs. Simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![MATLAB R2016b+](https://img.shields.io/badge/MATLAB-R2016b%2B-orange)
![GNU Octave 6+](https://img.shields.io/badge/GNU%20Octave-6%2B-0790c0)
![No toolboxes](https://img.shields.io/badge/toolboxes-none%20required-success)

This repository has MATLAB/Octave code for the **bit error rate (BER) of MIMO-OFDM with index modulation (OFDM-IM)**. It compares a Monte-Carlo simulation with joint maximum-likelihood (ML) detection against an analytical union bound. The number of transmit antennas `Nt`, receive antennas `Nr`, subcarriers per group `N` and active subcarriers per group `K` can all be changed. The modulation order `M` can also be changed.

<p align="center">
  <img src="results/BER_comparison.png" width="640" alt="Analytical vs simulated BER for 1x1, 1x2 and 2x2 OFDM-IM">
</p>

## Features

- **Configurable sizes.** Any `Nt`, `Nr`, `N = NF/G`, `1 ≤ K ≤ N`, M-PSK or square M-QAM.
- **Two schemes.** Classic **OFDM-IM**, where idle subcarriers are zero, and **DM-OFDM-IM**, where idle subcarriers carry a rotated second constellation.
- **Realistic channel.** Frequency-selective Rayleigh channel with `L` taps and any power-delay profile, plus block interleaving of the groups.
- **Joint ML detector.** It searches every combination across all transmit antennas at once (`Q^Nt` per group), using array operations instead of loops.
- **Fast analytical bound.** Pairs with the same difference vector are grouped, and each determinant is reduced to a polynomial in SNR. The 2×2 bound takes less than a second.
- **Adaptive Monte-Carlo.** Each SNR point stops once it has a target number of errors.
- **No toolboxes needed.** The Gray PSK/QAM mapper and the bit tables are written in plain code.
- **Tests included.** The tests run in MATLAB and in GNU Octave, and GitHub Actions CI runs them on every push.

## Quick start

```matlab
% in MATLAB or Octave, from the repository root
main_BER_comparison
```

The figure and the data are saved in `results/` (`BER_comparison.png`, `BER_results.mat`).

To choose which curves are compared, edit section 2 of `main_BER_comparison.m`:

```matlab
cfg = struct('Nt', {1, 1, 2}, ...   % transmit antennas
             'Nr', {1, 2, 2}, ...   % receive antennas
             'K',  {2, 2, 2});      % active subcarriers per group
```

Section 1 holds the main parameters:

| Parameter   | Default     | Meaning                                         |
|-------------|-------------|-------------------------------------------------|
| `P.NF`      | 64          | total subcarriers                               |
| `P.G`       | 16          | groups (`N = NF/G` subcarriers per group)       |
| `P.CP`      | 16          | cyclic prefix (must be `≥ L-1`)                 |
| `P.L`       | 10          | channel taps                                    |
| `P.pdp`     | uniform     | power-delay profile (sums to 1)                 |
| `P.M`       | 4           | modulation order                                |
| `P.modType` | `'PSK'`     | `'PSK'` or `'QAM'`                              |
| `P.scheme`  | `'OFDM-IM'` | `'OFDM-IM'` or `'DM-OFDM-IM'`                   |
| `P.snrMode` | `'EbN0'`    | x-axis: `'EbN0'` (CP energy included) or `'EsN0'` |

Section 3 holds the Monte-Carlo settings: `MC.minErr`, `MC.minFrames` and `MC.maxFrames`.

## Repository layout

```
├── main_BER_comparison.m   % entry point: parameters, run, plot
├── src/
│   ├── buildSystem.m       % codebook, bit labels, Hamming table, channel correlation
│   ├── grayConstellation.m % Gray-labelled M-PSK / square M-QAM
│   ├── noisePower.m        % N0 from Eb/N0 or Es/N0
│   ├── simulateBER.m       % Monte-Carlo link simulation
│   ├── mlDetect.m          % joint ML detector across all Tx antennas
│   ├── analyticalBER.m     % union bound on the BER
│   ├── charCoeffs.m        % det(I+qA) as a polynomial in q (Newton identities)
│   ├── pageMul.m           % page-wise matrix product
│   ├── plotBER.m           % figure (font 12, line width 2, markers)
│   └── snrLabel.m
├── tests/run_tests.m       % unit and sanity tests
├── docs/theory.md          % system model and BER derivation
└── results/                % example figure and BER table (CSV)
```

## Example results

Default settings: NF = 64, G = 16, N = 4, K = 2, 4-PSK, L = 10. Numbers below are BER at the given Eb/N0 (full table in [`results/BER_default.csv`](results/BER_default.csv)).

| Eb/N0 | 1×1 analytical | 1×1 simulation | 1×2 analytical | 1×2 simulation | 2×2 analytical | 2×2 simulation |
|------:|------:|------:|------:|------:|------:|------:|
| 10 dB | 5.09e-2 | 2.13e-2 | 1.48e-3 | 9.82e-4 | 2.25e-3 | 1.47e-3 |
| 14 dB | 1.36e-2 | 7.12e-3 | 2.07e-4 | 1.37e-4 | 2.95e-4 | 1.93e-4 |
| 18 dB | 4.32e-3 | 2.65e-3 | 3.22e-5 | 1.82e-5 | 4.54e-5 | 3.33e-5 |

At medium and high SNR the bound has the same slope (diversity order) as the simulation and sits slightly above it. At low SNR it is loose and can go above 1; this is expected for union bounds.

## Theory in brief

For each group, the BER is bounded by the union bound:

$$
P_b \le \frac{1}{N_t\,p\,Q^{N_t}} \sum_{\mathbf X}\sum_{\hat{\mathbf X}\ne\mathbf X} P(\mathbf X\to\hat{\mathbf X})\, e(\mathbf X,\hat{\mathbf X}),
$$

Using the Chiani approximation of the Q-function and the exact correlation $\mathbf K$ between the interleaved subcarriers of a group, each pairwise error probability is

$$
P(\mathbf X\to\hat{\mathbf X}) \approx \tfrac{1}{12}\det(\mathbf I+q_1\mathbf A)^{-N_r}+\tfrac14\det(\mathbf I+q_2\mathbf A)^{-N_r},\qquad
\mathbf A=\mathbf K\odot(\mathbf E\mathbf E^H),
$$

where $\mathbf E=[\mathbf x_1-\hat{\mathbf x}_1,\dots,\mathbf x_{N_t}-\hat{\mathbf x}_{N_t}]$, $q_1 = 1/(4N_0)$ and $q_2 = 1/(3N_0)$. The full derivation is in [docs/theory.md](docs/theory.md).

## Running the tests

```matlab
cd tests
run_tests
```

From the command line with Octave: `octave-cli --eval "cd tests; run_tests"`.

## Complexity notes

- The ML search grows as `2^(Nt·p)` per group, where `p` is the number of bits per group. For example, Nt = 2 with N = 4, K = 4 and 4-PSK already gives 65,536 combinations per group, and the code prints a warning when a setting will be slow.
- The analytical sum grows as `U^Nt`, where `U` is the number of distinct difference vectors (317 for the default). Above `AN.maxTerms` the sum is estimated from random samples instead (the estimate is unbiased).

## Citation

If you use this code, please cite it with the metadata in [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button for it).

## References

1. E. Başar, Ü. Aygölü, E. Panayırcı, H. V. Poor, "Orthogonal frequency division multiplexing with index modulation," *IEEE Trans. Signal Process.*, vol. 61, no. 22, pp. 5536–5549, 2013.
2. E. Başar, "Multiple-input multiple-output OFDM with index modulation," *IEEE Signal Process. Lett.*, vol. 22, no. 12, pp. 2259–2263, 2015.
3. T. Mao, Z. Wang, Q. Wang, S. Chen, L. Hanzo, "Dual-mode index modulation aided OFDM," *IEEE Access*, vol. 5, pp. 50–60, 2017.
4. M. Chiani, D. Dardari, M. K. Simon, "New exponential bounds and approximations for the computation of error probability in fading channels," *IEEE Trans. Wireless Commun.*, vol. 2, no. 4, pp. 840–845, 2003.

## License

[MIT](LICENSE) © 2026 Adnan Tariq
