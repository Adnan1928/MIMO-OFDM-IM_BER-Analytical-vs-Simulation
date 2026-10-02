# Changelog

## [1.0.0] - 2026-10-02
### Added
- Generic MIMO-OFDM-IM simulator: any `Nt`, `Nr`, `N`, `K`, M-PSK / M-QAM.
- DM-OFDM-IM option (idle subcarriers carry a rotated constellation).
- Joint ML detection across all transmit antennas, with array operations.
- Analytical union bound with exact channel correlation inside a group.
- Plot with font 12, line width 2 and markers; results saved to `results/`.
- Unit and sanity tests; GitHub Actions CI with GNU Octave.
