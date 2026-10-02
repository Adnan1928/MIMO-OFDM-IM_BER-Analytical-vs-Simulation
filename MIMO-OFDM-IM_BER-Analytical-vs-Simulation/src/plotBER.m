function fig = plotBER(snrdB, berAna, berSim, lbl, P)
% PLOTBER  Analytical vs. simulated BER curves.
%   Analytical = solid line + open marker, simulation = dashed line +
%   filled marker; the same colour / marker is used for one configuration.
%   Font size 12, line width 2.
%
%   fig = plotBER(snrdB, berAna, berSim, lbl, P)
%     snrdB  : 1 x nS SNR axis [dB]
%     berAna : nC x nS analytical BER
%     berSim : nC x nS simulated BER
%     lbl    : 1 x nC cell array of configuration labels
%     P      : parameter struct (uses P.snrMode, P.scheme, P.NF, P.G, P.L)
    nC  = size(berAna, 1);
    mk  = {'o','s','d','^','v','>','<','p','h','*'};
    col = lines(nC);
    fig = figure('Color', 'w'); hold on; grid on; box on;
    hA  = []; hS = [];                               % line handles for the legend
    for c = 1:nC
        m = mk{mod(c-1, numel(mk)) + 1};
        hA(c) = semilogy(snrdB, berAna(c,:), ['-'  m], 'Color', col(c,:), ...
                         'LineWidth', 2, 'MarkerSize', 9);
        hS(c) = semilogy(snrdB, berSim(c,:), ['--' m], 'Color', col(c,:), ...
                         'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', col(c,:));
    end
    set(gca, 'YScale', 'log', 'FontSize', 12, 'LineWidth', 1);
    xlabel(snrLabel(P.snrMode), 'FontSize', 12);
    ylabel('Bit Error Rate (BER)', 'FontSize', 12);
    title(sprintf('%s: analytical vs. simulated BER (N_F=%d, G=%d, L=%d)', ...
          P.scheme, P.NF, P.G, P.L), 'FontSize', 12);
    legend([hA hS], [strcat({'Analytical, '}, lbl), strcat({'Simulation, '}, lbl)], ...
           'Location', 'southwest', 'FontSize', 12);
    ylim([1e-6 1]); xlim([snrdB(1) snrdB(end)]);
end
