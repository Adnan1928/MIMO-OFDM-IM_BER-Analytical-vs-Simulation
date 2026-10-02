function idxHat = mlDetect(Yg, Hg, sys, maxElems)
% MLDETECT  Exhaustive joint ML search over Q^Nt hypotheses per group,
%   vectorised over hypotheses and over a block of groups.
%   Hypothesis layout: dim 3 -> antenna 1, dim 4 -> antenna 2, ...
    [N, G, Nr] = size(Yg);
    Nt = sys.Nt; Q = sys.Q; C = sys.C;
    nH     = Q^Nt;
    gBlk   = max(1, floor(maxElems / (N * nH)));
    idxHat = zeros(G, Nt);
    Cq     = reshape(C, N, 1, Q);                       % N x 1 x Q

    for g0 = 1:gBlk:G
        gs  = g0:min(G, g0+gBlk-1);
        met = 0;
        for r = 1:Nr
            R = Yg(:, gs, r);                           % residual y - sum H x
            for t = 1:Nt
                HS = Hg(:, gs, r, t) .* Cq;             % N x |gs| x Q
                R  = R - reshape(HS, [N, numel(gs), ones(1, t-1), Q]);
            end
            met = met + sum(abs(R).^2, 1);              % 1 x |gs| x Q x ...
        end
        [~, lin] = min(reshape(met, numel(gs), nH), [], 2);
        idxHat(gs, :) = mod(floor((lin - 1) ./ Q.^(0:Nt-1)), Q) + 1;
    end
end
