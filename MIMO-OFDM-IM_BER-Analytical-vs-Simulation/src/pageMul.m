function C = pageMul(A, B)
% PAGEMUL  Page-wise matrix product (uses pagemtimes when available).
    if exist('pagemtimes', 'builtin') || exist('pagemtimes', 'file')
        C = pagemtimes(A, B);
    else
        [N, ~, P] = size(A);
        C = reshape(sum(reshape(A, N, N, 1, P) .* reshape(B, 1, N, N, P), 2), N, N, P);
    end
end
