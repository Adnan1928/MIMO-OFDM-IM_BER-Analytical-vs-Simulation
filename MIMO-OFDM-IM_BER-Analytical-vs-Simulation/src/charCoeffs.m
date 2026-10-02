function ek = charCoeffs(A)
% CHARCOEFFS  Elementary symmetric polynomials e_0..e_N of the eigenvalues
%   of every page of A (N x N x P), so that det(I + q*A) = sum_k e_k q^k.
    [N, ~, P] = size(A);
    dIdx = 1:N+1:N^2;
    pk = zeros(N, P);                                   % traces of A^k
    Ak = A;
    for k = 1:N
        Av = reshape(Ak, N^2, P);
        pk(k, :) = real(sum(Av(dIdx, :), 1));
        if k < N, Ak = pageMul(Ak, A); end
    end
    ek = zeros(N+1, P); ek(1, :) = 1;
    for k = 1:N
        acc = zeros(1, P);
        for i = 1:k
            acc = acc + (-1)^(i-1) * ek(k-i+1, :) .* pk(i, :);
        end
        ek(k+1, :) = acc / k;
    end
    ek = max(ek, 0);                                    % PSD -> e_k >= 0
end
