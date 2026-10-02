function A = grayConstellation(M, type)
% GRAYCONSTELLATION  Gray-labelled M-PSK or square M-QAM, unit avg energy.
%   A(l+1) is the point carrying label l.
    k = 0:M-1;
    if strcmpi(type, 'PSK')
        g = bitxor(k, floor(k/2));
        A = zeros(1, M);
        A(g+1) = exp(1j*2*pi*k/M);
    else
        m = sqrt(M);
        assert(mod(m,1) == 0, 'QAM order must be a square (4, 16, 64, ...).');
        kk  = 0:m-1;
        pam = zeros(1, m);
        pam(bitxor(kk, floor(kk/2)) + 1) = 2*kk - m + 1;   % Gray PAM
        A = pam(floor(k/m)+1) + 1j*pam(mod(k,m)+1);
        A = A / sqrt(mean(abs(A).^2));
    end
end
