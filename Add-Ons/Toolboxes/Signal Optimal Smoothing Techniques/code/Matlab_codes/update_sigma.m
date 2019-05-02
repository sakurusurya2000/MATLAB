function usigma = update_sigma(sigma)
%  Your covariance matrix is non positive definite such that det(sigma)=0?
%   No worry, correct for its eigenvalues!
% INPUT: sigma:    non positive-definite covariance matrix
% OUTPUT: usigma:    positive-definite covariance matrix
%
EPS = 10^-3;
ZERO = 10^-5;
%
if cond(sigma)<=Inf
    % the covariance matrix is not positive definite!
    [v, d] = eig(sigma);
    % set any of the eigenvalues that are <= 0 to some small positive value
    for n = 1:size(d,1)
        if (d(n, n) <= ZERO)
            d(n, n) = EPS;
        end
    end
    % recompose the covariance matrix, now it should be positive definite.
    usigma = v*d*v';
end
%
end

