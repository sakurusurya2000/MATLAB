function [B,TS,S2,r] = olsols(y,x,cop)
% OLS estimation of vector y:(T,1) against regressor matrix x:(T,k) that may or may not include constatnt term, depending on option (cop=1)
% Outputs: coefficients (B), t-statistics (TS), squared S.E.s (S2), and vector of disturbances
% Function is: [B,TS,S2,r] = olsols(y,x,cop);
% 
T=numel(y);
if cop==1;x=[x,ones(T,1)];end
xx=x'*x; xy=x'*y;
Pop1=pinv(xx);
Pop=Pop1-(Pop1-Pop1')./2.;
if cond(Pop)<=Inf;Pop=update_sigma(Pop);end
B=Pop*xy;
r=y-x*B; 
S2=(r'*r)/(T-size(xx,2));
den=sqrt(S2)*sqrt(diag(pinv(xx)));
TS=B./den;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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