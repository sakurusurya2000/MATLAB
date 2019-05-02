function [SST,CST,CSI,VEC] = statcorr_1(X,frac,dop)
%
% Code for producing corrected stationarity test statistics from signal X:(T,1)
% It is an applied extension of statcorr.m, and frac is first subsample fraction expressed in percent, e.g. 0.67
% Dop stands for option to use raw series (0), demeaned series (1), detrended series (2), demeaned+detrended series (3), and standardized (4)
% To batch all four dops, given frac, do: for dop=0:3;[SST,CST,CSI,VEC] = statcorr_1(SSN,frac,dop);SCS(dop+1,:)=[SST,CST];end
%
% Function is:  [SST,CST,CSI,VEC] = statcorr_1(X,frac,dop)
% 
T=numel(X);
N=ceil(frac*T);
trend=(1:T)';
a = polyfit(trend,X, 1);LTX=a(1)*trend+a(2); % Fitting straight line to X & obtaining linear trend of X
if dop==1
    X=bsxfun(@minus, X,mean(X));
elseif dop==2
    X=bsxfun(@minus,X,LTX);
elseif dop==3
    X=bsxfun(@minus, X,(mean(X)+LTX));
elseif dop==4
    X=zscore(X);
end
% 
WALL= speckop(X); % Smoothing by Spectral Representation Theorem
ltratio=LTX(T)/LTX(1);
stratio=WALL(T)/WALL(1);
[~,pratio] = ansaribradley(X(1:N),X(N+1:end));
cratio=[cyclex(X(1:N))/cyclex(X(N+1:end))]; % Linear trend ratio, smoothed trend ratio, pvalue subsample, cycle ratio
VEC=[ltratio,stratio,pratio,cratio];
[~,ns]=linstat(X,0,dop,0);
SST=ns(2);
[CST,CSI]= statcorr(SST,VEC);
%
%
%
end

