function [outcon,mem,dec]= conper(y,dmt)
%
% Conditional persistence (Kapeitanos, 2002), a check of speed of decay of autocorr. If series is I(1), I(0) the coefficient is expected to be 
% small (fast decay) or large (slow decay)
% Function is [outcon,dec] = conper(y,dmt)
% INPUT: y:(T,1) signal; dmt: option for raw series (1), demeaned series (2), detrended series (3)
% OUTPUT: outcon:(.1*T,3) matrix of individual shares of lagged y's, cumulated shares & lags ranked, 
% mem=mean memory lags, dec: autocorrelation lag decay 
%
if size(y,1)==1;y=y';end; % Make sure endogenous variable is T x 1 vector
T=size(y,1);trend=(1:T)';
length=floor(.1*T);
a = polyfit(trend,y, 1);yhat=a(1)*trend+a(2); % Fitting straight line to x
if dmt==2;y=y-mean(y);elseif dmt==3;y=yhat;end; % dmtrend code is from Shimotsu (see fric.m)
LM=lagmatfill(y,floor(T/2),1,length);V=LM'*LM/T; 
[s,d]=svd(V);as=abs(s);for i=1:floor(T/2);[~,ii]=max(as(:,i));is(i)=ii;end;
dd=diag(d);shares=(dd(1:length)/sum(dd));cshares=cumsum(shares);lags=(is(1:length))';
outcon=[shares cshares lags];
mem=floor(mean(lags));
acf=autocorr(y);dec=ERE(acf);
%
%
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yfill = lagmatfill(y,q,m1,m2)
%   Code for filling lagmatrix with past means by preserving length of original vector series y.
%   Guido Travaglini, December 1/11
%   Function is yfill = lagmatfill(y,lags,m1,m2)
%
%   INPUTS: 
%   y: original series of size T x 1
%   q: lags in lagmatrix =>1
%   m1, m2: initial observations of y wherefrom compute mean, m2>m1=>2
% 
%  OUTPUT:
%  yfill: filled-in lagmatrix of size T x q
%
T=length(y); yfill=[y,zeros(T,q)];
for j=1:q
 ylak=lagmatrix(yfill(:,j),1);yfill(1,j+1)=mean(yfill(m1:m2,j));
 yfill(2:end,j+1)=ylak(2:end,:);
end
yfill=yfill(:,2:end);
%
end