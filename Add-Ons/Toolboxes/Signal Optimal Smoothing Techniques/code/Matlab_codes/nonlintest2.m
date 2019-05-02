function outs1= nonlintest2(y,dops)
%
% Testing for nonlinearity in mean & variance of univariate time series y. Uses my olsols.m routine.
%  Function is:  outs1= nonlintest2(y,dops)
% WARNING: if constant option is explicitly included in ols.m, the coefficient/tstat output puts in the first position
% REFERENCE 1: Harvey, D. I., Leybourne, S. J., and Xiao, B., 2008, "A powerful test for linearity when the order of integration is unknown", 
% Studies in Nonlinear Dynamics & Econometrics, 12, Article 2. Available at: http://www.bepress.com/snde/vol12/iss3/art2 [HLX]
% REFERENCE 2: Terasvirta, T., 2011, Nonlinear models for autoregressive conditional heteroskedasticity, CREATES Research Paper 2011-2, Aarhus University.
% Function: outs1=nonlintest2(y,dops)
%  INPUT: y:(Tx1) time series, dops: option for retaining original series (0), for demeaning (1) or for demeaning and detrending (2)
%  OUTPUT: HLX weighted WALD statistic for null hypothesis of linearity, TERASVIRTA test for GARCH stationarity
%
T=size(y,1);
trend=(1:T)';
lags=3;
cop=1;
a = polyfit(trend,y, 1);yhat=a(1)*trend+a(2); % Fitting straight line to x
if dops==1;y=y-mean(y);elseif dops==2;y=y-mean(y)-yhat;end
% 
%
% A. HLX testing for null of linearity of series y,  H0: y(t)=a+b*y(t-1)+u(t); H1:  y(t)=a+b*y(t-1)+c*y(t-2)^2+d*y(t-3)^3+e(t), a cubic-lag function
%  Test for c,d==0, such that H0 is linearity of y. Computed statistic
%  WALCOM is distributed as Chisq(2). If pval of WALCOM is <5% the null of linearity is rejected, series is significantly nonlinear
%
% 1. Testing for linearity of series expressed in levels, assumes y=I(0)
yfill=lagmatfill(y,lags,1,ceil(T/10));
yfills=[(yfill(:,2)).^2,(yfill(:,3)).^3,yfill(:,1)]; % Form linear, quadratic & cubic level single lags
R=yfills(:,3);[~,~,S2,]=olsols(y,R,cop);SS0=S2;
R=yfills;[~,~,S2,]=olsols(y,R,cop); SS1=S2; % Unrestricted OLS regression=H1: y(t)=a1*y(t-1)+a2*[y(t-1)]^2+a3*[y(t-1)]^3+constant+e(t)
WALDS=T*((SS0/SS1)-1)  ;%outfc(2,1); % Wald-type test
% 2. Testing for linearity of series expressed in differences, assumes y=I(1)
dy=diff(y);
dyfill=lagmatfill(dy,lags,1,ceil(T/10));
dyfills=[(dyfill(:,2)).^2,(dyfill(:,3)).^3,dyfill(:,1)]; %  Form linear, quadratic & cubic first-difference single lags
R=dyfills(:,3);[~,~,S2] = olsols(dy,R,cop); SS0=S2;%
R=dyfills;[~,~,S2] = olsols(dy,R,cop); SS1=S2; % Unrestricted OLS regression=H1: Dy(t)=a1*Dy(t-1)+a2*[Dy(t-1)]^2+a3*[Dy(t-1)]^3+constant+e(t)
% outfc = FCstats(B,Q,S2,T);
WALDU=T*((SS0/SS1)-1); % outfc(2,1); % Wald-type test
% Putting 1, 2 together and form weighted Wald statistic
lam=exp(-.1*(WALDU/WALDS)^2);
WALCOM= abs((1-lam)*WALDS+lam*WALDU);
%
% TERASVIRTA test for asymmetric GARCH
%
% p=1;q=1;  modelg = garch('Offset',NaN,'GARCHLags',p,'ARCHLags',q,'Distribution','Gaussian'); 
% fit = estimate(modelg,y,'print',false);0= infer(fit,y);  et=y-mean(y);
% e2=(et).^2;ev2=e2./V0;ev21=(ev2-1).^2;SSR0=sum(ev21);
% LL=-.5*ev2;GL=gradient(LL); % Gradient of timely Log Likelihood
% etp=et.*(et>0);etn=et.*(et<0); % Asymmetries
% R=[GL,etp,etn];[~,~,S2,]=olsols(ev2,R,cop);SSR1=S2;
%TTEST=((SSR0/SSR1)-1);
% outfc = FCstats(B,1,S2,T);TTEST=outfc(2,1); % Wald-type test
%
% wnames={'Composite Wald','TERASVIRTA'}; % Names 
wnames={'Composite'}; % Name
wtests=[WALCOM];
outs1= [wnames(1);num2cell(wtests)];
%
%
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [B,TS,S2,r] = olsols(y,x,cop)
% OLS estimation of vector y:(T,1) against regressor matrix x:(T,k) that may or may not include constatnt term, depending on option (cop=1)
% Function is: [B,TS,S2] = olsols(y,x,cop);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
