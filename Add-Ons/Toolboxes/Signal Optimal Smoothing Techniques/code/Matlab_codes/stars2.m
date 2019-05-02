function star1= stars2(y,dops)
%
% Testing for stationary ESTAR-AESTAR alternatives to linear UR: KSS ADF-type and  SOLLIS F-type tests. Uses my olsols.m routine.
% Has "dops" option to demean series (1), detrend (2), demean + detrend  (3);
% Uses Modified Akaike Information Criterion (MAIC), Ng, S., and P. Perron (2001): 
% Lag Length Selection and the Construction of Unit Root Tests with Good Size and Power, Econometrica 69, 1519-1554.
%
% REFERENCE 1: Kapetanios, G., Shin, Y., Snell, A., 2003, Testing for a Unit Root in the Non-linear STAR Framework, Journal of Econometrics 112: 359-379 [KSS]
% REFERENCE 2: Sollis, R., 2009, A Simple Unit Root Test Against Asymmetric STAR Nonlinearity With An Application To Real Exchange Rates In Nordic Countries,
% Economic Modelling 26, (2009): 118125 [SOLLIS]
% REFERENCE 3: Shintani, M., 2013, The inf-t test for a unit root against asymmetric ESTAR models, Japanese Economic Review, 64(1), 3-15 [SHINT]
%
%  INPUT: y:(Tx1) time series, dops: option for retaining original series (0), for demeaning (1), for detrending (2), for demeaning and detrending (3)
%  OUTPUT: KSS & SOLLIS tests for stationarity of respective STAR model
%
T=size(y,1);
trend=(1:T)';
maxlag=ceil(12*(ceil(T/100))^.25); % Rule fixed by Perron & Qu, 2006, A Simple Modification to Improve the Finite Sample Properties of Ng and Perrons Unit Root Tests
cop=1;
% 
% Demeaned and detrended original variable (raw data), required in both KSS and SOLLIS
%
a = polyfit(trend,y, 1);yhat=a(1)*trend+a(2); % Fitting straight line to x
if dops==1;y=y-mean(y);elseif dops==2;y=y-yhat;elseif dops==3;y=y-mean(y)-yhat;end
dy=diff(y);dyfillags=lagmatfill(dy,maxlag,1,ceil(T/10));
y2=y.^2;y3=y.^3;y4=y.^4;
%
%
% Testing for UR in the presence of nonlinear adjustments, KSS: untransformed first differences, conventional DF with single lag of level & symmetric ESTAR, basically a DF/ADF test with
% cubed lagged level & lagged first-differences in order to detect nonstationarity vs globally stationary exponential smooth transition autoregressive (ESTAR) nonlinearity,
%
for k=1:maxlag;
    R=[y3(1:end-1),dyfillags(:,1:k)];
    [B,~,~,r]=olsols(dy,R,cop);
    sig2=sum(r.^2)/(T-maxlag);
    LS=log(sig2);
    ttk=((B(1)^2)/sig2)*sum(y2(maxlag+1:end));
    MC=2*(ttk+k)/(T-maxlag);
    LOGL(k)=LS+MC;
end
%
[~,maiclag]=min(LOGL);
R=dyfillags(:,1:maiclag);[B,~,S2]=olsols(dy,R,cop);SS0=S2;L0= -.5*T*(1+log(2*pi)+log(S2/T)); DF=numel(B); % H0 stationary model
R=[y3(1:end-1),dyfillags(:,1:maiclag)];[~,TS,S2]=olsols(dy,R,cop);SS1=S2;L1= -.5*T*(1+log(2*pi)+log(S2/T)); % H0 symmetric ESTAR model
KADF=TS(1);  % ADF test of cubed lagged untransformed
WALDK=T*((SS0/SS1)-1);  %outfc(2,1); % Wald-type test
% [~,PKU,LRKU] = lratiotest(L1,L0,DF);
%
 % SOLLIS : dy(t)=a0+a1*y3(t-1)+a2*y4(t-1)+b(j)*Sum(dy(t-j))+v1(t), null of lags only model vs alternative of asymmetric ESTAR
%
for k=1:maxlag;
    R=[y3(1:end-1),y4(1:end-1),dyfillags(:,1:k)];
    [B,~,~,r]=olsols(dy,R,cop);
    sig2=sum(r.^2)/(T-maxlag);
    LS=log(sig2);
    ttk=((B(1)^2)/sig2)*sum(y2(maxlag+1:end));
    MC=2*(ttk+k)/(T-maxlag);
    LOGL(k)=LS+MC;
end
%        
[~,maiclag]=min(LOGL);
R=dyfillags(:,1:k);[B,~,S2]=olsols(dy,R,cop);SS0=S2;DF=numel(B); L0= -.5*T*(1+log(2*pi)+log(S2/T)); %H0 first-differenced lags only 
R=[y3(1:end-1),y4(1:end-1),dyfillags(:,1:maiclag)];[B,~,S2]=olsols(dy,R,cop);SS1=S2; L1= -.5*T*(1+log(2*pi)+log(S2/T)); DF=numel(B);%H1 lagged levels & first-differenced lags 
% outfc = FCstats(B,2,S2,T);FSOL=outfc(1,1);WALDS=outfc(2,1); % Wald- and F-type tests
WALDS=T*((SS0/SS1)-1); FSOL=WALDS/2; % F-test computed as WALDS/2, where 2 is the # restrictions
% [~,PSU,LRSU] = lratiotest(L1,L0,DF);  %  LLRatio test, store pvalue and computed statistic
%
%
% SHINT: dy(t)=a0+a1*y2(t-1)+b(j)*Sum(dy(t-j))+v1(t), null of UR vs alternative of  ESTAR ==> test of symmetric and nonlinear mean reverting properties.
%
K=-10;
p=mean(y2)^(-.5);yde=K*sqrt(y2).*(1-exp(-p^2.*y2)); % Eq. (4), Reference 3
R=yde(1:end-1);B=olsols(dy,R,cop); DF=numel(B); 
for k=1:maxlag;
    R=[yde(2:end),dyfillags(:,1:k)];
    [B,~,~,r]=olsols(dy,R,cop);
    sig2=sum(r.^2)/(T-maxlag);
    LS=log(sig2);
    ttk=((B(1)^2)/sig2)*sum(y2(maxlag+1:end));
    MC=2*(ttk+k)/(T-maxlag);
    LOGL(k)=LS+MC;
end
%    
[~,maiclag]=min(LOGL);
R=dyfillags(:,1:maiclag);[~,~,S2]=olsols(dy,R,cop);SS0=S2;  L0= -.5*T*(1+log(2*pi)+log(S2/T)); % H0: UR nonlinear model
R=[yde(2:end),dyfillags(:,1:maiclag)];[~,TS,S2]=olsols(dy,R,cop); SS1=S2;L1= -.5*T*(1+log(2*pi)+log(S2/T)); %H1 ESTAR model
ITEST=TS(1); % SHINT Inf-t-test. Cvs are tabulated in Reference 3
WALSH=T*((SS0/SS1)-1);  % Wald-type test
% [~,PSH,LRSH] = lratiotest(L1,L0,DF);  %  LLRatio test, store pvalue and computed statistic
%
% Critical values for untransformed (raw) and transformed series KSS
% t-tests, SOLLIS F-tests & for Shintani inf-test; 200<T=500
% UN=untransformed, DM=Demeaned only, DT= Demeaned + Detrended.
% KSS symmetric ESTAR at 1% and 5%: UN(-2.79,-2.51), DM(-3.48, -2.93), DT(-3.93, -3.40)
% Sollis AESTAR at 1% and 5%: UN(4.241, 2.505), DM(6.236, 4.557), DT(8.344, 6.292)
% Shintani symmetric ESTAR at 1% and 5%: UN(-2.85, -2.27), DM(-3.71, -3.13), DT (-4.18, -3.64)
% 
wnames={'ESTAR KSS','AESTAR SOLLIS','AESTAR SHINT'}; % Names of Wald-type tests for Nonlinearity
wtests=[WALDK,WALDS,WALSH];
unames={'KSS ADF','SOLLIS F-TEST','SHINT Inf-t-test'}; % Names of UR tests for Nonstationarity
utests=[KADF,FSOL,ITEST];
% ={'KSS LR','SOLLIS LR','SHINTANI LR'}; lrtests=[LRKU, LRSU,LRSH;PKU,PSU,PSH]; % LLR computed stats and pvalues
%
star1= [wnames;num2cell(wtests);unames;num2cell(utests)]; % star2= [lrnames;num2cell(lrtests(1,:))];
%
%
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,TS,S2,r] = olsols(y,x,cop)
% OLS estimation of vector y:(T,1) against regressor matrix x:(T,k) that may or may not include constant term, depending on option (cop).
% If cop=1, constant term features the last in B coefficient vector.
% 
T=numel(y);
if cop==1;x=[x,ones(T,1)];end
xx=x'*x; xy=x'*y;
B=pinv(xx)*xy;
r=y-x*B; 
S2=(r'*r)/(T-size(xx,2));
den=sqrt(S2)*sqrt(diag(pinv(xx)));
TS=B./den;
%
end
