function [ls,ns,crit] = linstat(X,rep,rop,cor)
%
%   Testing by 5% pvalue for linearity vs nonlinearity and for stationarity vs nonstationarity of signal X
%   Function: [ls,ns,crit] = linstat(X,rep,rop,cor)
%
%   INPUT: signal X, rep: option of verbose reply, rop: option in stars2 to
%   use  raws series (0), demeaned series (1), detrended series (2), standardized (3)
%   rop is   used in stars2,  and option for corrected stationarity test statistic: cor=1
%
%   OUTPUT: 
%  (1) ls: (2,1) reporting index of (non)linearity (1=linear, 2=nonlinear) and nonstationarity (2) or stationarity (1) at 5% pvalue 
%  (2) ns:  Harvey nonlinearity test, nots~chisq(2) and  ADF/KSS stationarity test at 5% in the presence of linearity/nonlinearity
%  (3) crit: critical Unit Root (UR) test value at 5% for ADF/KSS 
%
if size(X,1)==1;X=X';end; % Make sure signal X is entered columnwise
T=numel(X);
trend=(1:T)';
a = polyfit(trend,X, 1);LTX=a(1)*trend+a(2); % Fitting straight line to x & obtaining linear trend of X
if rop==1
    X=bsxfun(@minus, X,mean(X));
elseif rop==2
    X=bsxfun(@minus,X,LTX);
elseif rop==3
    X=zscore(X);
end
nlincrit=5.99; ksscrit=-2.91; % Fixing 5% critical values beyond which we have nonlinearity and nonstationarity
outs1= nonlintest2(X,rop);nots=cell2mat(outs1(2)); % Harvey nonlinearity test, nots~chisq(2);
star1= stars2(X,rop); % KSS stationarity test in the presence of nonlinearity
%
% Testing for stationarity in linear/nonlinear series; % 5.99 is Chisq(2) at 95%
% If signal is linear, such that nots<5.99, use ADF test for stationarity, elseif signal is nonlinear, such that nots>5.99, use KSS test for stationarity 
if  nots<nlincrit; 
    [~,p,ADF,CDF]=adftest(X); crit=CDF; % Critcal UR 5% ADF
    nas= abs(ADF)-abs(CDF);
    if nas<0;CADF=2;else;CADF=1;end; % CADF is index that establishes whether ADF is nonstationary (2) or stationary (1) at 5% pvalue 
    NSTEST=[num2cell(1),num2cell(CADF),num2cell(p)];ns=[nots,num2cell(ADF)];
else
    KSS=cell2mat(star1(4,1)); % if abs(KSS)> critvalue, X is stationary ==> p=1
    if abs(KSS)>abs(ksscrit);p=1;else;p=0;end; crit= -2.91; % Critcal UR 5% KSS
    kas=abs(KSS)-abs(ksscrit);
    if kas<0;KASS=2;else;KASS=1;end; % KASS is index that establishes whether KSS is nonstationary (2) or stationary (1) at 5% pvalue 
    NSTEST= [num2cell(2),num2cell(KASS),num2cell(p)];ns=[nots,num2cell(KSS)];
end
ls=cell2mat(NSTEST(1:end-1));ns=cell2mat(ns);
if cor==1
    [~,~,lookvar] = mulbreaks(X);N=lookvar(1,2);frac=N/T;ts.frac=N;
    [~,CST] = statcorr_1(X,frac,rop);
    ns(2)=CST;
    if CST<ksscrit;ls(2);ls(2)=1;end
end
if rep==1;
 if ls(1)==1;disp('series is linear');else;disp('series is nonlinear');end;
 if ls(2)==0;disp('series is nonstationary');else;disp('series is stationary');end;
end
%
%
%
end

