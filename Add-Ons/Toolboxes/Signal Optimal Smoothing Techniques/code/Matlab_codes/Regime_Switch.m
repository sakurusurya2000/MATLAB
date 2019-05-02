function [nreg,lookall,cldvs,cocos,met] = Regime_Switch(X,verb,opreg)
%
% SEQUENTIAL AND CATEGORICAL REGIME SWITCH TESTING (RST)
%
%  Formula: [nreg,lookall,cldvs,cocos,met] = Regime_Switch(X,verb,opreg)
%
% Inputs: (1) signal X:(T,1);  
% (2) verbose definition of series declaring number of sequential regime switches (0 to cnum).
% opreg is option for sequential subsample variance/total sample variance (1) or sequential variance ratio between neighboring subsamples (2).
%
% Output:  
% Sequential RST:
% nreg: number of regimes detected, can be zero or more
% lookall table reporting cluster dates of occurrence, value of X at median cluster period and corresponding value of variance of X, and p-values of test statistics
% regdates: dates at which reime switches occur, if any. Returns [] if switches are zero
% clvds: (2,1) matrix of cluster dates (row 1) and values (row 2)  corresponding to F-test statistics whose pval <=.05, 
% cocos: time-break dates and values obtained by code mulbreak.m
% met: method 1 is F-test, method 2 is Ansari-Bradley
% Categorical RST:
% poster: Gaussian posterior probability distribution of clusters, assumes as H0: 1 regime vs 2 regimes
% QLR: Quasi Likelihood Ratio, Cho and White 2007.
%
T=size(X,1);
Fcrit=.05; % F-test pval under which H0 cannot be rejected
nlincrit=5.99; % 5% critical value (CV) for testing linearity, if measured statistic is >than CV, H0 of linearity cannot be rejected, leading to acceptance of H1 of nonlinearity.  
[~,ns] = linstat(X,0,0,0);nots=ns(1); % Harvey nonlinearity test, nots~chisq(2);
[~,~,lookvar] = mulbreaks(X);
L1=lookvar(:,1);L2=lookvar(:,2);
bdates= (L2(L1>0.05))'; % Fixing breakdates which occur for local variance > variance of X by at least 5%
td=[bdates,T];
[~,tcs,td1] = tchunks(td);if size(td1,1)==1;td1=td1';end
num=numel(td1);
%
if opreg==1; % Comparing tchunk variance to entire sample variance ==> F test and its pvalue in case of linearity, Ansari-Bradley test for variance equality in case of nonlinearity 
 if num>=1
  if nots<=nlincrit
    for i=1:num;nn(i)=numel(tcs(i).a);Variance(i)=var(X(tcs(i).a));vv(i)=var(X(tcs(i).a))/var(X);Pvals(i)=1-cdf('f',vv(i),nn(i),T-1);Median(i)=median(X(tcs(i).a));end
    nups=sum(Pvals<Fcrit);if nups>0;nreg=nups;else;nreg=0;end;met=1; % Method 1 is F-test
    look=[td1,Median',sqrt(Variance'),Pvals'];lookall=look(1:end-1,:);
  else
    for i=1:num;[h,Pv]=ansaribradley(X(tcs(i).a),X);Variance(i)=var(X(tcs(i).a));Pvals(i)=Pv;Median(i)=ceil(median(tcs(i).a));hh(i)=h;end
    nups=sum(hh);if nups>0;nreg=nups;else;nreg=0;end;met=2; % Method 2 is Ansari-Bradley
    look=[td1,Median',sqrt(Variance'),Pvals'];lookall=look(1:end-1,:);
  end
    cldates=lookall(lookall(:,4)<=Fcrit);L22=lookall(:,2);clvals=L22(lookall(:,4)<=Fcrit);cldvs=[cldates(1:end),clvals(1:end)]; nreg=(numel(cldvs))/2;
 else 
     cldates=0;clvals=0;cldvs=0;nreg=0;
 end
end
%
if opreg==2; % Comparing tchunk variance sequentially, F & Ansari-Bradley tests as above
  if num>1;
 if nots<=nlincrit
    for i=1:num-1;nn(i)=numel(tcs(i).a);Variance(i)=var(X(tcs(i).a));vv(i)=var(X(tcs(i).a))/var(X(tcs(i+1).a));Pvals(i)=1-cdf('f',vv(i),nn(i),T-1);Median(i)=ceil(median(tcs(i).a));end
    nups=sum(Pvals<Fcrit);if nups>0;nreg=nups;else;nreg=0;end;met=1; % Method 1 is F-test
    lookall=[td1(1:end-1),Median',sqrt(Variance'),Pvals'];
 else
    for i=1:num-1;[h,Pv]=ansaribradley(X(tcs(i).a),X(tcs(i+1).a));Variance(i)=var(X(tcs(i).a));Pvals(i)=Pv;Median(i)=ceil(median(tcs(i).a));hh(i)=h;end
    nups=sum(hh);if nups>0;nreg=nups;else;nreg=0;end;met=2; % Method 2 is Ansari-Bradley
    lookall=[td1(1:end-1),Median',sqrt(Variance'),Pvals'];
 end
   cldates=lookall(lookall(:,4)<=Fcrit);L22=lookall(:,2);clvals=L22(lookall(:,4)<=Fcrit);cldvs=[cldates(1:end),clvals(1:end)];nreg=(numel(cldvs))/2;
  else
     cldates=0;clvals=0;cldvs=0;nreg=0;
  end
end
%
if verb==1;if isempty(cldvs)==0;disp(['series has   ' num2str(nreg)    '   regime switches']);else;disp(['series has   ' num2str(0) '   regime switches']);end;end
[tbdates,tbvals] = mulbreaks(X);
cocos=[tbdates;tbvals];
%
%
%
end
%