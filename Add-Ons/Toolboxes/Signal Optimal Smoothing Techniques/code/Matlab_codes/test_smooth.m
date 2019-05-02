function ts= test_smooth( X)
% 
%   Performs all possible tests & smoothers on signal X
%   Function is ts=test_smooth(X);
%   HHT: Hilbert-Huang Transorm; SRT: Spectral Representation Theorem, simple; SRH: Spectral Representation Theorem, Hamiltonian correction
%   Uses one HHT and two outputs of SRT (see speckoph.m)
%  Guido Travaglini, August 25/15
%
% INPUT: 
% vector signal X:(T,1)
%
% OUTPUT: 
% object matrix ts including the following statistics/data:
%
% 1. ts.HHT time series, n=1
% 2. ts.SRT, time series, n=2, first is standard SRT, second is envelope-adjusted
% 3. ts.errs, error time series: X-HHT, X-SRT
% 4. ts.ersmev, mean and variance of error and squared error time series: first row errors, second row squared errors
% 5. ts.errstat, conventional stationarity and corrected stationarity test statistics of errors; rows are HHT and two SRT, first three columns are coventional  ADF/KSS tests, other columns are corrected tests 
% 6. ts,serrs, squared errors time series: (X-HHT)^2, (X-SRT)^2
% 7. ts.linerser,  Harvey value (first row) and index number (second row) of linearity (1), nonlinearity (2) of errors and squared errors
% 8. ts.serrstat, conventional ADF/KSS stationarity and corrected stationarity test statistics of squared errors; rows are HHT and two SRT, first three columns are coventional  tests, other columns are corrected tests 
% 9. ts.prmse, PRMSEs of HHT-X and of SRT-X squared errors
% 10. ts.cross, Number of crossings of HHT, SRT wrt X
% 11. ts.ltrend, Linear trend of X
% 12. ts.normalp, Test for Normality of X, p-value of Kolmogorov-Smirnov test
% 13. ts.linstat1, Index # of Harvey linearity and ADF/KSS stationarity tests of X, (1): linear, nonstationary, (2) nonlinear, stationary
% 14. ts.linstat2, Values of Harvey linearity and ADF/KSS stationarity tests of X
% 15. ts.corrstat, Corrected stationarity test statistic of X
% 16. ts.lratio, Linear trend ratio last-to-first observation
% 17. ts.sratio, Smoother ratios
% 18. ts.pratio, p-value of the F-test (linear X) or Ansari-Bradley test (nonlinear X)
% 19. ts.cycles, FFT cycle of X and three smoothers (frequency)
% 20. ts.cratio, Subsample cycle ratio of X
% 21. ts.conper, Conditional persistence of X
% 22. ts.mimas, Squared value and locational distances of X wrt HHT, SRT attained at respective min/max value and location 
% 23. ts.turnp, # of turning points from phase analysis of X and of three smoothers
% 24. ts.rsnum, Number of Regime Switches
% 25. ts.rsdates, Regime Switch date(s) - if any - sorted in ascending  chronological order and entered as vertical vector
% 26. ts.metsel, prefered smoother method (1 if HHT, 2,3 if SRT)
% 27. ts.xsels, (T,2) Matrix reporting X and minimum PRMSE-selected smoother
% 28. ts.nums, # of time chunks of full sample X, if any
% 29. ts.chunks, subsample time series structure of chunks, if any
%
T=numel(X);
[~,~,lookvar] = mulbreaks(X);N=lookvar(1,2); if N>0;frac=N/T;else;frac=.5*T;end;
ts.frac=frac;
%
%  HHT, selecting optimal smoother by means of optimal siftings found via minimum PRMSE
er=eemd(X,0,1); % EMDs: first er==X, then highest to lowest frequency 
ser=size(er,2);for i=1:ser;SER(:,i)=sum((er(:,i:ser))');end
SEX=SER(:,3:end-1);sz=size(SEX,2); % Remove unwanted series
for k=2:sz;SS2=(SEX(:,k-1)-SEX(:,k)).^2;end;SD=sum(SS2)./(sum(SEX.^2));
b1=ERE(SD);
ts.HHT=SEX(:,b1); 
%
% SRT, simple speckoph
 [sop,som]= speckoph(X); ts.SRT=[sop,som];
%
%  Compare the 2 methods HHT and SRT: Prmes, error series, their stationarity and corrected stationarity test statistics, and crossovers of smoothed series wrt X
 XS=[ts.HHT,ts.SRT];
 XE=bsxfun(@minus,X,XS);
 ts.errs=XE;
 XE2=(bsxfun(@minus,X,XS)).^2; 
 ts.serrs=XE2; 
 ts.ersmev=[mean(ts.errs) var(ts.errs);mean(ts.serrs) var(ts.serrs)];
 %
 for i=1:3;y1=XE(:,i);[ls1,ns1] = linstat(y1,0,0,0);LSE1(i)=ls1(1);NSE1(i)=ns1(1);y2=XE2(:,i);[ls2,ns2] = linstat(y2,0,0,0);LSE2(i)=ls2(1);NSE2(i)=ns2(1);end; % ts.linerser=[NSE1, NSE2;LSE1,LSE2]; 
 % Harvey value and index number of linearity (1), nonlinearity (2) of errors and squared errors
 % for i=1:3;y=XE(:,i);[SST,CST] = statcorr_1(y,frac,1);CCC(i,:)=[SST,CST];end;ts.errstat=CCC;
 % for i=1:3;y=XE2(:,i);[SST,CST] = statcorr_1(y,frac,1);CCS(i,:)=[SST,CST];end;ts.serrstat=CCS;
 ts.prmse=sqrt(mean(XE2)/mean(X.^2));
% for i=1:3;y=XS(:,i);[~,ny] = crossovers(X,y);NY(i)=ny;end;ts.cross=NY;
 %
% Linear trend
trend=(1:T)';
a = polyfit(trend,X, 1);xhat=a(1)*trend+a(2); ts.ltrend=xhat; % Fitting straight line to x & obtaining linear trend of X
%
% Test for Normality  by Kolmogorov-Smirnov test, pvalue
[~,PL] = kstest(X); ts.normalp=PL;
%
% Harvey nonlinearity test and  ADF/KSS stationarity test 
[ls,ns] = linstat(X,0,0,0);ts.linstat1=ls;ts.linstat2=ns;
%
% Corrected Stationarity statistic
[~,CST]= statcorr_1(X,frac,0);ts.corrstat=CST;
%
% Linear trend, Smoother & other ratios, HHT & SRT
gox = gradient(xhat); ts.lratio=1-mean(gox);
% for i=1:3;y=XS(:,i);st(i)=y(T)/y(1);end;ts.sratio=st; % Smoother ratios
if ls(1)==2;[~,pratio] = ansaribradley(X(1:N),X(N+1:end));else;vs=[var(X(1:N)),var(X(N+1:end))];ft=max(vs)/min(vs);pratio=1-cdf('f',ft, T-N,N);end; ts.pratio=pratio; %  p-value of the Ansari-Bradley test
ts.cycleX=cyclex(X);ts.cratio=[cyclex(X(1:N))/cyclex(X(N+1:end))]; % Total cycle and subsample cycle ratio
% for i=1:3;cs(i)=cyclex(XS(:,i));end; ts.cycle=[ts.cycleX,cs];
%
% Conditional persistence (Kapeitanos, 2002), a check of speed of decay of autocorr. If series is I(1), I(0) the coefficient is expected to be small (fast decay) or large (slow decay)
[~,dec] = conper(X,1);ts.conper=dec;
%
% Computing turning points from phase of X and the three smoothers
px=phase(X);dpx=diff(px);nzx=nonzeros(dpx);TPX=numel(nzx); % TPX: turning points of X
% for k=1:3;ps(:,k)=phase(XS(:,k));dps(:,k)=diff(ps(:,k));nzs(k).a=nonzeros(dps);TPS(k,:)=numel(nzs(k).a);end; % TPS: turning points of the three smoothers
% ts.turnp=[TPX,TPS'];
%
% Comparing value and locational distance of X wrt HHT, SRT. Squared results put greater penalty on outlier distances
HS=[ts.HHT,ts.SRT];
% [MAVX,MALX]=max(X);[MIVX,MILX]=min(X);
for i=1:3;
y=HS(:,i);[a,b]=max(y);maxvals(i)=a;maxlocs(i)=b;[c,d]=min(y);minvals(i)=c;minlocs(i)=d;
end
% ts.mimas=[(minvals-MIVX).^2+(maxvals-MAVX).^2;(minlocs-MILX).^2+(maxlocs-MALX).^2];
%
%  Sequential cluster-based regime switches wrt entire sample 
[nreg,~,cldvs] = Regime_Switch(X,0,1);
if nreg>=1;cld=cldvs(:,1);ts.rsnum=numel(cld);ts.rsdates=sort(cld);else;ts.rsnum=0;ts.rsdates=0;end; 
%
%  Comparing full-sample with broken sample RS results, if any
[~,b]=min(ts.prmse);
ts.metsel=b; % Prefered smoother method (1 if HHT, 2,3 if SRT)
ts.xsels=[X,HS(:,b),ts.ltrend]; 
if isempty(ts.rsdates)==0;
   td=[ts.rsdates; numel(X)];  
   [~,tcs] = tchunks(td);
   ts.nums=nreg+1;ts.chunks=tcs; % Number and time strecth of X chunks
else
    ts.nums=0;ts.chunks=[];
end
%
if ts.nums>=2; for i=1:ts.nums;CX(i).a=X(ts.chunks(i).a);end;end; % Producing X chunks, where ts.nums >=2
%

    % ts=test_smooth(CX(i).a);SM(i).a=[ts.HHT,ts.SRT];PS(i,:)=ts.ersmev(2,4:end);PR(i,:)=ts.prmse;end; 
   % [~,b]=min(mean(PR));i=1:num;small=vertcat(SM(i).a);sels=small(:,b);ts.xselp=[X,sels];
   % [~,b]=min(mean(PS));i=1:num;small=vertcat(SM(i).a);selv=small(:,b);ts.xselv=[X,selv];
   % else
    % [~,b]=min(ts.ersmev(2,4:end));sels=XS(:,b);ts.xsels=[X,sels];end
%
%
%
end

