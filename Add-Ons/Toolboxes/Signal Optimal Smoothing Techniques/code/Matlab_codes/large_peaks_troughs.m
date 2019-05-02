function lpt = large_peaks_troughs(X,dates,freq,grop,trop,serop)
%
%  Code that produces largest peaks and troughs of time series X, stationary or nonstationary.
%  The series y is first filtered by two-step Savitsky-Golay/Daubechies Wavelet Transform and then detrended. 
%  Peaks (troughs) are found in the positive (negative) region(s). Minimal distance is established within routine, no need to supply it from outside.
%  Utilizes Matlab my codes:  peaks_and_troughs.m, rempeaks.m, remtroughs.m, Sgoldeb.m
%  GT,08.16.2013, revised 09.07.2013.
%
%  Function:  lpt = large_peaks_troughs(X,dates,freq, grop,trop,serop);
%
%   INPUT:
%   X:(Tx1) original time series
%   dates:(Txn) dates in Excel format, daily, weekly, monthly, quarterly or  yearly, e.g.  '01-Dec-2010' or  '2010'
%   freq: time frequency of data vector i.e. daily, weekly, etc. If yearly 1, elsewise 2
%   grop: option to produce graph (1) 
%   trop: option to include returns of X in graph
%   serop: option to use X raw series (1) or smoothed X series (2)
%
%  OUTPUT:
%  A structure containing 14 different results.
%  lpt.maxmin = [mav mad;miv mid;mavd madd;mivd;midd]; Maximum & minimum values and dates of original series (X) and of its smoothed transform  (dser)
%  lpt.peaks =  [pval, pdate];  pval & pdate: value and date of largest peaks of demeaned series, descending order of value magnitude
%  lpt.troughs = [tval,tdate]; tval & tdate: value and date of largest troughs of demeaned series, increasing order of value magnitude
%  lpt.plox, lpt.tlox = Location of peaks & troughs, respectively
%  lpt.cycle = mean of mean peak-to-peak and trough-to-trough cycle durations (if any) of demeaned series
%  lpt.sgolx:(Tx1)= SG-filtered series
%  lpt.xnoise:(Tx1) difference between original series and SG-filtered series
%  lpt.wnoises:((Tx2)=[wnoisex,wnoiseX]; noise wavelet detailed decomposition noise
%  lpt.xhat:(Tx1) trend line of SGF/DWT-filtered original series
%  lpt.dser:(Tx1) SG-filtered and demeaned series
%  lpt.phase:(Tx1) series of phase of lpt.dser [0,pi]; Is nonzero only with serop=1
%  lpt.omega: real expressing average frequency of demeaned series
%  lpt.timelines: (1x6) vector including (1) timelength of X, (2) halflife average cycle length, (3) cycle length of normalized phase<1, (4) (2)+(3),  (5) cycle length of normalized phase>1, (6) (2)+(5)
%  lpt.betas: constant and trend slope terms in OLS that determines xhat, and related tstatistics
%  lpt.sgstats: optimal polynomial degree & window size of SG filtered series x
%
% MEMO: w/serop=1, X=lpt.dser; w/serop=2, X=lpt.dser+lpt.xhat
%
if freq==1;ds=datenum(dates,'yyyy');else;ds=datenum(dates);end
T=size(X,1);trend=(1:T)';
savop=2;
%
% Optionally use returns
if trop==1;rets=[price2ret(X);mean(price2ret(X(1:24)))];end
%  Optimal Savitsky-Golay & DWT procedure
[~,WA,WX,~,cycles] = Sgoldeb(X,freq,dates,grop,savop);
%[SX,pars] = sgolaut(X,freq,dates,pmax,dmax,grop,exop);% yf=abs(fft(SX));[~,ii]=max(yf(2:end)); cycle=ceil(T/ii);% wpol=pars(1);dpol=pars(2);x=X-SX; % wpol=pars(1);dpol=pars(2);
if serop==1;x=X;else;x=WA;end; % Select among raw or smoothed X
res=ols_mine(x,trend,.05);xhat=res.betas(2)+res.betas(1)*trend; % Fitting linear trend to x
%
if serop==1;dser=x;else; dser=x-xhat;end;
%
% Maximum & minimum values and dates of original series X
[mav,ip]=max(X);[miv,it]=min(X);
mads=datevec(ds(ip));mids=datevec(ds(it));
if freq==1;mad=mads(1);mid=mids(1);else;mad=mads(1:3);mid=mids(1:3);end
% Maximum & minimum values and dates of smoothed X
[mavd,ipd]=max(dser);[mivd,itd]=min(dser);
madsd=datevec(ds(ipd));midsd=datevec(ds(itd)); 
if freq==1;madd=madsd(1);midd=midsd(1);else;madd=madsd(1:3);midd=midsd(1:3);end
%
phase=unwrap(angle(dser)); % Find phase [0,pi]
%
[pe,tr]=peaks_and_troughs(dser); % Peaks and troughs of function (dser)
%
% Retain in dser only largest peaks and troughs. There must be a minimum of 2 each.
%
if serop==1
pes=[dser(pe),pe];
trs=[dser(tr),tr];
else
cp=find((dser(pe))>0);pes=[dser(pe(cp)), pe(cp)];
cn=find((dser(tr))<0);trs=[ dser(tr(cn)), tr(cn)];
end
%
pp=[1;pes(:,2)];
Pdp=max(abs(diff(pp))); 
pval=pes(:,1);ploc=pes(:,2);
[pval,ploc] = rempeaks(pval,ploc,Pdp);
dvp=datevec(ds(ploc));
pdate=dvp(:,1:3);
%
tt=[1;trs(:,2)];
Pdt=max(abs(diff(tt))); %Pdt=ceil(mean(diff(tt)));
tval=trs(:,1);tloc=trs(:,2);
[tval,tloc] = remtroughs(tval,tloc,Pdt);
dvt=datevec(ds(tloc));
tdate=dvt(:,1:3);
%
[~,Loc]= ismember(trend,ploc);plox=(Loc>0)'; % Location of peaks of function (x)
[~,Loc]= ismember(trend,tloc);tlox=(Loc>0)'; % Location of troughs of function (x)
% [~,~,~,wnoiseX,cycX] = Sgoldeb(X,freq,dates,0,2); % Optimal Daubechies transform of X;[~,~,~,wnoisex,cycx] = Sgoldeb(x,freq,dates,0,2); % Optimal Daubechies transform of x
%
% Results
%
lpt.xhat=xhat;
lpt.maxmin=[mav mad;miv mid;mavd madd;mivd,midd];
lpt.peaks=[pval, pdate];
lpt.troughs=[tval,tdate];  
lpt.plox=plox;lpt.tlox=tlox;
lpt.cycle=cycles;
lpt.xnoise=WX; % MEMO: X==x+xhat+noise
lpt.dser=dser;
lpt.phase=phase;
lpt.omega=2*pi/lpt.cycle(1);  
% 
if phase(1)==0;dp1=0;else;dp1=pi;end;
dphase=[dp1;diff(lpt.phase)];
cneg=find(dphase<0,1,'last');cpos=find(dphase>0,1,'last');
cyc2=lpt.cycle(2)/2; % Half-life long-term cycle
lpt.timelines=ceil([T, cyc2, cneg, cneg+cyc2, cpos, cpos+cyc2]);
% Prepare the figure, if grop=1
%
if grop==1
 M=max(dser);
 mindate=min(ds); maxdate=max(ds);
 hf=figure('Units','normalized','Position',[0.05 0.05 0.9 0.7]);
 ha=get(hf,'CurrentAxes');
 set(ha,'XTickLabelMode','manual','XTickMode','manual','XLimMode','manual',...
    'XLim',[mindate maxdate],'XTick',linspace(mindate,maxdate,10));
%
  subplot(2,2,1)
  plot(ds,X,'-k',ds,dser,'--k',ds,zz,'-k','MarkerSize',12), grid on
  xlabel('Time');ylabel('Series');
  ht=title('Pane 1') ;set(ht,'Fontsize', 10);
  legend('Original series','Detrended series','Location','NorthWest')
  datetick('x','yyyy')
% 
 subplot(2,2,2)
  plot(ds,X,'-k',ds,wnoiseX,'--k',ds,zz,'-k','MarkerSize',12), grid on
  xlabel('Time');ylabel('Series');
  ht=title('Pane 2') ;set(ht,'Fontsize', 10);
  legend('Original series','Noise','Location','NorthWest')
  datetick('x','yyyy')
 %
  subplot(2,2,3)
  plot(ds,dser, '-.k',ds,plox*M,'-k',ds,-tlox*M,'-k',ds,haar*M,'-k',ds, zz,'-k','MarkerSize',12), grid on
  xlabel('Time');ylabel('Series');
  ht=title('Pane 3') ;set(ht,'Fontsize', 10);
  legend('Demeaned series','Daubechies wavelet coefficients','Peaks (+),Troughs (-)','Location','NorthWest')
  datetick('x','yyyy')
%
  subplot(2,2,4)
  if trop==1
     plot(ds,rets*10,'-.k',ds,phase,'-k',ds,zz,'-k','MarkerSize',12), grid on
     xlabel('Time');ylabel('Series');
     ht=title('Pane 4') ;set(ht,'Fontsize', 10);
     legend('Returns','Phase of demeaned series','Location','SouthWest')
     datetick('x','yyyy')
  else
     plot(ds,dser,'-.k',ds,phase*M,'-k',ds,zz,'-k','MarkerSize',12), grid on
     xlabel('Time');ylabel('Series');
     ht=title('Pane 4') ;set(ht,'Fontsize', 10);
     legend('Demeaned series','Phase of demeaned series','Location','SouthWest')
     datetick('x','yyyy')
  end
end


   