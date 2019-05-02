function lpt = large_peaks_troughs_2(X,dates,freq,sgop)
%
%  Code that produces largest peaks and troughs of time series X, stationary or nonstationary.
%  The series X is first filtered by two-step Savitsky-Golay/Daubechies Wavelet Transform (SGF/DWT) and then detrended. 
%  Peaks (troughs) are found in the positive (negative) region(s). Minimal distance is established within routine, no need to supply it from outside.
%  Utilizes Matlab my codes:  peaks_and_troughs.m, rempeaks.m, sgoldeb.m
%  GT,08.16.2013, revised 09.07.2013.
%
%  Function:  lpt = large_peaks_troughs_2(X,dates,freq,sgop);
%
%   INPUT:
%   X:(Tx1) original time series
%   dates:(Txn) dates in Excel format, daily, weekly, monthly, quarterly or  yearly, e.g.  '01-Dec-2010' or  '2010'
%   freq: time frequency of data vector i.e. daily, weekly, etc. Yearly=1, Monthly=2, Quarterly=3
%   sgop: option for SGF/DWT parameter (==savop in sgoldeb)
%
%  OUTPUT:
%  A structure containing 15 different results.
% lpt.peaX=[pval,pdate];
% lpt.troX=[tval,tdate];
% lpt.peax=[pval,pdate];
% lpt.trox=[tval,tdate];
% lpt.ltrend= long-run broken trend=DWT approximation
% lpt.dser=dser;
% lpt.cycle=cycles;
% lpt.phase=phase;
% lpt.xhat=xhat;
% lpt.betas=[res.betas(1),res.betas(2)];
% lpt.xnoise=WX; % MEMO: X=lpt.ltrend+lpt.xnoise
% lpt.omega=2*pi/lpt.cycle(1);  
% lpt.maxmin=[mav mad;miv mid;mavd madd;mivd,midd];
% lpt.shares: PCA shares & rank
% lpt.name: Daubechies combo (see 'sgoldeb.m')
%
if freq==1;
    ds=datenum(dates,'yyyy');
elseif freq==2
    ds=datenum(dates, 'mmmyyyy');
else
    ds=datenum(dates, 'QQ-YYYY');
end
% ds=datenum(dates,'yyyy');
T=size(X,1);trend=(1:T)';
%
%  Optimal Savitsky-Golay & DWT procedure
[~,x,WX,~,cycles,names] = Sgoldeb(X,freq,dates,0,sgop);
lpt.ltrend=x; % Long-run broken trend=DWT approximation
a = polyfit(trend,x, 1);xhat=a(1)*trend+a(2); % Fitting straight line to x
dser=x-xhat; % Obtained linearly detrended series
phase=unwrap(angle(dser)); %phase=unwrap(angle(gradient(dser))); % The angle function can be expressed as angle(z) = imag(log(z)) = atan2(imag(z),real(z)); phase is included in [0,pi]
%
% Maximum & minimum values and dates of original series X
[mav,ip]=max(X);[miv,it]=min(X);
mads=datevec(ds(ip));mids=datevec(ds(it));
if freq==1;mad=mads(1);mid=mids(1);else;mad=mads(1:3);mid=mids(1:3);end
% Maximum & minimum values and dates of dser
[mavd,ipd]=max(dser);[mivd,itd]=min(dser);
madsd=datevec(ds(ipd));midsd=datevec(ds(itd)); 
if freq==1;madd=madsd(1);midd=midsd(1);else;madd=madsd(1:3);midd=midsd(1:3);end
%
if cycles(2)/cycles(1)>1
    Pd=cycles(2);
else
    Pd=cycles(1);
end
%
% Peaks and troughs of X
[pe,tr]=peaks_and_troughs(X); 
pval=X(pe);ploc=pe;[pval,ploc] = rempeaks(pval,ploc,Pd);
tval=X(tr);tloc=tr;[tval,tloc] = rempeaks(tval,tloc,Pd);
dvp=datevec(ds(ploc));pdate=dvp(:,1);
dvt=datevec(ds(tloc));tdate=dvt(:,1);
lpt.peaX=[pval,pdate];lpt.troX=flipud([tval,tdate]);
%
% Peaks and troughs of dser+xhat
xd=dser+xhat;
[pe,tr]=peaks_and_troughs(xd); 
pval=xd(pe);ploc=pe;[pval,ploc] = rempeaks(pval,ploc,Pd);
tval=xd(tr);tloc=tr;[tval,tloc] = rempeaks(tval,tloc,Pd);
dvp=datevec(ds(ploc));pdate=dvp(:,1);
dvt=datevec(ds(tloc));tdate=dvt(:,1);
lpt.peax=[pval,pdate];lpt.trox=flipud([tval,tdate]);
%
% Results
%
lpt.dser=dser;
lpt.cycle=cycles;
lpt.phase=phase;
lpt.xhat=xhat;
lpt.betas=[a(1),a(2)];
lpt.xnoise=WX; % MEMO: X==x+xhat+noise
lpt.omega=2*pi/lpt.cycle(1);  
lpt.maxmin=[mav mad;miv mid;mavd madd;mivd,midd];
lpt.name=names;
% 
if phase(1)==0;dp1=0;else;dp1=pi;end;
dphase=[dp1;diff(lpt.phase)];
cneg=find(dphase<0,1,'last');cpos=find(dphase>0,1,'last');
cyc2=lpt.cycle(2)/2;
lpt.timelines=ceil([T, cyc2, cneg, cneg+cyc2, cpos, cpos+cyc2]);
%
% PCA of X: (1) cycle, (2) linear trend, (3)  noise
xx=[dser xhat X-dser-xhat];XX=xx'*xx/T;
[u,s]=svd(XX);
i=1:3; [~,rank]=max(abs(u(:,i))); % Rank
lpt.shares=[rank',diag(s)./trace(s)]; % Percent shares & rank of linearly detrended broken trend, linear trend & noise

%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS %%%%%%%%%%%%%%%%%%
     
 function [pe,tr]=peaks_and_troughs(x)
% Finds peaks and throughs of time series x, after Nagi Hatoum, Matlab File Exchange, copyright 2005
ds=diff(x);
ds=[ds(1);ds];%pad diff
filter=find(ds(2:end)==0)+1;%%find zeros
ds(filter)=ds(filter-1);%%replace zeros
ds=sign(ds);
ds=diff(ds);
tr=find(ds>0); % Troughs
pe=find(ds<0); % Peaks   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function [pval,ploc] = rempeaks(pval,ploc,Pd)
% Remove Peaks Separated By Less Than Min Peak Distance.
%  Start with the larger peaks to make sure we don't accidentally keep a small peak and remove a large peak in its neighborhood. 
%
%  INPUT: 
%  (1) pval: values of peaks available, obtained by findpeaks.m or  peaks_and_troughs.m
%  (2) ploc: location of peaks in (1)
%  (3) Pd: minimum affordable distance between peaks 
%
%  OUTPUT:
%  (1) pval: values of peaks found sorted in descending order
%  (2): ploc: location of peaks found
%
if isempty(pval) || Pd==1,
    return
end
% Order peaks from large to small
[pval, idx] = sort(pval,'descend');
ploc = ploc(idx);
idelete = ones(size(ploc))<0;
for i = 1:length(ploc),
    if ~idelete(i),
        % If the peak is not in the neighborhood of a larger peak, find secondary peaks to eliminate.
        idelete = idelete | (ploc>=ploc(i)-Pd)&(ploc<=ploc(i)+Pd); 
        idelete(i) = 0; % Keep current peak
    end
end
pval(idelete) = [];
ploc(idelete) = [];

