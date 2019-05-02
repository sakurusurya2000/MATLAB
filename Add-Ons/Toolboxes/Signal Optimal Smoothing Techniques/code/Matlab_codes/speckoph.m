function [sop,som,lop,opn,res,ltrend,mes,SX] = speckoph(x)
%
%   Very fast method to find two-step Hamiltonian minimizing SRT smoother and its associated lag. 
%   Produces also mean of SRT+envelope of x:(T,1). Uses my routines specrep.m,envelope_mine.m, crossovers.m, hamilton.m 
%   function: [sop,som,lop,opn,res,ltrend,mes] = speckoph(x);
%
%   Input:
%   x: time series
%  
%   Output:
%   sop: (T,1) time series of optimal smoother
%   som: (T,1) time series of mean of optimal smoother+envelope(x), has  higher resolution than sop
%   lop: optimal lag found
%   opn: procedure chosen among the six provided in specrepop.m
%   res: lookup table (2,2), 1st row has PRMSEs of sop & som, 2nd row has #  crossovers of  sop, som wrt original series x. 
%   Useful for selecting, if desired, output sop or som depending on PRMSE  and # of crossovers 
%   ltrend: linear trend of original series x
%   mes: (T,1) time series of means of all lag-admissible smoothers
%
T=size(x,1);trend=(1:T)';
TT=[50,100,200,300,500,1000,1500,2500];
ST=ceil(sqrt(TT)); [~,bb]=min(abs(TT-T)); % Selecting maximum lags depending on size of T, excess lags may cause grid search process to freeze
SX=zeros(T,ST(bb));
for k=1:floor(ST(bb));sx=specrep(x,k,1); SX(:,k)=sx;[H0,H1,H10,HALL] = hamilton(x,sx,1);hs(k,:)=[H0,H1,H10,HALL];end
j=floor(ST(bb));mes=(mean(SX(:,1:j)'))';
[~,lop]=min(hs(:,3));sop=SX(:,lop);
xx=bsxfun(@minus,x,SX);
v1=mean(xx.^2);vxs=sqrt(v1./var(x)); % PRMSE
[~,opn]=min(vxs);
som=(sop+envelope_mine(x,1))/2;
sos=[sop,som];
xss=bsxfun(@minus,x,sos);
[~,csop] = crossovers(x,sop);[~,csom] = crossovers(x,som); % Computing # of crossovers of fitted wrt original series
sox=sum(xss.^2);crs=[csop,csom];
res=[sox;crs];
a = polyfit(trend,x, 1);ltrend=a(1)*trend+a(2); 
 %
 %
 %
end
%
function [mx,s1,s2]= envelope_mine(x,exop)
% Huang's envelope by cubic splining over the local extrema
% Has 'exop' option to extend x left- and rightwards (1) by 15% of T
% Function is  [mx,s1,s2]= envelope(x,exop)
%
% Inputs: x(T,1) signal, exop: option
% Outputs: mx: mean of the two extrema cubic splines, which are s1 (upper) and s2 (lower)
%
% Suggested use as a loop: for i=1:6;x= envelope(x,exop);XX(:,i)=x;end;
% Creates XX(T x 6) matrix of mean envelopes sequentially embedded, in fact, the EMDs produced by the code HHT_mine2.m
%
%   
T=numel(x);
if size(x,1)==1;x=x';end; % Make sure data vector is columnwise
if exop==1;Tex=ceil(.15*T);T=T+2*Tex;x=wextend('1D','sym',x,Tex);end;
[~,p1] = findpeaks(x);s1 = (spline([0; p1; T],[0; x(p1); 0],1:T))'; % Finding upper envelope
[~,p2] = findpeaks(-x);s2 = (spline([0; p2; T],[0; x(p2); 0],1:T))'; % Finding lower envelope
mx=(mean([s1,s2]'))';
if exop==1;mx=mx(Tex+1:end-Tex);s1=s1(Tex+1:end-Tex);s2=s2(Tex+1:end-Tex);end;
%
end
%
function [sx,cycles,phase,sm,ns,ss,omega,NSTEST,BB,BBS]= specrep(x,lags,exop)
% 
%   Code that produces Spectral Representation Theorem (SRT) function of signal x for select number of sincos lags. 
%   Deals with stationary/nonstationary and linear/nonlinear data. 
%   Works great as a smoother! Ref.: J.D. Hamilton, Time Series Analysis, PUP, 1994, Ch. 6
%   Optimal lags, if requested, is obtainable by means of SPECREPOP.M which iterates across lags this selfsame code
%
%   Function: = [sx,cycles,phase,sm,ns,ss,omega,NSTEST,BB,BBS]=specrep(x,lags,exop); Uses olsols.m & norser.m routines
%   Equation:  x(t)=SUM(A(j)*sin(omega*t-j))+SUM(D(j)*cos(omega*t-j)+constant;  j=1,...,lags
%                        
%   INPUTS:
%   (1) x:(T x 1) time series; 
%   (2) selected 1<= lags <=T/2
%   (3) exop: option to extend x left- and rightwards (1) by 15% of T, works bad with non normal series
%
%   OUTPUTS:
%   (1) sx:(T x 1) fitted time series
%   (2) cycles: cycles of raw and of smoothed series sx
%   (3) phase: phase angle of untrended sx
%   (4) sm: Portion of variance of x that can be attributed to cycles with frequency equal to lags
%   (5) ns: RMSE of sx / RMSE of x
%   (6) ss: sum of squared deviations between x and sx
%   (7) omega: frequency vector coresponding to given lags 
%   (8) NSTEST: (4,1) vector including info log: [index,names,ADF/KSS,pvalue];
%   (9) BB: (2*lags+1,1) vector of sine/cosine coefficients + constant. First (1 to lags) coeffs. are related to sine, 
%   subsequent (lags+1 to  2*lags) coeffs. are related to cosine, and last term is a constant 
%   (10) BBS: portion of the sample variance of x that can be attributed to cycles with frequency omega(j)
%   
if size(x,1)==1;x=x';end
T=size(x,1);trend=(1:T)';
if exop==1;Tex=ceil(.15*T);T1=T+2*Tex;x=wextend('1D','sym',x,Tex);end;
tnames={'ADF_test','KSS_test'};
%
%KSS symmetric ESTAR at 1% and 5%: UN(-2.79,-2.51), DM(-3.48, -2.93), DT(-3.93, -3.40)
outs1= nonlintest2(x,1);nts=cell2mat(outs1(2)); % Harvey nonlinearity test, nt~chisq(2);
star1= stars2(x,1); % KSS stationarity test in the presence of nonlinearity
%
% Testing for stationarity in linear/nonlinear series; % 5.99 is Chisq(2) at 95%
%  MEMO: p=1, series x is stationary, p=0 nonstationary
% If signal is linear, such that nts<5.99, use ADF test for stationarity, elseif signal is nonlinear, such that nts>5.99, use KSS test for stationarity 
if  nts<5.99; 
    [~,p,ADF]=adftest(x); % p-->0 x is stationary, p>0.05 x is nonstationary
    NSTEST=[num2cell(1),tnames(1),num2cell(ADF),num2cell(p)];
else
    KSS=cell2mat(star1(4,1)); % if abs(KSS)> critvalue, x is stationary ==> p=0
    if abs(KSS)>2.51;p=1;else;p=0;end; % -2.51 is UN5%
    NSTEST= [num2cell(2),tnames(2), num2cell(KSS),num2cell(p)];
end
%
j=1:lags;h=cell2mat(NSTEST(end));
if h==1;
     omega=(2*pi*j/T1)'; % Compute frequency for linear signal
else
     omega=(pi*j/T1)'; % Compute frequency for nonlinear signal
end
%
t=1:T1;SA=sin(omega*t)';SC=cos(omega*t)';
SAC=[SA,SC];
BB = real(olsols(x,SAC,1));
BBS=.5*(real(BB(1))^2+real(BB(2))^2); % Portion of the sample variance of x that can be attributed to cycles with frequency omega(j)
alpha=BB(1:lags);
delta=BB(lags+1:end-1);
SS=real([SAC,ones(T1,1)]);
sx=SS*BB;
sm=(.5*(sum(alpha.^2)+sum(delta.^2)))/var(x);
ns=rms(sx)/rms(x);
ss=norm(bsxfun(@minus,x,sx))/T;    
if exop==1;x=x(Tex+1:end-Tex);sx=sx(Tex+1:end-Tex);end;
a = polyfit(trend,sx,1);sxhat=a(1)*trend+a(2); % Fitting straight line to x
dser=sx-sxhat; % Obtained linearly detrended series
phase=unwrap(angle(dser)); % Phase of detrended series
%
cycle1= FFT_cycle(x);cycle2= FFT_cycle(sx);
cycles=[cycle1, cycle2];
end
%
function [cy,ny] = crossovers(x,y)
% Check to see the number and values of the crossovers of y wrt x
%  Abridged from Matlab Exchange File ID: #19642 by J. Millin
%  Function is: [cy,ny]=crossovers(x,y)
%
%  input: two signals x,y of equal length T
%  output: (1) cy: signal of length T1<T reporting exact values of crossovers; (2)  # of exact crossovers
%
[~,cy] = lineintersect(x,y);
ny=numel(cy);
end
%
function [H0,H1,H10,HALL,ES] = hamilton(x1,x2,hop)
%
% Finding indexes and the stable-branch root of the two-step LQ Hamiltonian
%
% The LQ Hamiltonian is: J=.5*H(x(t) s.t. s(t+1)=B(1)*s(t)+B(2)*x(t); where H is some quadratic index to be minimized in the actual variable x(t) 
% and in the estimated/unobservable variable s(t). The constraint, i.e. the dynamics of s(t), is s(t+1)=(.).
% Solution after finding Euler equations of J is in s(t) and in the costate lambda such that: z(t+1)*PHI=z(t)*PSI ==> z(t+1)=inv(PHI)*PSI*z(t), 
% where z(t+1)=[s(t+1),lambda(t+1)]; z(t)=[s(t),lambda(t)]
%
% The inputs respectively are x(t) and s(t), and hop=option to keep raw series (1) or detrend the series (2) or demean & detrend the series (3), trend is linear 
% The code returns:
% (1) the first-step index H0= the standard squared distance between x(t) and s(t); 
% (2) the second-step index H1=the squared distance between the variables' second moments; 
% (3) the summed two steps' index H10=H0+H1; 
% (4) HALL=H10+the squared fitted constraint, 
% (5) ES=the stable-branch root for solving z(t+1)
%
if size(x1,1)==1;x1=x1';end;if size(x2,1)==1;x2=x2';end; % Making sure the two series are entered as column vectors
T=numel(x1);trend=(1:T)';
if hop==2;
    a = polyfit(trend,x1, 1);ltrend=a(1)*trend+a(2);x1=x1-ltrend;b = polyfit(trend,x2, 1);ltrend=b(1)*trend+b(2);x2=x2-ltrend;
elseif hop==3;
    a = polyfit(trend,x1, 1);ltrend=a(1)*trend+a(2);x1=x1-mean(x1)-ltrend;b = polyfit(trend,x2, 1);ltrend=b(1)*trend+b(2);x2=x2-mean(x2)-ltrend;
end
%
xl=lagmatfill(x2,2,1,floor(T/4));
H0=mean((x1-x2).^2);
H1=mean((mean((x2-xl(:,1)).^2)-mean(((xl(:,1)-xl(:,2))).^2).^2)); 
H10=H0+H1;
[B,~,H3] = olsols(x2(2:end),[x2(1:end-1),x1(2:end)],1);
HALL=H10+H3;
ES=1/(2*(2*B(1) + B(2))); 
end
%