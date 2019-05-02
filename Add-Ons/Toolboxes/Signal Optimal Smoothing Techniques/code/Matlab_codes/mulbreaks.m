function [tbdates,tbvals,lookvar] = mulbreaks(x)
%
%   Multiple break dates of linear or nonlinear signal x based on finding local minima of squared residuals of a dynamic regression with unknown
%   time breaks in constant and trend, and trimming. Break dates and values are exhibited in order of appearance. A lookup table of descending 
%   variance ratios together with their break dates  is also produced, where variance ratio is Var(X(tbdate))/Var(X)
%   Code allows no more than five time breaks.
%   Ref.: Bai J. and Perron P. (2003) Computation and Analysis of Multiple Structural Change Models, Journal of Applied Econometrics, 18, 1-22.
%   Guido Travaglini, June.05.14
%   Function is: [tbdates,tbvals,lookvar] = mulbreaks(x);
%
T=numel(x);trend=(1:T-1)';
%
% Harvey nonlinearity test, H0 is nts~chisq(2);
outs1= nonlintest2(x,1);nts=cell2mat(outs1(2)); 
if nts>5.99; NT=2;else;NT=1;end
%
% Preparing input series for dynamic regression w/breaks
x21=x.^2;x2=x21(1:end-1);
x31=x.^3;x3=x31(1:end-1);
maxlag=ceil(12*(ceil(T/100))^.25); 
bct=triu(ones(T,T)); % Matrix of breaks in constant
bcmat=bct(2:T,1:T-1); % Matrix of breaks in trend
bmat=zeros(T,T);
for j=1:T;for i=1:T;bmat(j,i)=(i-j+1);end;end
btmat=triu(bmat(2:T,1:T-1));
dx=diff(x);
dxfillags=lagmatfill(dx,maxlag,1,ceil(T/10));
%
for k=1:maxlag;
    if NT==2;
      R=[x3,dxfillags(:,1:k)];
    else
      R=[x(1:end-1),dxfillags(:,1:k)];
    end
    [B,~,~,r]=olsols(dx,R,1);
    sig2=sum(r.^2)/(T-maxlag);
    LS=log(sig2);
    ttk=((B(1)^2)/sig2)*sum(x2(maxlag+1:end));
    MC=2*(ttk+k)/(T-maxlag);
    LOGL(k)=LS+MC;
end
%
[~,maiclag]=min(LOGL);
trim=floor(.1*T);
%
% Produce dynamic regression of centered and scaled variables w/breaks
for i=trim:T-trim
  if NT==2
         [~,TS1,S21] = olsols(zscore(dx),zscore([x3,dxfillags(:,1:maiclag),trend,bcmat(:,i),btmat(:,i)]),1);
  else
         [~,TS1,S21] = olsols(zscore(dx),zscore([x(1:end-1),dxfillags(:,1:maiclag),trend,bcmat(:,i),btmat(:,i)]),1);
  end
  if isreal(TS1)==0;TS1=imag(TS1);end;
  TSB1(i)=TS1(3); SS1(i-trim+1)=S21;
end
if size(SS1,1)==1;SS1=SS1';end
%
[~,tlx]=peaks_and_troughs(SS1);
Pd= FFT_cycle(SS1);
if Pd/T>.1;Pd=ceil(Pd/2);end
tks=SS1(tlx);
[~,tlx] = remtroughs(tks,tlx,Pd);
tbdates= sort(tlx); tbvals=SS1(tbdates);
if size(tbdates,1)<=5;
    tbdates=tbdates+trim;tbvals=x(tbdates);
else
   svals=sort(tbvals,'ascend');
   for i=1:5;a=find(SS1==svals(i));tbdates5(i)=a;end
   tbdates=tbdates5+trim;tbvals=x(tbdates);
end
if isempty(tbdates)==0;
 if size(tbdates,1)>1;tbdates=tbdates';end;if size(tbvals,1)>1;tbvals=tbvals';end
else
    tbdates=[];tbvals=[];
end
%
% Finding variances corresponding to reported break dates
VX=(bsxfun(@minus,x,mean(x))).^2; 
las= [VX(tbdates)/var(x),tbdates'];
if  isempty(tbdates)==0;
sal=sortrows(las,1);
lookvar= flipud(sal);
else
    lookvar=[];
end
%
%
%
%
end

