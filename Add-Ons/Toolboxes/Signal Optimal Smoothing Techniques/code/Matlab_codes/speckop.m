function [sop,som,lop,opn,res,ltrend,mes] = speckop(x)
%
%   Very fast method to find relative-variance minimizing SRT smoother and its associated lag. 
%   Produces also mean of SRT+envelope of x:(T,1). Uses my routine specrepop.m which calls specrep.m
%   function: [sop,som,lop,opn,res,ltrend,mes] = speckop(x,hop);
%
%   Input: signal x:(T,1)
%
%   Output:
%   sop: (T,1) time series of optimal smoother
%   som: (T,1) time series of mean of optimal smoother+envelope(x), usually with higher resolution than sop 
%   lop: optimal lag found
%   opn: procedure chosen among the six provided in specrepop.m
%   res: lookup table (2,2), 1st row has PRMSEs of sop & som, 2nd row has #  crossovers of  sop & som wrt original series x. 
%   Useful for selecting, if desired, output sop or som depending on PRMSE  and # of crossovers 
%   ltrend: linear trend of original series x
%   mes: mean of smoother (specrep) over all feasible lags 
%
T=size(x,1);trend=(1:T)';
a = polyfit(trend,x, 1);ltrend=a(1)*trend+a(2); 
[oplag1,oplag2,oplag3,oplag4,oplag5,oplag6,oplag7,mes]= specrepop(x,1);
ops= [oplag1,oplag2,oplag3,oplag4,oplag5,oplag6,oplag7];
for i=1:7;SX(:,i)=specrep(x,ops(i),1);end
xx=bsxfun(@minus,x,SX);
v1=mean(xx.^2);vxs=sqrt(v1./var(x)); % PRMSE
[~,opn]=min(vxs);
lop= ops(opn);
sop=SX(:,opn);
som=(sop+envelope(x,1))/2;
sos=[sop,som];
xss=bsxfun(@minus,x,sos);
[~,csop] = crossovers(x,sop);[~,csom] = crossovers(x,som); % Computing # of crossovers of fitted wrt original series
sox=sum(xss.^2);crs=[csop,csom];
res=[sox;crs];
 %
 %
 %
end

