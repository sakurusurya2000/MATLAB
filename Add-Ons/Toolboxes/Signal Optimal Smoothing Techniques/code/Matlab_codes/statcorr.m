function [CST,CSI]= statcorr(SST,VEC)
%
%  Correction for ADF/KSS stationarity test. 
%
%  Inputs: 
%  SST: Original ADF/KSS stationarity test and VEC=[SLR,SSR,SVR,SCR], where:
%  1. SLR: Subsample linear trend ratio; Under null of non-growing linear trend: E(SLR)=1
%  2. SSR: Subsample smoother ratio; Under null of non-growing smoother: E(SSR)=1;
%  3. SVR: p-value of subsample variance ratio; Under null of equality of variances: E(SVR)=1
%  4. SCR: Subsample cycle length ratio; Under null of subsample equal cycle length: E(SCR)=1
%  
%  Output:
%  CST: corrected stationarity test statistic
%  CSI: corrected stationarity index CST/SST
%
CST= -[abs(SST)-sum(abs(abs(VEC)-1))];
CSI=abs(CST)/abs(SST);
%
%
%
end