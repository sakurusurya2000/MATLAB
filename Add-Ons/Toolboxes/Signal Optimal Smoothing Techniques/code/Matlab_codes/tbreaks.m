function [bc,bt,tbc,tbt,otb,ott,bands] = tbreaks(y,x,z,trimin,trimen,graph,equ)
%
% Finds time series of coefficients and t-statistics of constant & trend terms in regression of y=f(x) with initial trimming.
% Note: Trimming is asymmetric here: initial is customary, final is minimal e.g. 2 or 3 observations. Requires OLS_MINE.M and GMM_MINE.M.
% 
% Function is: [bc,bt,tbc,tbt,otb,ott,bands] = tbreaks(y,x,z,trimin,trimen,graph,equ)
% 
% INPUT: 
% y=endogenous variable time series, Tx1.
% x=exogenous variable(s) time series, Txk, k>=1.
% z=instruments time series, Txq, q>=k.
% WARNING: if only y is available, x,z may be set as lags/leads of y.
% trimin  = percent initial trimming, usually set to: .10-.15
% trimen = integer, few will suffice
% graph=1 option to plot output, 0 noplot
% Option to choose equation method: OLS is EQU=0; GMM is EQU=1.
%
% OUTPUT (interval is Tt=trimin*T:T+trimen):
% bc, bt = time series of coefficients of breaks of constant & trend, Ttx1 each
% tbc, tbt = time series of t-statistics of breaks of constant & trend, Ttx1 each
% otb = other coefficients: (Tt x k+2); ott = other coefficient tstats : (Tt x k+2).
% bands= optimal bandwidth for each tbreak when GMM is used. Else, bands=0.
% graph=1 plots tbc & tbt; graph =0 no plot.
%
%
T=length(y);trend=1:T;LBEG=round(trimin*T);LEND=round(T-trimen*T); %trend must be Tx1 vector: [1,2,...T]'
c=ones(length(y),1);alpha=.05;
betaser=zeros(LEND-LBEG+1,size(x,2)+4);tstatser=zeros(LEND-LBEG+1,size(x,2)+4);
%
% Creating TxT matrices of bc & bt.
%
bcmat=tril(ones(T,T),-1);
ra=repmat(trend,T,1)';
rap=zeros(T,T);
for i=1:T; rap(:,i)=ra(:,i)-i;end;
btmat=tril(rap);
%
% Looping over the trimmed interval.
%
if equ==0   % get OLS results:
 for i=LBEG:LEND
   tst=[bcmat(:,i),btmat(:,i),trend'];
   j=i-LBEG+1;
   res=ols_mine(y,[tst,x],alpha); % Constant term is automatically included in invoked procedure.
   betaser(j,:)=res.betas(:,2)'; tstatser(j,:)=res.tstats(:,2)';bands(j)=NaN;
   bc(j)=betaser(j,1);bt(j)=betaser(j,2);tbc(j)=tstatser(j,1);tbt(j)=tstatser(j,2);
   otb(j,:)=betaser(j,3:end);ott(j,:)=tstatser(j,3:end);
  end
 %
else   % if equ==1, get GMMH results:
 for i=LBEG:LEND
   tst=[bcmat(:,i),btmat(:,i),trend'];
   j=i-LBEG+1;
   res= gmm_mine(y,[tst,x],[tst,c,z]);  % Make sure that both x and z include tst and the constant term.
   betaser(j,:)=res.betas(:,2)'; tstatser(j,:)=res.tstats(:,2)';bands(j)=res.hac;
   bc(j)=betaser(j,1);bt(j)=betaser(j,2);tbc(j)=tstatser(j,1);tbt(j)=tstatser(j,2);
   otb(j,:)=betaser(j,3:end);ott(j,:)=tstatser(j,3:end);
 end
% end
%
bc=bc';bt=bt';tbc=tbc';tbt=tbt';otb=otb';ott=ott';
if graph==1
 plot(1:length(tbc'),tbc','*',1:length(tbc'),tbt','+')
 title('Time series of t-statistics of breaks in constant and trend','FontWeight','bold')
 legend('t-statistic of break in constant','t-statistic of break in trend',2);
end
%
%
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ADD-IN FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%

function resols = ols_mine(y,x,alpha)
% OLS estimator. Do not add constant: it is already included at end of regressor list. 
% Function is:  resols = ols_mine(y,x,alpha)
%
%  INPUT: 
%  y=  LHS endogenous variable, Tx1
%  x= matrix of RHS exogenous variable(s) : Txk, k=>1 
%  alpha=confidence region [usually .05]
%
% OUTPUT:
%  resols, a structure containing:
%  1) beta_OLS, 2) tstats_OLS ,3) bint_OLS, respectively the coefficients, the t-stats and their CIs.
%  4) v_OLS: the variance-covariance matrix, 5) r_OLS and 6) rint_OLS: the residual and their CIs.
%  7) gstats_OLS includes general statistics:  Rsquared (r2_OLS); Fstat and pval (F, probF) 
%   Error variance (s2); Sum of squared errors(SSE); LogLik(Log-likelihood); Durbin-Watson coefficient of error autocorrelation (DW_OLS),
%  8) se_OLS: the coefficients'  standard errors.
%  9) fitted values
%
%
isc= iscolumn(y); if isc==0;y=y';end;                                                                      % Determine whether the input vector y is column. If not, it is automatically redressed 
T=length(y);if T~=length(x);x=x';end;                                                                      % Make sure vector/matrix x has same length as y
c=ones(T,1);x=[x,c];                                                                                              % Including constant in regressor list
beta_OLS=pinv(x'*x)*(x'*y);                                                                                  % OLS parameter vector
px=size(beta_OLS,1);nu=T-px;                                                                               % Residual degrees of freedom
%
% Find confidence intervals for each component of beta_OLS           
 r_OLS = y-x*beta_OLS;                                                                                         % Residuals.
 normr = norm(r_OLS);                                                                                            % Sqrt(sum(r.^2))
 rmse = normr/sqrt(nu);s2 = rmse^2;                                                                         % Estimator of RMSE and error variance. 
%
% Find the standard errors of the residuals.
%
ser=sqrt(sum(r_OLS.^2)/nu);tval = tinv((1-alpha/2),nu);                                             % tinv: inverse t distribution
rint_OLS = [(r_OLS-tval*ser) (r_OLS+tval*ser)];
v_OLS=ser^2*pinv(x'*x);                                                                                         %   TxT Covariance matrix of OLS (e'e)
se_OLS=sqrt(diag(v_OLS));                                                 
tstats_OLS=beta_OLS./se_OLS;  
bint_OLS = [beta_OLS-tval*se_OLS, beta_OLS+tval*se_OLS];
%
% Calculate R-squared and the other statistics.
%
SSE = normr.^2;                                                                                                         % Error sum of squares: r_OLS'*r_OLS
RSS = norm(y-x*beta_OLS-mean(y))^2;                                                                     % Regression sum of squares.
TSS = norm(y-mean(y))^2;                                                                                         % Total sum of squares.
r2_OLS = 1 - SSE/TSS;                                                                                             % R-square statistic.
F = (RSS/(px-1))/s2;                                                                                                   % F statistic for joint null of regression
probF = 1 - fcdf(F,px-1,nu);                                                                                       % Significance probability for regression
LogLik = -.5*T*(1+log(2*pi)+log(SSE/T));
DW_OLS=sum((r_OLS(2:T)-r_OLS(1:T-1)).^2)/sum(r_OLS.^2);
gstats_OLS = [r2_OLS F probF s2 SSE LogLik DW_OLS];
resols.betas= beta_OLS;resols.tstats=tstats_OLS; resols.bint=bint_OLS; resols.vols=v_OLS;
resols.rols=r_OLS;resols.rint=rint_OLS; resols.gstats=gstats_OLS;resols.ses=se_OLS;
% resols.fit=x*(resols.betas(1:end-1))+resols.betas(end);
%
%
%
end


function res= gmm_mine(y,x,z)
% Produces results of TSLS and ordinary 2-step, HAC and Jacknife 2-step GMM estimation of  multiple linear 
% regression equation with (possibly) endogenous and/or mismeasured regressors and selected instruments.
% WARNING: CONSTANT TERM IS INCLUDED AS LAST REGRESSOR!
% Function:  res= gmm_mine(y,x,z)
%   
% INPUTS:
% y=endogenous variable (Tx1)
% x=regressors (Txk)
% z=instruments (Txq), q=>k
%
% OUTPUT:
% res, a structure, includes:
%
% res.hac: # of HAC lags
% res.mxe: Mean covariance between x and structural equation residuals
% res.granger: Granger causality F stats running from residuals to x and viceversa + critical value (c_v)
% res.bze: : Fstat and pvalue of coefficient vector of relationship beween structural OLS equation residuals and  instruments z
% res.RF: R^2, Reduced-Form F statistics and their pvals; 
% res.PIE: RF coefficient matrix; 
% res.erCD: Excess rank of PIE and Cragg-Donald statistic;
% res.inst: Mineig of instrument orthog. and GPg.
% res.betas, res.tstats: betas and tstats
% res.genstats: standard error of estimation, J-stat and its p-value, Asymptotic J-stat and its p-value, Anderson-Rubin statistic and its pvalue,
%  Asymptotic Anderson-Rubin statistic and its pvalue, Durbin-Watson statistic, ARCH of residuals-1 lag,  First-order Autocorrelation coeff.
% res.haus: standard Hausman test & pvalue
% res.aicbic: Akaike and Schwartz ICs on # of regressors (k)
% res.bias: Different forms of bias (Newey & Smith 2004)
% res.foxes: FOCS
% res.resids: time series of estimated residuals for each of the three methods
%
%
n=length(y);c=ones(n,1); x=[x,c];        % Puts constant term at end of regressor list
k=length(x(1,:));q=length(z(1,:));
xx=x'*x; xy=x'*y;zzi = pinv(z'*z); zx = z'*x; zy = z'*y;P=z*zzi*z';M=eye(n)-P;
beta_OLS=pinv(xx)*xy;r_OLS=y-x*beta_OLS; sigmasq_OLS=sum(r_OLS.^2)/n;v_OLS=sigmasq_OLS*pinv(xx);
% A. Find relationship between regressors and OLS disturbance i.e. 'endogeneity'. 2 ways: x'e and Granger causality.
xe=bsxfun(@times,x,r_OLS);meanxe=mean(sqrt(xe(2:n,:).^2));
max_lag=6; alpha=.10; % Granger causality test
for i=1:k-1
  [Fe,c_v] =  granger_cause(x(:,i),r_OLS,alpha,max_lag);  % Past errors causes current x
  [Fx,c_v] =  granger_cause(r_OLS,x(:,i),alpha,max_lag);  % Past x causes current errors
ecx(i)=Fe;xce(i)=Fx;
end
% B. Find relationship between instruments and OLS disturbance i.e. 'exogeneity':
% [~,~,~,~,stats]=regress(r_OLS,[z,c]); 
% bz_fstats=stats(2);bz_fpvals=stats(3);
alpha=.05;resols=ols_mine(r_OLS,z,alpha); 
bz_fstats=resols.gstats(2);bz_fpvals=resols.gstats(3);
%for i=1:q
 %   [b,bint,r,rint,stats]=regress(r_OLS,[z(:,i),c]);bze(i)=b(1); 
  %  bz_rhosq(i)=stats(1);bz_fstats(i)=stats(2);bz_fpvals(i)=stats(3);
% end
%
%  C. Find relationship between instruments and regressors, i.e. 'relevance'. Three ways to obtain it: RF matrix, Concentration Parameter and Cragg-Donald statistic
%  C1. Compute RF (PIE) of x wrt all z, and compute V matrix, first-stage R^2, F statistics and their pvals
%for i=1:k-1
%[b,~,r,~,stats]= regress(x(:,i),[z,c]);BRF(:,i)=b;v(:,i)=r;
%BRF(i,:)=b;v(i,:)=r; % xhat(:,i)=x(:,i)-[z,c]*b;
%fs_rhosq(i)=stats(1);fs_fstats(i)=stats(2);fs_fpvals(i)=stats(3);
%end
for i=1:k-1
 resols=ols_mine(x(:,i),z,alpha);  
 BRF(:,i)=resols.betas;v(:,i)=resols.rols;
fs_rhosq(i)=resols.gstats(1);fs_fstats(i)=resols.gstats(2);fs_fpvals(i)=resols.gstats(3);
end
PIE=BRF(1:end-1,:)';
ErankPIE=sprank(PIE)-size(PIE,1); % PIE:(k-1 x q-1) is coefficient matrix that connects x to z, it must be of full rank to ensure relevance. ErankPIE computes excess rank (0 if full).
[~,SP]=svd(PIE);minsvPIE=min(diag(SP));           % Find minimum singular value of Pie
% C2. Compute minimum eigenvalue of concentration parameter (Stock & Yogo, 2004)
sigma=v'*v/(n-k-q);CP=(PIE*(z'*z)*PIE')/(sigma+eps*eye(size(sigma,1)));CP=min(eig(CP));
%sigma=v*v'/(n-k-q);zd = bsxfun(@minus, z, mean(z));CP=(PIE*(zd'*zd)*PIE')/sigma;CP=min(eig(CP)); %use of zd suggested by Imbens, NBER, 2007
% C3. Compute Cragg-Donald statistic for instrument weakness, excludes constant
fs_sigma=v'*v/((n-q-k+2)^2);gt=(fs_sigma^-.5*x(:,1:end-1)'*P*x(:,1:end-1)*fs_sigma^-.5)/q;CD=min(eig(gt));
% fs_sigma=v*v'/((n-q-k+2)^2);gt=(fs_sigma^-.5*x(:,1:end-1)'*P*x(:,1:end-1)*fs_sigma^-.5)/q;CD=min(eig(gt));%
beta_TSLS=pinv(zx'*zzi*zx)*(zx'*zzi*zy);r_TSLS=y-x*beta_TSLS; sigmasq_TSLS=sum(r_TSLS.^2)/n;v_TSLS=sigmasq_TSLS*pinv(zx'*zzi*zx);
Haus_TSLS=(beta_OLS-beta_TSLS)'*pinv(v_TSLS-v_OLS)*(beta_OLS-beta_TSLS);
% Producing moments  u: (Txq), E(0,sigmasq).... 
u=bsxfun(@times,z,r_TSLS);  meanu=mean(u); 
% .... & jacobians G: (qxk), and mineig of G'*G (k x k) and min singular value of G
G=-zx/n;mineigG=min(eig(G'*G)); [~,SJ]=svd(G');minsvdG=min(diag(SJ));
% Second step (Ordinary GMM)
%
w=pinv(u'*u);v_GMM=pinv(zx'*w*zx); beta_GMM=v_GMM*(zx'*w*zy);
r_GMM=y-x*beta_GMM;se_GMM=sqrt(sum(r_GMM.^2)/n);
Haus_GMM=(beta_OLS-beta_GMM)'*pinv(v_GMM-v_OLS)*(beta_OLS-beta_GMM);pval_Haus_GMM=1-cdf('chi2',Haus_GMM,k);
AIC_GMM=log(sum(r_GMM.^2/n))+2*k/n;BIC_GMM=log(sum(r_GMM.^2/n))+log(n)*k/n;
DW_GMM = sum(diff(r_GMM).^2)/sum(r_GMM.^2); [~,ARCH_GMM_pval,ARCH_GMM,critval] = archtest(r_GMM,1,.05);  %Engle ARCH test of residuals, 1 lag
CR=corrcoef(r_GMM(2:end),r_GMM(1:end-1));AU_GMM=CR(1,2); % First-order autocorrelation coefficient
betaseg=sqrt(diag(v_GMM));tstats_GMM=beta_GMM./betaseg;
%j_GMM=n*mean(u)*w*mean(u)'; 
sigmasq_GMM=sum(r_GMM.^2)/n;j_GMM=r_GMM'*P*r_GMM/sigmasq_GMM;j_GMM_ASY=(j_GMM-q)/sqrt(2*q); % Standard and asymptotic Jstat (Andrews & Stock, 2007)
pj_GMM=1.-cdf('chi2',j_GMM, q-k);pj_GMM_ASY=1-cdf('Normal',j_GMM_ASY,0,1);  
AR_GMM=(n-q)*(r_GMM'*P*r_GMM)/(r_GMM'*M*r_GMM); p_AR_GMM=1.-chi2cdf(AR_GMM, q);  % Standard and asymptotic AR (Anderson-Rubin) test for betas=beta_GMM.
AR_GMM_ASY=sqrt(q)*((AR_GMM/q)-1.);pAR_GMM_ASY=1-cdf('Normal',AR_GMM_ASY,0,2);
fox_GMM=G'*w*meanu';
% Computing H,P,a and biases BI, BO, BG (Newey & Smith 2004). Uses their notation.
omega=u'*u/n;sigma=pinv(G'*pinv(omega)*G);CP_GMM=min(eig(G'*pinv(omega)*G)); %zx'*w*zx=n*G'*pinv(omega)*G; sigma(Long Run variance)=n*v_GMM; n*w=inv(omega)
H_GMM=sigma*G'*pinv(omega);P_GMM=pinv(omega)-pinv(omega)*G*H_GMM;GPg_GMM=meanu*P_GMM*meanu';
avec=zeros(q,1);xzj=zeros(k,q);for j=1:q;xzj(:,j)=-(x/n)'*z(:,j);avec(j)=trace(sigma*xzj(:,j)*xzj(:,j)')/(2*n);end
BIG=H_GMM*(-avec+G*H_GMM*meanu')/n;BOG=H_GMM*(meanu*meanu'*P_GMM*meanu')/n;BGG=-sigma*(((q-k)*meanxe/(se_GMM^2)))'/n;
betastats_gmm=[beta_GMM,tstats_GMM];
stats_gmm=[se_GMM,j_GMM,pj_GMM,j_GMM_ASY,pj_GMM_ASY,AR_GMM,p_AR_GMM,AR_GMM_ASY,pAR_GMM_ASY,DW_GMM,ARCH_GMM,ARCH_GMM_pval,AU_GMM];
bias_gmm=[BIG,BOG,BGG];
%
% GMM with HAC, Newey-West (NW) Bartlett window (uses two codes by Kyriacou's GMM-TOOLBOX: optbandw-Bartlett and kernelest).
%
%
% Perform prewhitening of moments (u) by means of VAR(1) [uses LeSage Spatial Econometrix var_resid.m code]
%  resid = var_resid(u,1);ut=u(2:end,:)-resid;
b = u*ones(q,1); nyu = 1; nn = fix((n^(1/9))^2);cgamma = 1.1447;c = toeplitz(b,b(1:nn+1));
transfb = repmat(b,1,nn+1);
sigma = (1/n)*(sum(c.*transfb))';                                  % NW bandwidth selection, step 2
s0 = sigma(1,1)+2*sum(sigma(2:end,1));       
jj = (1:nn)';                                                                    % NW bandwidth selection, step 3 
snu = 2*sum((jj.^nyu).*sigma(2:end,1));  
gammahat = cgamma*((snu/s0)^2)^(1/(2*nyu+1));           % NW bandwidth selection, step 4
bw = fix(gammahat*(n^(1/(2*nyu+1))));
W2 = zeros(n-bw-1,1); a = (1:bw)/(bw+1);
% if bw==0;D=eye(n);end;
% if bw~0; W1 = [1;(1-a)'];W = [W1;W2]; D = (toeplitz(W))';end;
if bw==0;D=eye(n);
else
    W1 = [1;(1-a)'];W = [W1;W2]; D = (toeplitz(W))';end;
%
if size(D,1)>size(u,1);D=D(1: size(u,1),1: size(u,1));end;
% else;D=D;end;
%
wh = pinv(u'*D*u); v_GMMH=pinv(zx'*wh*zx); beta_GMMH=v_GMMH*(zx'*wh*zy);   
r_GMMH=y-x*beta_GMMH;se_GMMH=sqrt(sum(r_GMMH.^2)/n);sigmasq_GMMH=sum(r_GMMH.^2)/n;
Haus_GMMH=(beta_OLS-beta_GMMH)'*pinv(v_GMMH-v_OLS)*(beta_OLS-beta_GMMH);pval_Haus_GMMH=1-cdf('chi2',Haus_GMMH,k);
u=bsxfun(@times,z,r_GMMH);  meanu=mean(u);
xe=bsxfun(@times,x,r_GMMH);meanxe=mean(sqrt(xe(2:n,:).^2));
AIC_GMMH=log(sum(r_GMMH.^2/n))+2*k/n;BIC_GMMH=log(sum(r_GMMH.^2/n))+log(n)*k/n;
DW_GMMH = sum(diff(r_GMMH).^2)/sum(r_GMMH.^2); [~,ARCH_GMMH_pval,ARCH_GMMH,critval] = archtest(r_GMMH,1,.05);  %Engle ARCH test of residuals, 1 lag
CR=corrcoef(r_GMMH(2:end),r_GMMH(1:end-1));AU_GMMH=CR(1,2); % First-order autocorrelation coefficient
betaseh=sqrt(diag(v_GMMH));tstats_GMMH=beta_GMMH./betaseh;
%j_GMMH=n*meanu*wh*meanu'; 
sigmasq_GMMH=sum(r_GMMH.^2)/n;j_GMMH=r_GMMH'*P*r_GMMH/sigmasq_GMMH;j_GMMH_ASY=(j_GMMH-q)/sqrt(2*q); % Standard and asymptotic Jstat (Andrews & Stock, 2007)
pj_GMMH=1.-cdf('chi2',j_GMMH, q-k);pj_GMMH_ASY=1-cdf('Normal',j_GMMH_ASY,0,1);  
AR_GMMH=(n-q)*(r_GMMH'*P*r_GMMH)/(r_GMMH'*M*r_GMMH); p_AR_GMMH=1.-chi2cdf(AR_GMMH, q);  % Standard and asymptotic AR (Anderson-Rubin) test for betas=beta_GMMH.
AR_GMMH_ASY=sqrt(q)*((AR_GMMH/q)-1.);pAR_GMMH_ASY=1-cdf('Normal',AR_GMMH_ASY,0,2);
fox_GMMH=G'*wh*meanu';
% Computing H,P,a and biases BI, BO, BG (Newey & Smith 2004). Uses their notation.
omega=u'*D*u/n;sigma=pinv(G'*pinv(omega)*G); CP_GMMH=min(eig(G'*pinv(omega)*G));    
H_GMMH=sigma*G'*pinv(omega);P_GMMH=pinv(omega)-pinv(omega)*G*H_GMMH;GPg_GMMH=meanu*P_GMMH*meanu';
avec=zeros(q,1);xzj=zeros(k,q);for j=1:q;xzj(:,j)=-(x/n)'*z(:,j);avec(j)=trace(sigma*xzj(:,j)*xzj(:,j)')/(2*n);end
BIGH=H_GMMH*(-avec+G*H_GMMH*meanu')/n;BOGH=H_GMMH*(meanu*meanu'*P_GMMH*meanu')/n;BGGH=-sigma*(((q-k)*meanxe/(se_GMMH^2)))'/n;
betastats_gmmh=[beta_GMMH,tstats_GMMH];
stats_gmmh=[se_GMMH,j_GMMH,pj_GMMH,j_GMMH_ASY,pj_GMMH_ASY,AR_GMMH,p_AR_GMMH,AR_GMMH_ASY,pAR_GMMH_ASY,DW_GMMH,ARCH_GMMH,ARCH_GMMH_pval,AU_GMMH];
bias_gmmh=[BIGH,BOGH,BGGH];
%
% Jackknife 2-step GMM. Ref. Newey-Windmeijer, Econometrica, 2009.
%
% First step (Jackknife TSLS: JIVE2). Ref. Angrist, Imbens & Krueger, AIK (1999).
%
xr=repmat(x,1,n);zr=repmat(z,1,n);
for i=1:n
A=xr(:,(i-1)*k+1:i*k); A(i,:)=[]; xj=A(:,1:end);      % x stripped off of one row for each i=1:n
B=zr(:,(i-1)*q+1:i*q); B(i,:)=[]; zj=B(:,1:end);      % z stripped off of one row for each i=1:n
% xhats2: nxk matrix (eq. 6, AIK) where each row includes xhat's ith.observations
xhats2(i,:)=z(i,:)*pinv(z'*z)*(z'*x-zj'*xj)/(1-1/n);
end
beta_JACK2=pinv(xhats2'*x)*xhats2'*y;
%
% Second step
%
r_JTSLS=y-x*beta_JACK2; % r_JTSLS=r_GMM;
n2=n^2;sigmasq_JTSLS=sum(r_JTSLS.^2)/n;
u=bsxfun(@times,z,r_JTSLS);
xr=repmat(x,1,n);zr=repmat(z,1,n);ur=repmat(u,1,n);rr=repmat(r_JTSLS,1,n);
% For x,z,u produce time series of length n-1 that, for each i.th observation (1=1,...n), exclude the i.th row:
for i=1:n
A=xr(:,(i-1)*k+1:i*k); A(i,:)=[]; xj=A(:,1:end);      % x stripped off of one row for each i=1:n
B=zr(:,(i-1)*q+1:i*q); B(i,:)=[]; zj=B(:,1:end);      % z stripped off of one row for each i=1:n
C=ur(:,(i-1)*q+1:i*q); C(i,:)=[]; uj=C(:,1:end);      % u stripped off of one row for each i=1:n
D=rr(:,(i-1)+1:i); D(i)=[];rj=D(:,1:end);
% xhatsg: nxk matrix (eq. 6 of AIK GMM-adjusted) where each row includes xhat's ith.observations
xhatsg(i,:)=z(i,:)*pinv(u'*u)*(z'*x-zj'*xj)/(1-1/n);
hj=xj'*zj*pinv(u'*u)*zj'*xj;  % HJ(i).a=hj;
%pj=zj*pinv(zj'*zj)*zj';
pj=z*pinv(zj'*zj)*z';PJ(i).a=pj;MJ(i).a=eye(n)-PJ(i).a;
%ej=(rj'*zj*pinv(zj'*zj)*zj'*rj)/sum(rj.^2/nu);EJ(i)=ej;
ej=(rj'*zj*pinv(zj'*zj)*zj'*rj)/sigmasq_JTSLS;EJ(i)=ej;
wj=pinv(uj'*uj);WJ(i).a=wj;
Gj=-xj'*zj/n;lambdaj=(Gj*pinv(u'*u)*(uj'*uj)*pinv(u'*u)*Gj')/(n2*(n-1));
vj=pinv(hj)*(x'*z*pinv(u'*u)*z'*x+lambdaj)*pinv(hj);  VJ(i).a=vj/n;   %Asymptotic variance of betas
end
beta_JGMM=pinv(xhatsg'*x)*xhatsg'*y;
r_JGMM=y-x*beta_JGMM; se_JGMM=sqrt(sum(r_JGMM.^2)/n);sigmasq_JGMM=sum(r_JGMM.^2)/n;
for i=1:n;PP(i,:)=reshape(PJ(i).a',1,numel(PJ(i).a));end;meansPJ=reshape(mean(PP),size(PJ(1).a));
for i=1:n;CM(i,:)=reshape(MJ(i).a',1,numel(MJ(i).a));end;meansMJ=reshape(mean(CM),size(MJ(1).a));
for i=1:n;CV(i,:)=reshape(VJ(i).a',1,numel(VJ(i).a));end;meansVJ=reshape(mean(CV),size(VJ(1).a));
for i=1:n;CW(i,:)=reshape(WJ(i).a',1,numel(WJ(i).a));end;meansWJ=reshape(mean(CW),size(WJ(1).a));
% clear PP CM CV CW
v_JGMM=pinv(zx'*meansWJ*zx);
Haus_JGMM=(beta_OLS-beta_JGMM)'*pinv(v_JGMM-v_OLS)*(beta_OLS-beta_JGMM);pval_Haus_JGMM=1-cdf('chi2',Haus_JGMM,k);
betaseh=sqrt(diag(v_JGMM));tstats_JGMM=beta_JGMM./betaseh;
AIC_JGMM=log(sum(r_JGMM.^2/n))+2*k/n;BIC_JGMM=log(sum(r_JGMM.^2/n))+log(n)*k/n;
DW_JGMM= sum(diff(r_JGMM).^2)/sum(r_JGMM.^2); [~,ARCH_JGMM_pval,ARCH_JGMM,critval] = archtest(r_JGMM,1,.05);  %Engle ARCH test of residuals, 1 lag
CR=corrcoef(r_JGMM(2:end),r_JGMM(1:end-1));AU_JGMM=CR(1,2); % First-order autocorrelation coefficient
u=bsxfun(@times,z,r_JGMM);  %meanu=mean(u);j_JGMM=n*meanu*meansWJ*meanu'; 
j_JGMM=mean(EJ); % instead of j_JGMM=r_JGMM'*P*r_JGMM/sigmasq_JGMM;
j_JGMM_ASY=(j_JGMM-q)/sqrt(2*q); % Standard and asymptotic Jstat (Anatolyev-Gospodinov, 2008)
pj_JGMM=1.-cdf('chi2',j_JGMM, q-k);pj_JGMM_ASY=1-cdf('Normal',j_JGMM_ASY,0,1);  
AR_JGMM=(n-q)*(r_JGMM'*meansPJ*r_JGMM)/(r_JGMM'*meansMJ*r_JGMM); p_AR_JGMM=1.-chi2cdf(AR_JGMM, q);  % Standard and asymptotic AR (Anderson-Rubin) test for betas=beta_JGMM.
AR_JGMM_ASY=sqrt(q)*((AR_JGMM/q)-1.);pAR_JGMM_ASY=1-cdf('Normal',AR_JGMM_ASY,0,2);
fox_JGMM=G'*meansWJ*meanu';
% Computing H,P,a and biases BI, BO, BG (Newey & Smith 2004). Uses their notation.
omega=u'*u/n;sigma=pinv(G'*pinv(omega)*G); CP_JGMM=min(eig(G'*pinv(omega)*G));    
H_JGMM=sigma*G'*pinv(omega);P_JGMM=pinv(omega)-pinv(omega)*G*H_JGMM;GPg_JGMM=meanu*P_JGMM*meanu';
avec=zeros(q,1);xzj=zeros(k,q);for j=1:q;xzj(:,j)=-(x/n)'*z(:,j);avec(j)=trace(sigma*xzj(:,j)*xzj(:,j)')/(2*n);end
BIGJ=H_JGMM*(-avec+G*H_JGMM*meanu')/n;BOGJ=H_JGMM*(meanu*meanu'*P_JGMM*meanu')/n;BGGJ=-sigma*(((q-k)*meanxe/(se_JGMM^2)))'/n;
betastats_jgmm=[beta_JGMM,tstats_JGMM];
stats_jgmm=[se_JGMM,j_JGMM,pj_JGMM,j_JGMM_ASY,pj_JGMM_ASY,AR_JGMM,p_AR_JGMM,AR_JGMM_ASY,pAR_JGMM_ASY,DW_JGMM,ARCH_JGMM,ARCH_JGMM_pval,AU_JGMM];
bias_jgmm=[BIGJ,BOGJ,BGGJ];
%
res.hac=bw;
res.mxe=[meanxe,GPg_GMM,GPg_GMMH,GPg_JGMM];
res.granger= [ecx,c_v;xce,c_v];
res.bze=[bz_fstats,bz_fpvals];
res.RF=[fs_rhosq,fs_fstats,fs_fpvals];res.PIE=PIE;
res.erCD=[ErankPIE,minsvPIE,CP,CD];
res.jacob=[mineigG,minsvdG];
res.inst=[CP_GMM,CP_GMMH,CP_JGMM];
res.betas=[beta_GMM, beta_GMMH, beta_JGMM];
res.tstats=[tstats_GMM, tstats_GMMH, tstats_JGMM];
res.haus=[Haus_GMM,Haus_GMMH,pval_Haus_GMM,pval_Haus_GMMH];
res.genstats=[stats_gmm;stats_gmmh;stats_jgmm];
res.aicbic=[AIC_GMM,BIC_GMM;AIC_GMMH,BIC_GMMH;AIC_JGMM,BIC_JGMM];
res.bias=[bias_gmm,bias_gmmh,bias_jgmm];
res.foxes=[fox_GMM,fox_GMMH,fox_JGMM];
res.resids=[r_GMM,r_GMMH,r_JGMM];
%
%
%
end

