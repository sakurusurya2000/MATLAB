function [mdate,ert]= ERE(X)
%  Optimal stopping time (OST) of screeplot applying a modification of Ahn-Horenstein's ERE method to a single signal X(T,1) 
%  Function is [mdate, ert]=ERE(X)
%  Mdate is OST which maximizes the normed ratio of two adjacent screeplot values, ert is the screeplot 
%  Reference: Ahn S.C. and Horenstein A.R., 2013, "Eigenvalue Ratio Test for the Number of Factors", Econometrica, 81, 12031227.
%
%   
T=numel(X);
ert=[];
for j=3:T;ert(j-2)=norm(X(j-2)-X(j-1))/norm(X(j-1)-X(j));end;
[~,mdate]=max(ert);
%
end

