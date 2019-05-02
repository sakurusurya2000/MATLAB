function [coins, coind] = tbreakclust(tbdates,clusterdates)
%
%  Finds how many and which time breaks coincide with clusterdates
%  Function is: [coins, coind] = tbreakclust(tbdates,clusterdates)
%
% INPUTS: time breakdates and clusterdate vectors, may differ in length
% OUTPUTS: # of coincident dates and their location, may be zero
%
for j=1:numel(tbdates);t1plus=tbdates(j)+(1:5);t1mis=fliplr(tbdates(j)-(1:5));t1all=[t1mis,tbdates(j),t1plus];TALL(j,:)=t1all;end
TALL1=TALL';TALL1(:);
for i=1:numel(clusterdates);[~,b]=find(clusterdates(i)==TALL1);bs(i)=isempty(b);end
b2=(bs-1).^2;
coins=sum(b2>0);
coind=clusterdates(b2>0);
%
end

