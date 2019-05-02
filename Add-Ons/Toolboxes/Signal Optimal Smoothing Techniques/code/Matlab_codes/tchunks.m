function [num,tcs,td1] = tchunks(td)
%
% Finds time intervals among a given vector of random time chunks, including last observation which must be an input.
% Function: [num,tcs] = tchunks(td)
% INPUT: td:(num x 1) vector of time chunks including end of sample
% OUTPUT: (1) # of time chunks, (2) 1 x num struct array
%  Function is: [num,tcs] = tchunks(td);
% EXAMPLE: find tchunks lying between td=[45 78 99 110] where 110 is last observation of underlying signal:
% tcs(1).a=(1,2,...45), tcs(2).a=(46,47,..,78),...,tcs(num).a=(100,101,...,110)
%   
if size(td,1)==1;td=td';end; %Make sure data are entered columnwise
td1=td+1;
TD=[td,td1];
TD2=TD(:,1);
TD1=[1;TD(1:end-1,2)];
TTD=[TD1,TD2];
a=find(TTD(:,1)==TTD(:,2));
if isempty(a)==1;TTD(a,:)=[];end
num=(size(TTD,1));
for i=1:num;tcs(i).a=(TTD(i,1):TTD(i,2));end
td1=(TTD(:,2))';
%
end