 function [pe,tr]=peaks_and_troughs(x)
 %
% Finds peaks and throughs of time series x, after Nagi Hatoum, Matlab File Exchange, copyright 2005
% Both samples must be of equal length. If not so, the longer sample is stripped off of rando location elements
ds=diff(x);
ds=[ds(1);ds];%pad diff
filter=find(ds(2:end)==0)+1;% find zeros
ds(filter)=ds(filter-1);% replace zeros
ds=sign(ds);
ds=diff(ds);
tr=find(ds>0); % Troughs
pe=find(ds<0); % Peaks   
spe=size(pe,1);str=size(tr,1);sd=abs(spe-str); % Measuring length of pe, tr
 if spe>str;
     rdd=sort(randsample(pe,sd));
     [~,IA] = intersect(pe,rdd);
     pe(IA)=[];
 else
     rdd=sort(randsample(tr,sd));
     [~,IA] = intersect(tr,rdd);
     tr(IA)=[];
 end
 %
 %
 %
  end
