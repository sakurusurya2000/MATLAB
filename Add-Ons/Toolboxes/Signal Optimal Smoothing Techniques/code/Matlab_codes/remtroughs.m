function [tks,tlx] = remtroughs(tks,tlx,Pd)
%  Remove Troughs Separated By Less Than Min Trough Distance.
%  Function: [tks,tlx] = remtroughs(tks,tlx,Pd);
%  Start with the larger troughs to make sure we don't accidentally keep a small trough and remove a large trough in its neighborhood. 
%  INPUT: 
%  (1) tks: values of troughs available, obtained by peaks_and_troughs.m 
%  (2) tlx: location of troughs in (1)
%  (3) Pd: minimum affordable distance between troughs 
%
%  OUTPUT:
%  (1) tks: values of troughs found sorted in descending order
%  (2): tlx: location of troughs found
% 
if isempty(tks) || Pd==1,
    return
end
% Order troughs from large to small
[tks, idx] = sort(tks,'ascend');
tlx = tlx(idx);
idelete = ones(size(tlx))<0;
for i = 1:length(tlx),
    if ~idelete(i),
        % If the trough is not in the neighborhood of a larger trough, find
        % secondary troughs to eliminate.
        idelete = idelete | (tlx>=tlx(i)-Pd)&(tlx<=tlx(i)+Pd); 
        idelete(i) = 0; % Keep current trough
    end
end
tks(idelete) = [];
tlx(idelete) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


