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

