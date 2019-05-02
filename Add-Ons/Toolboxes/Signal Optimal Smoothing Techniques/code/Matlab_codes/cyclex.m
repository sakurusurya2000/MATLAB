function cycle= cyclex(x)
% Finds HI-RES mean cycle of series x by means of FFT
% Formula: cycle=cyclex(x)
%   
T=numel(x);
fy=abs(fft(x));
screep=fy(1:ceil(.5*T)); % Build screeplot of FFT of x
cycle= ERE(screep);
%
end

