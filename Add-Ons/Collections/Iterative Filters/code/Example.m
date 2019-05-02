
% Example 1
% 
% Equation (37) on page 14
%
% Ref: A. Cicone, J. Liu, H. Zhou. 'Adaptive Local Iterative Filtering for 
%      Signal Decomposition and Instantaneous Frequency analysis'. 
%      http://arxiv.org/abs/1411.6051
%

dt=0.001;

t=0:0.001:1;

x=(2*(t-0.5).^2+0.2).*sin(20*pi*t+0.2*cos(40*pi*t));

y=4*(t-0.5).^2;

z=x+y+1;

figure
plot(t,z)

%% Using a constant extension of the signal

opts=Settings_IFs('IFs.delta',0.001,'IFs.NIMFs',1,'IFs.Xi',1.9,'IFs.alpha',1);
IMF = IFs(z,opts);


plot_imf_v8(IMF);

%% Using a periodic extension of the signal

opts=Settings_IFs('IFs.delta',0.001,'IFs.NIMFs',1,'IFs.Xi',1.9,'IFs.alpha',1,'IFs.extensionType','p');
IMF = IFs(z,opts);


plot_imf_v8(IMF);

%% Using reflection

opts=Settings_IFs('IFs.delta',0.001,'IFs.NIMFs',1,'IFs.Xi',1.9,'IFs.alpha',1,'IFs.extensionType','r');
IMF = IFs(z,opts);


plot_imf_v8(IMF);




