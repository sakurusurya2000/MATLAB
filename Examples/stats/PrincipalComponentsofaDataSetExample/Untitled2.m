clc
clear all
close all
t=0:.0001:5;
f=1;
x=sin(2*pi*f*t);
%fs=2;
%ts=0:1/fs:5;
%xs=sin(2*pi*f*ts);
%plot(x)
imf = tvf_emd(x);
%for i=1:50001
%a1(1, i)=imf(1, i);
%a2(1, i)=imf(2, i);
%a3(1, i)=imf(3, i);
%a4(1, i)=imf(4, i);
%a5(1, i)=imf(5, i);
%a6(1, i)=imf(6, i);
%a7(1, i)=imf(7, i);
%end
%figure
%plot(a1)
%figure
%plot(a2)
%figure
%plot(a3)
%figure
%plot(a4)
%figure
%plot(a5)
%figure
%plot(a6)
%figure
%plot(a7)
x=transpose(imf);
level = 5;
wname = 'sym4';
npc = 'kais';
[x_sim, qual, npc] = wmspca(x ,level, wname, npc);
kp = 0;
for i = 1:7
    subplot(7,2,kp+1), plot(x (:,i)); axis tight;
    title(['Original signal ',num2str(i)])
    subplot(7,2,kp+2), plot(x_sim(:,i)); axis tight;
    title(['Simplified signal ',num2str(i)])
    kp = kp + 2;
end

coeff = pca(imf);
for i=1:50001
    b1(i, 1)=coeff(i, 1);
    b2(i, 1)=coeff(i, 2);
    b3(i, 1)=coeff(i, 3);
    b4(i, 1)=coeff(i, 4);
    b5(i, 1)=coeff(i, 5);
    b6(i, 1)=coeff(i, 6);
end
figure
plot(b1)
figure
plot(b2)
figure
plot(b3)
figure
plot(b4)
figure
plot(b5)
figure
plot(b6)
