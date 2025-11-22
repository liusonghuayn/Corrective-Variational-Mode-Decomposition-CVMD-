clc
clear
close all

Fs = 100;
Ts = 1/Fs;
t  = 1/Fs:1/Fs:710/Fs;
T  = length(t);

x2 = cos(2*pi*7*t);
x1 = cos(2*pi*(20-0.2*t).*t);
x3 = cos(2*pi*(2+0.1*t).*t);
xn = x1 + x2 + x3 + randn(1,T).*0.3;

kmax = 10;
[Kpro, IMF_best, omega_best] = CVMD(xn, Fs, kmax);

% 画出原信号 + 各 IMF，所有 IMF 的 y 轴范围一致
plotCVMD_IMFs(t, xn, IMF_best, 'CVMD decomposition');