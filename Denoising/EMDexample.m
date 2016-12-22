clc
clear all
close all
load ecg1;
xc=ecg;
N=length(xc);
x=xc+sqrt(60)*randn(1,N);
Fs=1000;
imf = EMD(x);
plot_hht(x,imf,1/Fs);
[y,imfp]=EMDFilter(x,imf);
