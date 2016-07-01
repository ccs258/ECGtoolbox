clear all
close all
clc

load('ecg.mat');
s0=ecg(1:2000);
N=length(s0);
sigma=60;
s2=s0+sqrt(sigma)*randn(1,N);

N=5;
c=medfilt1(s2,N);
ecgRecMED=c;
[P1,Q1,R1,S1,T1,Pon1,Poff1,QRSon1,QRSoff1,ST1,Toff1]=pqrstDetect(ecgRecMED,500);

lev=5;
c=wden(s2,'minimaxi','s','sln',lev,'sym8');
ecgRecWAV=c;
[P2,Q2,R2,S2,T2,Pon2,Poff2,QRSon2,QRSoff2,ST2,Toff2]=pqrstDetect(ecgRecWAV,500);
