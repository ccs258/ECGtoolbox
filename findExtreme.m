function [nmax,nmin] = findExtreme(x)
%% Description
%find the extreme points of the signal x
%returns the index of the extreme points of x

%% Inputs
% x : signal

%% Outputs
%nmax returns the index of the max extreme points of x
%nmin returns the index of the min extreme points of x

%% Author : Lisha Chen    
% contact : lishachen00@gmail.com , lisha_chen@hust.edu.cn
% Dont forget to reference if you found this script usefull

%%
nmax=[];
nmin=[];


for i=2:length(x)-1
if x(i-1)<=x(i)&&x(i+1)<=x(i)&&((x(i-1)-x(i))^2+(x(i+1)-x(i))^2)~=0
    nmax=[nmax,i];
end

if x(i-1)>=x(i)&&x(i+1)>=x(i)&&((x(i-1)-x(i))^2+(x(i+1)-x(i))^2)~=0
    nmin=[nmin,i];
end
end

%% plot

% plot(x);
% hold on
% plot(nmax,x(nmax),'r*');
% 
% figure,
% plot(x);
% hold on
% plot(nmin,x(nmin),'r*');