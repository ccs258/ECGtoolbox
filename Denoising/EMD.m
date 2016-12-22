function imf = EMD(x)
% Empiricial Mode Decomposition (Hilbert-Huang Transform)
%
% Refs:
% N. Huang, Z. Shen, S. Long, M. Wu, H. SHIH, Q. ZHENG, N. Yen, C. Tung, and H. Liu, 
%¡°The empirical mode decomposition and the Hilbert spectrum 
% for nonlinear and non-stationary time series analysis,¡± 
% Proc. R. Soc. A Math. Phys. Eng. Sci., vol. 454, no. 1971, pp. 995, 903, 1998.
%%
% Input: original signal
% Output: first IMF, second IMF, ... last residual
%         in cell class
%%
x   = transpose(x(:));
imf = [];
while ~ismonotonic(x)
    x1 = x;
    sd = Inf;
    while (sd > 0.1) || ~isimf(x1)
        s1 = getspline(x1);         % peak spline
        s2 = -getspline(-x1);       % valley spline
        x2 = x1-(s1+s2)/2;
       
        sd = sum((x1-x2).^2)/sum(x1.^2);
        x1 = x2;
    end
   
    imf{end+1} = x1;
    x          = x-x1;
end
imf{end+1} = x;
% decide if x is monotonic
function u = ismonotonic(x)
u1 = length(findpeaks(x))*length(findpeaks(-x));
if u1 > 0
    u = 0;
else
    u = 1;
end
% decide if x is IMF
function u = isimf(x)
N  = length(x);
u1 = sum(x(1:N-1).*x(2:N) < 0);                  % number of zero-crossings
u2 = length(findpeaks(x))+length(findpeaks(-x)); % number of extreme points
if abs(u1-u2) > 1
    u = 0;
else
    u = 1;
end
% create spline based on peaks
function s = getspline(x)
N = length(x);
p = findpeaks(x);
s = spline([0 p N+1],[0 x(p) 0],1:N);

function n = findpeaks(x)
% Find peaks. 
n    = find(diff(diff(x) > 0) < 0); % second derivative<0
u    = find(x(n+1) > x(n));
n(u) = n(u)+1;                      

