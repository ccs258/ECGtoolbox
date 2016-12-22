function plot_hht(x,imf,Ts)
% Plot the HHT.
% :: Syntax
%    The array x is the input signal and Ts is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd
% imf = emd(x);
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});
    th   = unwrap(angle(hilbert(imf{k})));  % 相位
    d{k} = diff(th)/Ts/(2*pi);          % 瞬时频率
end

% display IMF
M = length(imf);
N = length(x);
c = linspace(0,(N-1)*Ts,N);             % 0:Ts:Ts*(N-1)
for k1 = 0:4:M-1
    figure
    for k2 = 1:min(4,M-k1)
        subplot(4,2,2*k2-1)
        plot(c,imf{k1+k2})
        set(gca,'FontSize',8,'XLim',[0 c(end)]);
        title(sprintf('No.%d IMF', k1+k2))
        xlabel('Time/s')
        ylabel(sprintf('IMF%d', k1+k2));
       
        subplot(4,2,2*k2)
        [yf, f] = FFTAnalysis(imf{k1+k2}, Ts);
        plot(f, yf)
        set(gca,'FontSize',8,'XLim',[0 c(end)]);
        title(sprintf('No.%d IMF frequency spectrum', k1+k2))
        xlabel('f/Hz')
        ylabel('|IMF(f)|');
    end
end
figure
subplot(211)
plot(c,x)
set(gca,'FontSize',8,'XLim',[0 c(end)]);
title('original signal')
xlabel('Time/s')
ylabel('Origin');
subplot(212)
[Yf, f] = FFTAnalysis(x, Ts);
plot(f, Yf)
set(gca,'FontSize',8,'XLim',[0 c(end)]);
title('frequency spectrum of original signal')
xlabel('f/Hz')
ylabel('|Y(f)|');

end

function [Y, f] = FFTAnalysis(y, Ts)
Fs = 1/Ts;
L = length(y);
NFFT = 2^nextpow2(L);
y = y - mean(y);
Y = fft(y, NFFT)/L;
Y = 2*abs(Y(1:NFFT/2+1));
f = Fs/2*linspace(0, 1, NFFT/2+1);
end


% % Hilbert analysis
% function [yenvelope, yf, yh, yangle] = HilbertAnalysis(y, Ts)
% yh = hilbert(y);
% yenvelope = abs(yh);                % envelop
% yangle = unwrap(angle(yh));         % phase
% yf = diff(yangle)/2/pi/Ts;          % instantaneous frequency
% end