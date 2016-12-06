function [signal_WWF]=WWF(input_signal)

% wavelet wiener filter
% filter is designed for filtering signals with sampling frequency of 500 Hz
% Input: the noisy signal
% Output: filtered signal

% Refs:
% Wavelet Wiener filter of ECG signals
% Eva Sedlackova
% Bc. Joshua Jan??2014 VUT FEKT Brno

% Author: Lukas Smital
% Access: http://www.ubmi.feec.vutbr.cz/vyzkum-a-vyvoj/produkty
% Commented by: Lisha Chen
 
% filtering parameters to calculate the estimate input SNR
  stupen_rozkladu_WT1=3; % the degree of decomposition of WT1
  stupen_rozkladu_WT2=3; % the degree of decomposition of WT2
  vlnka_WT1='bior3.9';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.3;  

[vstupni_signal_prodlouzeny,delka_prodlouzeni]=prodlouzeni_automaticke(input_signal,5);   
% extended signal using mirrored cone to eliminate transient phenomenon
% and lengths for achieving the desired function swt

% /// Prost wavelet filtering - a prior signal estimation u?ite¨¨n¨¦ho///

matice_rozklad_po=swt(vstupni_signal_prodlouzeny,stupen_rozkladu_WT1,vlnka_WT1);  
% WT1: performs discrete wavelet transformation of za?um¨¬n¨¦ho signal 
% with the selected degree in decomposition and the chosen type of filter

for i=1:stupen_rozkladu_WT1 
   threshold=prah_adaptivni_empiricky(matice_rozklad_po(i,:))*empiricka_konstanta; 
   % calculate the threshold
   threshold=conv(threshold,ones(1,25)*0.04,'same'); 
   % convolution window size 25 samples
   matice_rozklad_po(i,:)=prahovani_hybridni(matice_rozklad_po(i,:),threshold); 
   % thresholding matrix decomposition
end

% IWT1: Inverse wavelet transformation
% of each band with WT coefficients after thresholding
signal_PO=iswt(matice_rozklad_po,vlnka_WT1);

% /// Calculation and application of Wiener korek¨¨n¨ªho factor ///

matice_rozklad_wf_po=swt(signal_PO,stupen_rozkladu_WT2,vlnka_WT2); 
%WT2 prior signal estimation u?ite¨¨n¨¦ho
matice_rozklad_wf_os=swt(vstupni_signal_prodlouzeny,stupen_rozkladu_WT2,vlnka_WT2); 
%WT2 zaru?en¨¦ho signal to obtain a noise variance

for i=1:stupen_rozkladu_WT2
 
koeficienty_pasmo_i_po=matice_rozklad_wf_po(i,:);
koeficienty_pasmo_i_os=matice_rozklad_wf_os(i,:);

% Calculating the standard deviation of noise in a floating window
s_o_sumu=zeros(length(koeficienty_pasmo_i_os),1);
koeficienty_pasmo_i_os=prodlouzeni(koeficienty_pasmo_i_os,500); %extend?(*)
for z=1:(length(koeficienty_pasmo_i_os)-500)
    s_o_sumu(z)=smerodatna_odchylka_sumu(koeficienty_pasmo_i_os(z:z+500));
end
koeficienty_pasmo_i_os=zkraceni(koeficienty_pasmo_i_os,500); % truncate (*)
 
%HW: calculation and application of korek¨¨n¨ªho factor
 for x=1:length(koeficienty_pasmo_i_po)   
     koeficienty_pasmo_i_os(x)=koeficienty_pasmo_i_os(x)*...
         (((koeficienty_pasmo_i_po(x))^2)/...
         ((((koeficienty_pasmo_i_po(x))^2)+(s_o_sumu(x))^2)));
 end
     matice_rozklad_wf_os(i,:)=koeficienty_pasmo_i_os;
end

signal_WWF=iswt(matice_rozklad_wf_os,vlnka_WT2); 
%IWT2: inverse wavelet transformation after application of korek¨¨n¨ªho factor
signal_WWF=zkraceni(signal_WWF,delka_prodlouzeni); 
% reduced signal after filtering WWF to its original length
                                  
SNR_odhad=10*log10(sum(signal_WWF.^2)/sum((input_signal-signal_WWF).^2)); 
% SNR achieved after filtration WWF
vstupni_SNR=5*round(SNR_odhad/5); % Input SNR rounded to multiples of 5

% thresholding the input SNR
% when the input SNR did not fall within the range for which the filter is ready
if       vstupni_SNR<-5
         vstupni_SNR=-5;
elseif   vstupni_SNR>50
         vstupni_SNR=50;
end

% if the estimate of the input SNR~=45 will perform filtration
% the appropriate parameters for the estimated input level disturbed
if vstupni_SNR~=45

% define filter parameter for each level of input SNR -5 db and 50 db 
switch vstupni_SNR
  case {-5}
  stupen_rozkladu_WT1=5; % the degree of decomposition of WT1
  stupen_rozkladu_WT2=5; % the degree of decomposition of WT2
  vlnka_WT1='rbio3.3';
  vlnka_WT2='rbio4.4';
  empiricka_konstanta=3.1;
  case {0}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=5;
  vlnka_WT1='rbio1.3';
  vlnka_WT2='rbio4.4';
  empiricka_konstanta=3.7;
  case {5}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=5;
  vlnka_WT1='rbio1.3';
  vlnka_WT2='rbio4.4';
  empiricka_konstanta=3.2;
  case {10}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=5;
  vlnka_WT1='rbio1.3';
  vlnka_WT2='rbio4.4';
  empiricka_konstanta=3.1;
  case {15}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=5;
  vlnka_WT1='db4';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.6;
  case {20}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=4;
  vlnka_WT1='db4';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.8;
  case {25}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=4;
  vlnka_WT1='bior4.4';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.5;
  case {30}
  stupen_rozkladu_WT1=4;
  stupen_rozkladu_WT2=4;
  vlnka_WT1='bior4.4';
  vlnka_WT2='sym4';
  empiricka_konstanta=1.9;
  case {35}
  stupen_rozkladu_WT1=3;
  stupen_rozkladu_WT2=4;
  vlnka_WT1='bior4.4';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.5;
  case {40}
  stupen_rozkladu_WT1=3;
  stupen_rozkladu_WT2=3;
  vlnka_WT1='bior3.9';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.6;
  case {45}
  stupen_rozkladu_WT1=3;
  stupen_rozkladu_WT2=3;
  vlnka_WT1='bior3.9';
  vlnka_WT2='sym4';
  empiricka_konstanta=2.3;
  case {50}
  stupen_rozkladu_WT1=2; 
  stupen_rozkladu_WT2=3;
  vlnka_WT1='sym6';
  vlnka_WT2='bior3.3';
  empiricka_konstanta=2.5;
end


% /// Prost wavelet filtering - a prior signal estimation u?ite¨¨n¨¦ho///

matice_rozklad_po=swt(vstupni_signal_prodlouzeny,stupen_rozkladu_WT1,vlnka_WT1);  
% WT1: performs discrete wavelet transformation of za?um¨¬n¨¦ho signal 
% with the selected degree in decomposition 
% and chosen type of ripples on the prior signal estimation u?ite¨¨n¨¦ho

% Shortening the length of the belt wavelet
% coefficients and threshold values ??stored in the output variables


for i=1:stupen_rozkladu_WT1 
 
% calculate the threshold
  threshold=prah_adaptivni_empiricky(matice_rozklad_po(i,:))*empiricka_konstanta;
  threshold=conv(threshold,ones(1,25)*0.04,'same'); 
  % destroyed during the threshold pr¨´m¨¬rov¨¢n¨ªm window size 25

% thresholding wavelet coefficients
  matice_rozklad_po(i,:)=prahovani_hybridni(matice_rozklad_po(i,:),threshold);
    
end
 
% IWT1: inverse wavelet transformation of individual zones with the WT coefficients after thresholding?
signal_PO=iswt(matice_rozklad_po,vlnka_WT1);

% /// Calculation and application of Wiener koekeniho factor ///

matice_rozklad_wf_po=swt(signal_PO,stupen_rozkladu_WT2,vlnka_WT2); 
% WT2 prior estimate uziteeneho signal
matice_rozklad_wf_os=swt(vstupni_signal_prodlouzeny,stupen_rozkladu_WT2,vlnka_WT2); 
% WT2 zaruseneho signal

for i=1:stupen_rozkladu_WT2
 
koeficienty_pasmo_i_po=matice_rozklad_wf_po(i,:);
koeficienty_pasmo_i_os=matice_rozklad_wf_os(i,:);

% calculating standard deviation of noise in a floating window
s_o_sumu=zeros(length(koeficienty_pasmo_i_os),1);
koeficienty_pasmo_i_os=prodlouzeni(koeficienty_pasmo_i_os,500); % extend?(*)
for z=1:(length(koeficienty_pasmo_i_os)-500)
    s_o_sumu(z)=smerodatna_odchylka_sumu(koeficienty_pasmo_i_os(z:z+500));
end
koeficienty_pasmo_i_os=zkraceni(koeficienty_pasmo_i_os,500); % truncat?(*)
 
% HW: calculate and apply Wiener korekeniho factor
 for x=1:length(koeficienty_pasmo_i_po)   
     koeficienty_pasmo_i_os(x)=koeficienty_pasmo_i_os(x)*...
         (((koeficienty_pasmo_i_po(x))^2)/...
         ((((koeficienty_pasmo_i_po(x))^2)+(s_o_sumu(x))^2)));
 end
     matice_rozklad_wf_os(i,:)=koeficienty_pasmo_i_os;
end

signal_WWF=iswt(matice_rozklad_wf_os,vlnka_WT2); 
% IWT2: inverse wavelet transformation after apply korekeniho factor

signal_WWF=zkraceni(signal_WWF,delka_prodlouzeni);  
% after filter WWF, truncate the signal

end
end

function [threshold]=prah_adaptivni_empiricky(waveletCoefficient)
% empirical adaptive threshold
% calculation of adaptive empirical floating window sill
waveletCoefficient=prodlouzeni(waveletCoefficient,300);

delka_sig = length(waveletCoefficient);

threshold=ones(delka_sig,1);

for i=1:(delka_sig-300)
    threshold(i+150)=smerodatna_odchylka_sumu(waveletCoefficient(i:300+i)); 
    % estimate of the standard deviation
end

threshold=zkraceni(threshold,300);
end

function [waveletCoefficient]=prahovani_hybridni(waveletCoefficient,threshold)
% Hybrid thresholding (non-negative garrotte) wavelet coefficients

pocet_vzorku=length(waveletCoefficient); % number of samples

for i=1:pocet_vzorku
   
    % realize relationship thresholding for hybrid / non-negative garrotte /
    
    if abs(waveletCoefficient(i))<=threshold(i)  
       waveletCoefficient(i)=0;
    else
       waveletCoefficient(i)=waveletCoefficient(i)-(threshold(i))^2/waveletCoefficient(i);
    end
    
end
end

function [output_signal]=prodlouzeni(input_signal,extendLen)
% Extend signal at beginning and end with the desired length
%
% Extend signal by mirrored ends in suppressing PURPOSE
% eg. prodlouzeni(signal,18)
% signal will be extended by 9 samples at the beginning and at the end
% must be a barrel number

output_signal=zeros(1,length(input_signal)+extendLen);
for i=1:(extendLen/2)
output_signal(i)=input_signal((extendLen/2)+2-i);
output_signal(length(output_signal)+1-i)=input_signal(length(input_signal)-(extendLen/2)-1+i);
end

for i=1:length(input_signal)
output_signal((extendLen/2)+i)=input_signal(i);
end
end

function [signal_korekce,prodlouzeni]=prodlouzeni_automaticke(signal,stupen_rozkladu,minimalni_prodlouzeni)
% Extended signal with certain length
% and the degree of decomposition required minimum extension

% Extended signal with certain length? 
% This technique is used for transient suppression effect
% when applying filter
% before discrete wavelet transform swt
% the signal length depends on the number of required bands decomposition.

% Prodlouzeni_automaticke function itself spo¨¨¨ªt about how the pattern signal needs
% at the beginning and end of the extensions due to a selected degree of decomposition
% But at the same time zaru¨¨¨ª that this extension will be shorter than the value of x
% APPLICATION This feature is especially In case where zapot?eb make large
% Quantity wavelet transformation signals of various lengths. the last input
% argument is optional if the minimum length is extended not choose then
% this value is set at 120 vzrok (60 at the beginning, 60 at the end).

%% Input:
% signal: the signal you want to extend, 
% stupen_rozkladu: the degree of decomposition,
% minimalni_prodlouzeni: the minimum length which requires extended
%%

if nargin~=3
  minimalni_prodlouzeni=120;  
end

soucasna_delka=length(signal);

pow = 2^stupen_rozkladu;
potrebna_delka=ceil((soucasna_delka+minimalni_prodlouzeni)/pow)*pow;
prodlouzeni=potrebna_delka-soucasna_delka;

if  rem(prodlouzeni,2)==0
    odeber=0;
else
    prodlouzeni=prodlouzeni+1;
    odeber=1;
end

signal_korekce=zeros(length(signal)+prodlouzeni,1);
for i=1:(prodlouzeni/2)
signal_korekce(i)=signal((prodlouzeni/2)+2-i);
signal_korekce(length(signal_korekce)+1-i)=signal(length(signal)-(prodlouzeni/2)-1+i);
end

for i=1:length(signal)
signal_korekce((prodlouzeni/2)+i)=signal(i);
end

signal_korekce=signal_korekce(1:(length(signal_korekce)-odeber));
end

function [std]=smerodatna_odchylka_sumu(waveletCoefficient)
% Calculating a robust estimate of the standard deviation of noise floating window
% Input: a vector of wavelet coefficients?

std=(median(abs(waveletCoefficient)))/0.6745; 
end

function [output_signal]=zkraceni(input_signal,truncateLen)
%% truncate cone signal
% cut off signal at the beginning and end
% eg. zkraceni(signal,24)
% signal will be cut 12 samples off at the beginning and the end
% Input value must be barrel number

output_signal=input_signal((truncateLen/2)+1:length(input_signal)-(truncateLen/2));
end

