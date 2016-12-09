function R=rDetect(ecg,fs)
%% Description
% detect ecg R wave
% This function is better than the one shown on 
% http://www.librow.com/cases/case-2
% tested on EKG1 to EKG5
%% Inputs
% ecg : raw ecg vector
% fs : sampling frequency
%% Outputs
% index of R wave
%% usage
% for example after loading the the ecg mat files in matlab call the
% function as below ;
% R=rDetect(EKG1,250)
%% Author : Lisha Chen    
% contact : lishachen00@gmail.com , lisha_chen@hust.edu.cn
% Dont forget to reference if you found this script usefull
%%
    if nargin <2
       fs = 360; %default Sampling frequency
    end    

%% ECG signal preprocess-remove noise
ecg = ecg (:); % make sure its a vector

%Noise cancelation(Filtering)

f1=0.5; %cuttoff low frequency to get rid of baseline wander
f2=45; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecgButter = filtfilt(a,b,ecg);

ecg00=conv(ecg,fspecial('average',[10 1]),'same');
ecg01=cwt(ecg00,6,'mexh'); 

f1=0.5; % cuttoff low frequency to get rid of baseline wander
f2=20; % cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cut off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); % bandpass filtering
ecg1 = filtfilt(a,b,ecg01);

maxTT=sort(ecgButter,'descend')-mean(ecgButter);
maxTT=mean(maxTT(1:ceil(length(ecg)/fs)));
minTT=sort(ecgButter,'ascend')-mean(ecgButter);
minTT=mean(minTT(1:ceil(length(ecg)/fs)));

if abs(minTT)>abs(maxTT)
    ecgRaw=-ecg;
    maxT=abs(minTT);
else ecgRaw=ecg;
    maxT=abs(maxTT);
end

%% detect R
R=[];
for i=2:length(ecg1)-1
if ecg1(i-1)<=ecg1(i)&&ecg1(i+1)<=ecg1(i) ...
    &&(ecgButter(i)-mean(ecgButter))>0.3*maxT % max extreme points
    R=[R,i];
end
end

lenR=length(R);
% remove too near-----------------------------
i=2;
while i<=lenR
      if (R(i)-R(i-1))<0.3*fs 
          if ecgButter(R(i))>ecgButter(R(i-1)) 
              R(i-1)=[]; 
          else 
              R(i)=[]; 
          end
          lenR=length(R);
          i=i-1; 
      end
      i=i+1; 
end

% remove too near-----------------------------
i=2;
while i<=lenR
      if (R(i)-R(i-1))<0.5*mean(R(2:end)-R(1:end-1)) 
          if ecgButter(R(i))>ecgButter(R(i-1)) 
              R(i-1)=[]; 
          else 
              R(i)=[]; 
          end
          lenR=length(R);
          i=i-1; 
      end
      i=i+1; 
end

% find the right R wave based on local max--------
tTemp=mean(R(2:end)-R(1:end-1));
for i=2:length(R)-1
 if (ecgButter(R(i))-mean(ecgButter))<0.7*maxT
     [num,id]=max(ecgButter(max([R(i)-floor(tTemp/2),1]):...
         min([length(ecg),R(i)+floor(tTemp/2)])));
     if num>ecgButter(R(i))&&(R(i)-floor(tTemp/2)+id-1)>R(i-1)&&...
             (R(i)-floor(tTemp/2)+id-1)<R(i+1)
         R(i)=R(i)-floor(tTemp/2)+id-1;
     end
 end       
end

% remove too near-----------------------------
i=2;
while i<=lenR
      if (R(i)-R(i-1))<0.5*mean(R(2:end)-R(1:end-1)) 
          if ecgButter(R(i))>ecgButter(R(i-1)) 
              R(i-1)=[]; 
          else 
              R(i)=[]; 
          end
          lenR=length(R);
          i=i-1; 
      end
      i=i+1; 
end

% add missing R wave------------------------------
tTemp=mean(R(2:end)-R(1:end-1));
maxT=mean(ecgButter(R)-mean(ecgButter));
for i=2:length(R)
 if (R(i)-R(i-1))>1.5*tTemp
     [num,id]=max(ecgButter(R(i-1)+floor(tTemp/2):R(i)-floor(tTemp/2)));
     if (num-mean(ecgButter))>0.7*maxT&&~isempty(id)
         R=[R(1:i-1),R(i-1)+floor(tTemp/2)+id-1,R(i:end)];
     end
 end       
end
    
% find max based on ecgRaw------------------
for i=1:length(R)
    %k=max([R(i)-floor(tTemp/2),1]):min([R(i)+floor(tTemp/2),length(ecg)]); 
    k=max([R(i)-5,1]):min([R(i)+5,length(ecg)]); 
    [~,b]=max(ecgRaw(k));
    R(i)=k(1)+b-1;  
end
R=R(:);
% plot---------------
% figure,
% plot(ecg);
% hold on
% plot(R,ecg(R),'r*');
end

