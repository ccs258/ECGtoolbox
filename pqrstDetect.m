function [P,Q,R,S,T,Pon,Poff,QRSon,QRSoff,ST,Toff]=pqrstDetect(ecg,fs)
%% Description
% detect ecg signal feature
% including P wave, QRS wave, T wave
%% Inputs
% ecg : raw ecg vector
% fs : sampling frequency

%% Outputs
%their estimated index ratio of the sample number within one
%heart impulse are shown below
%
%% usage
% for example after loading the the ecg mat files in matlab call the
% function as below ;
% [P,Q,R,S,T,Pon,Poff,QRSon,QRSoff,ST,Toff]=featureDetect(ecg,500)


%% Author : Lisha Chen    
% contact : lishachen00@gmail.com , lisha_chen@hust.edu.cn
% Dont forget to reference if you found this script usefull
%%

    if nargin <2
       fs = 360; %default Sampling frequency
    end    

%% parameters

Ponr=0.187;
Poffr=0.279;
Pr=0.233;

QRSonr=0.308; 
QRSoffr=0.456; 
Qr=0.349;
Rr=0.382;
Sr=0.422;

STr=0.534; 
Tr=0.619; 
Toffr=0.709;

%% ECG signal preprocess-remove noise
ecg = ecg (:); % make sure its a vector

maxTT=sort(ecg,'descend')-mean(ecg);
maxTT=mean(maxTT(1:20));
minTT=sort(ecg,'ascend')-mean(ecg);
minTT=mean(minTT(1:20));

if abs(minTT)>abs(maxTT)
    ecgRaw=-ecg;
    maxT=abs(minTT);
else ecgRaw=ecg;
    maxT=abs(maxTT);
end

%Noise cancelation(Filtering)
f1=0.5; %cuttoff low frequency to get rid of baseline wander
f2=45; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg000 = filtfilt(a,b,ecg);

h=fspecial('average',[10 1]);
ecg00=conv(ecg,h,'same');

wtsig1=cwt(ecg,6,'mexh'); 
ecg01=wtsig1; 

f1=0.5; %cuttoff low frequency to get rid of baseline wander
f2=20; %cuttoff frequency to discard high frequency noise
Wn=[f1 f2]*2/fs; % cutt off based on fs
N = 3; % order of 3 less processing
[a,b] = butter(N,Wn); %bandpass filtering
ecg1 = filtfilt(a,b,ecg01);

diffECG00=diff(ecg00);
%% detect R
R=[];
for i=2:length(ecg1)-1
if ecg1(i-1)<=ecg1(i)&&ecg1(i+1)<=ecg1(i) ...
    &&(ecg1(i)-mean(ecg1))>0.7*maxT % max extreme points
    R=[R,i];
end
end

lenR=length(R);
while i<=lenR
      if (R(i)-R(i-1))<0.4*fs 
          if ecg1(R(i))>ecg1(R(i-1)) 
              R(i-1)=[]; 
          else 
              R(i)=[]; 
          end
          lenR=length(R);
          i=i-1; 
      end
      i=i+1; 
end

% find max based on ecgRaw
for i=1:lenR
    tR=R(i);
        k=(tR-5):(tR+5); 
        [~,b]=max(ecgRaw(k));
        R(i)=tR-6+b;  
end
R=R(:);

%% initialize and calculate parameters

S = zeros(lenR,1);%index of S wave
T = zeros(lenR,1);%index of T wave
Q = zeros(lenR,1);%index of Q wave
Pon=zeros(lenR,1);%begin of P wave
P=zeros(lenR,1);%index of P wave
Poff=zeros(lenR,1);%end of P wave
QRSon=zeros(lenR,1);%begin of QRS wave
QRSoff=zeros(lenR,1);%end of QRS wave
ST=zeros(lenR,1);%corner of T wave
Toff=zeros(lenR,1);%end of T wave

wtsig2=cwt(ecgRaw,8,'mexh'); 
Tsample=length(ecg)/lenR;

trans=mean(((R./Tsample)'-(0:lenR-1)))-Rr;

Ponr=Ponr+trans;
Poffr=Poffr+trans;
Pr=Pr+trans;

QRSonr=QRSonr+trans;
QRSoffr=QRSoffr+trans;
Qr=Qr+trans;

Sr=Sr+trans;

STr=STr+trans;
Tr=Tr+trans;
Toffr=Toffr+trans;

%% detect Q S QRSon QRSoff
%---------------Q---------------------------------------------------------
for i=1:lenR 
    for j=ceil(R(i)-1:-1:(R(i)-0.07*Tsample)) 
            if wtsig2(j)<=wtsig2(j-1)&&wtsig2(j)<=wtsig2(j+1) %min extreme points
                tempQvalue=ceil(j-0.04*Tsample); %detection window starting index
            break
            else if j==ceil(R(i)-0.07*Tsample)
                    tempQvalue=ceil((Qr+i-1)*Tsample);
                 end
            end
    end
    x1=tempQvalue; %(x1,y1)--Q wave candidates coordinate
    y1=ecgRaw(tempQvalue); 
    x2=R(i); %(x2,y2)--R wave coordinate
    y2=ecgRaw(R(i)); 
    a0=(y2-y1)/(x2-x1); 
    b0=-1; 
    c0=-a0*x1+y1; 
    dist=[]; 
    for k=tempQvalue:R(i) 
        tempdist=(abs(a0*k+b0*ecgRaw(k)+c0))/sqrt(a0^2+b0^2); 
        dist=[dist;tempdist]; 
    end   %distance between 2 points
    [~,b]=max(dist); %max of distance
    tempQvalue=tempQvalue+b-1;
    l=(tempQvalue-5):R(i); 
    [~,d]=min(ecgRaw(l)); 
    tempQvalue=tempQvalue-6+d; %correct coordinates based on ecgRaw
    Q(i)=tempQvalue; %Q--Q wave index           
end  

%---------------S---------------------------------------------------------
for i=1:lenR
    for j=ceil(R(i)+1:1:(R(i)+0.07*Tsample)) 
           if (wtsig2(j)<=wtsig2(j-1))&&(wtsig2(j)<=wtsig2(j+1)) %min extreme points
          tempSvalue=ceil(j+0.04*Tsample); %detection window starting index
               break
               else if j==ceil(R(i)+0.07*Tsample)
                    tempSvalue=ceil((Sr+i-1)*Tsample);
                    end
           end
    end
    
    x1=tempSvalue; %(x1,y1)--S wave candidates coordinate
    y1=ecgRaw(tempSvalue); 
    x2=R(i); %(x2,y2)--R wave coordinate
    y2=ecgRaw(R(i)); 
    a0=(y2-y1)/(x2-x1); 
    b0=-1; 
    c0=-a0*x1+y1; %line ax+by+c=0 
    dist=[]; 
    for k=R(i):tempSvalue 
        tempdist=(abs(a0*k+b0*ecgRaw(k)+c0))/sqrt(a0^2+b0^2); 
        dist=[dist;tempdist]; 
    end  %distance between 2 points           
    [~,b]=max(dist);  
    tempSvalue=R(i)+b-1; 
   l=R(i):(tempSvalue+10); 
   [~,d]=min(ecgRaw(l)); 
   tempSvalue=R(i)+d-1; %correct coordinates based on ecgRaw
   S(i)=tempSvalue; %S--S wave index      
end 

%---------------QRSon-----------------------------------------------------
for i=1:lenR 
    for j=ceil(Q(i)-1:-1:(Q(i)-0.04*Tsample)) 
            if wtsig2(j)>=wtsig2(j-1)&&wtsig2(j)>=wtsig2(j+1) %max extreme points
                QRSon(i)=j;
            break
               else if j==ceil(Q(i)-0.04*Tsample)
                    QRSon(i)=ceil((QRSonr+i-1)*Tsample);
                    end
            end
    end
end

for i=1:lenR
    tQRSon=QRSon(i);
        for k=ceil((tQRSon-0.02*Tsample):1:min((tQRSon+0.04*Tsample),Q(i)-1))
        if ecg00(k-1)<=ecg00(k)&&ecg00(k+1)<=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0
        QRSon(i)=k;        
        end
        end
end

%---------------QRSoff----------------------------------------------------
for i=1:lenR 
    for j=ceil(S(i)+1:1:(S(i)+0.06*Tsample)) 
            if wtsig2(j)>=wtsig2(j-1)&&wtsig2(j)>=wtsig2(j+1) %max extreme points
                QRSoff(i)=j;
            break
               else if j==ceil(S(i)+0.06*Tsample)
                    QRSoff(i)=ceil((QRSoffr+i-1)*Tsample);
                    end
            end
    end
end

for i=1:lenR
    tQRSoff=QRSoff(i);
        for k=ceil(max((tQRSoff+0.02*Tsample),S(i)+1):-1:(tQRSoff-0.06*Tsample)) 
        if ecg00(k-1)<=ecg00(k)&&ecg00(k+1)<=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0
        QRSoff(i)=k;        
        end
        end
end

%% detect P Pon Poff
%---------------Poff------------------------------------------------------
for i=1:lenR 
    for j=ceil(QRSon(i)-1:-1:(QRSon(i)-0.05*Tsample)) 
            if wtsig2(j)<=wtsig2(j-1)&&wtsig2(j)<=wtsig2(j+1) %min extreme points
                Poff(i)=j;
            break
            else if j==ceil(QRSon(i)-0.04*Tsample)
                    Poff(i)=ceil((Poffr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tPoff=Poff(i);
    tempPoff=Poff(i);
        for k=ceil(min((tPoff+0.03*Tsample),QRSon(i)-1):-1:(tPoff-0.03*Tsample)) 
        if ecg00(k-1)>=ecg00(k)&&ecg00(k+1)>=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0 ...
                &&(ecg00(k-5)-ecg00(k))/5<=-maxT/500
            tempPoff=k;
        break       
        end        
        [~,b]=min(ecgRaw((tempPoff-0.02*Tsample):min((tempPoff+0.07*Tsample),QRSon(i)-1))); 
        Poff(i)=ceil(tempPoff-0.02*Tsample+b-1);
        end       
end

%---------------P---------------------------------------------------------
for i=1:lenR 
    for j=ceil(Poff(i)-0.02*Tsample:-1:(Poff(i)-0.1*Tsample)) 
            if wtsig2(j)>=wtsig2(j-1)&&wtsig2(j)>=wtsig2(j+1) %max extreme points
                P(i)=j;
            break
            else if j==ceil(Poff(i)-0.1*Tsample)
                    P(i)=ceil((Pr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tP=P(i);
    tempP=P(i);
        for k=ceil((tP-0.05*Tsample):1:min((tP+0.05*Tsample),Poff(i)-1))
        if ecg00(k-1)<=ecg00(k)&&ecg00(k+1)<=ecg00(k) ...
             &&(ecg00(i)-mean(ecg00))>0.1*maxT                      
        tempP=k;
        break
        end
        [~,b]=max(ecgRaw((tempP-0.1*Tsample):min((tempP+0.1*Tsample),Poff(i)-1)));
        P(i)=ceil(tempP-0.1*Tsample+b-1);
        end
end

%---------------Pon-------------------------------------------------------
for i=1:lenR 
    for j=ceil(P(i)-0.02*Tsample:-1:(P(i)-0.1*Tsample)) 
            if wtsig2(j)<=wtsig2(j-1)&&wtsig2(j)<=wtsig2(j+1) %min extreme points
                Pon(i)=j;
            break
            else if j==ceil(P(i)-0.1*Tsample)
                    Pon(i)=ceil((Ponr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tPon=Pon(i);
        for k=ceil((tPon-0.03*Tsample):min((tPon+0.03*Tsample),P(i)-0.02*Tsample))
        if ecg00(k-1)>=ecg00(k)&&ecg00(k+1)>=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0 ...
                &&diffECG00(k+1)>=0&&abs(diffECG00(k-1))<1           
        Pon(i)=k;        
        end
        end
end

%% detect T ST Toff
%---------------ST--------------------------------------------------------
for i=1:lenR 
    for j=ceil(QRSoff(i)+1:1:(QRSoff(i)+0.15*Tsample)) 
            if wtsig2(j)<=wtsig2(j-1)&&wtsig2(j)<=wtsig2(j+1) 
                ST(i)=j; % min extreme points
            break
            else if j==ceil(QRSoff(i)+0.15*Tsample)
                    ST(i)=ceil((STr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tST=ST(i);
        for k=ceil(max((tST-0.03*Tsample),QRSoff(i)+1):(tST+0.03*Tsample)) 
        if ecg00(k-1)>=ecg00(k)&&ecg00(k+1)>=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0 ...
                &&diffECG00(k+1)>maxT/500            
        ST(i)=k;        
        end
        end
end

%---------------T---------------------------------------------------------
for i=1:lenR 
    for j=ceil(ST(i)+0.02*Tsample:1:(ST(i)+0.15*Tsample)) 
            if wtsig2(j)>=wtsig2(j-1)&&wtsig2(j)>=wtsig2(j+1) ...
                    && j>=0.55*Tsample+Tsample*(i-1)&&j<=0.65*Tsample+Tsample*(i-1)
                T(i)=j;
            break % max extreme points
            else if j==ceil(ST(i)+0.15*Tsample)
                    T(i)=ceil((Tr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tT=T(i);
    tempT=T(i);
        for k=ceil(max((tT-0.03*Tsample),ST(i)+1):(tT+0.03*Tsample)) 
        if ecg00(k-1)<=ecg00(k)&&ecg00(k+1)<=ecg00(k) ...
             &&(ecg00(i)-mean(ecg00))>0.2*maxT              
        tempT=k;
        break
        end
        [~,b]=max(ecgRaw(max((tempT-0.04*Tsample),ST(i)+1):(tempT+0.04*Tsample))); 
        T(i)=ceil(max((tempT-0.04*Tsample),ST(i)+1)+b-1);
        end
end

%---------------Toff------------------------------------------------------
for i=1:lenR 
    for j=ceil(T(i)+0.05*Tsample:1:(T(i)+0.15*Tsample)) 
            if wtsig2(j)<=wtsig2(j-1)&&wtsig2(j)<=wtsig2(j+1) 
                Toff(i)=j;
            break  % min extreme points
            else if j==ceil(T(i)+0.15*Tsample)
                    Toff(i)=ceil((Toffr+i-1)*Tsample);
                end
            end
    end
end

for i=1:lenR
    tToff=Toff(i);
    tempToff=Toff(i);
        for k=ceil(max((tToff-0.03*Tsample),T(i)+0.05*Tsample):(tToff+0.03*Tsample)) 
        if ecg00(k-1)>=ecg00(k)&&ecg00(k+1)>=ecg00(k) ...
                &&((ecg00(k-1)-ecg00(k))^2+(ecg00(k+1)-ecg00(k))^2)~=0
        tempToff=k;
        break
        end
        [~,b]=min(ecgRaw(max((tempToff-0.04*Tsample),T(i)+0.05*Tsample):(tempToff+0.04*Tsample))); 
        Toff(i)=ceil(max((tempToff-0.04*Tsample),T(i)+0.05*Tsample)+b-1);
        end
end

%% plot

figure,
plot(ecgRaw);
hold on
plot(Toff,ecgRaw(Toff),'r*');
plot(ST,ecgRaw(ST),'r*');
plot(T,ecgRaw(T),'r*');

plot(Q,ecgRaw(Q),'r*');
plot(R,ecgRaw(R),'r*');
plot(S,ecgRaw(S),'r*');
plot(QRSoff,ecgRaw(QRSoff),'r*');
plot(QRSon,ecgRaw(QRSon),'r*');

plot(Poff,ecgRaw(Poff),'r*');
plot(Pon,ecgRaw(Pon),'r*');
plot(P,ecgRaw(P),'r*');
hold off
