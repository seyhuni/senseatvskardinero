clear; clc; close all;
fid = fopen('huseyinkardinero16122020.bin');
% fid = fopen('deney2.bin');
% fid = fopen('deney3.bin');
data = fread(fid, '*int16');
fclose(fid);
uzunluk=length(data);

j=2:12:uzunluk;
kanal2=data(j);

% Kardinero
d1=-kanal2;
d1=d1-mean(d1);
[c,l]=wavedec(d1,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
d1=double(d1)-a8;
d1=d1/max(d1);
blp= fir1(100,20/500,'low');
ecg1 = filter(blp,1,d1);

%ecg1 QRS 
t = 1:length(ecg1);  % aVR de 0.35 diðerinde 0.5
% [~,r_peak_pos1] = findpeaks(ecg1,'MinPeakHeight',0.35,...
%                                     'MinPeakDistance',200);
figure
hold on 
plot(t,ecg1)
load huseyinsenseat16122020.mat;

%ECG Sensor
d2=datas3;
d2=d2-mean(d2);
[c,l]=wavedec(d2,8,'db6');
a8 = wrcoef('a',c,l,'db6',8);
d2=double(d2)-a8;
d2=d2/max(d2);
blp= fir1(100,20/180,'low');
ecg2 = filter(blp,1,d2);

%ecg1 QRS 
tt = 1:length(ecg2);  % aVR de 0.35 diðerinde 0.5
% [~,r_peak_pos2] = findpeaks(ecg2,'MinPeakHeight',0.35,...
%                                     'MinPeakDistance',200);
figure
hold on 
plot(tt,ecg2)
% plot(r_peak_pos2,ecg2(r_peak_pos2),'rv','MarkerFaceColor','r')
% hv=ones(1,length(ecg2));

%ecg sensor points
%975 to 62002
%cardinero points
%3994 to 163953

% ecgSensorData = resample(ecg2,13,5);
figure
plot(ecg1)
title('Filtered Kardinero')
figure
plot(ecg2)
title('Filtered sensEAT')
figure
plot(datas3)
title('SensEAT')
figure
plot(-kanal2)
title('Kardinero')

senseat=ecg2(260:61608);
kardison=ecg1(2247:162616);
%260 to 619 
%2247 to 3188
% (3188-2247)/(619-260)
ecgSensorData = resample(senseat,941,359);
figure
plot(kardison)
hold on
plot(ecgSensorData,'r')

kardineroData=kardison(9007:end);
ecgSensorDatason=ecgSensorData(9162:end);
kardineroData=kardineroData(1:137700);
ecgSensorDatason=ecgSensorDatason(1:137733);

%ecgSensorData Characteristic Points
fs=1000;
%% Wavelet Transform
% Using Wavelet Transform to Decompose signal
[c,l] = wavedec(ecgSensorDatason,8,'db6');
%% Wavelet Reconstruction with Coefficients
% DESCRIPTIVE TEXT
for t=1:8
    D(:,t)=wrcoef('d',c,l,'db6',t);
end
a2=wrcoef('a',c,l,'db6',2);
%% R-Peak Detection
% Use selected coefficients to R Peak Detection
e1=D(:,3)+D(:,4)+D(:,5);
e2=(D(:,4).*(D(:,3)+D(:,5)))/2.^8;
R_Peak_Detect_Ecg=e1.*e2;
R_Peak_Detect_Ecg_Positive = zeros(1,length(R_Peak_Detect_Ecg));
for k=1:length(R_Peak_Detect_Ecg)   
if R_Peak_Detect_Ecg(k)>0
R_Peak_Detect_Ecg_Positive(k)=R_Peak_Detect_Ecg(k);
end
end
last_ecg=R_Peak_Detect_Ecg_Positive;
threshold=max(last_ecg);
threshold=threshold*0.01;
r_peak=0;
x=ceil(fs/9);

for k=x:(length(last_ecg)-x)
        gecici=last_ecg(k-x+1:k+x);
        if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>threshold)
        r_peak(k) = last_ecg(k);
    end
end
r_peak_pos=find(r_peak>0);
for j=1:length(r_peak_pos)-1
    if(abs(r_peak_pos(j)-r_peak_pos(j+1))<=53) 
        if(r_peak(r_peak_pos(j))>r_peak(r_peak_pos(j+1)))
        r_peak(r_peak_pos(j+1))=0;
        elseif(r_peak(r_peak_pos(j))<r_peak(r_peak_pos(j+1)))
            r_peak(r_peak_pos(j))=0;
        end
    end
        
end
t=find(r_peak>0);
r_peak_pos=t;

r_peak_last=zeros(1,length(last_ecg));
for t=1:length(r_peak_pos)
    mt3=ecgSensorDatason(r_peak_pos(t)-60:r_peak_pos(t)+60);
    mn3=max(mt3);
    r_peak_last(find(mt3==mn3)+r_peak_pos(t)-61)=mn3;
end
r_peak_pos_last=find(r_peak_last>0);

for i=1:length(r_peak_pos_last)-5
if((r_peak_pos_last(i+1)-r_peak_pos_last(i))<500) & (ecgSensorDatason(r_peak_pos_last(i+1))>ecgSensorDatason(r_peak_pos_last(i)))
        r_peak_pos_last(i)= [];
    else if(((r_peak_pos_last(i+1)-r_peak_pos_last(i))<500) & (ecgSensorDatason(r_peak_pos_last(i+1))<ecgSensorDatason(r_peak_pos_last(i))))
                r_peak_pos_last(i+1)= [];
        end
    end
end
%% Q & S Detection
hv=ones(1,length(ecgSensorDatason));
s_peak=sqrt(-1)*hv;
q_peak=s_peak;

for k=1:length(r_peak_pos_last)
    mt=ecgSensorDatason(r_peak_pos_last(k):r_peak_pos_last(k)+40);
    mn=min(mt);
    s_peak(find(mt==mn)+r_peak_pos_last(k)-1)=mn;
end
s_peak_pos_last=find(s_peak~=sqrt(-1));

for m=1:length(r_peak_pos_last)
    mt2=ecgSensorDatason(r_peak_pos_last(m)-50:r_peak_pos_last(m));
    mn2=min(mt2);
    q_peak(r_peak_pos_last(m)-50-1+find(mt2==mn2))= mn2;
end
q_peak_pos_last=find(q_peak~=sqrt(-1));

%% P and T Detection   ******* Burda kaldým Buraya kadar Süper 
e4=D(:,4)+D(:,5)+D(:,6)+D(:,7)+D(:,8);

t_peak=hv*sqrt(-1);
for t=1:length(s_peak_pos_last)-1
    mt6=e4(s_peak_pos_last(t):s_peak_pos_last(t)+250);
    mn6=max(mt6);
    t_peak(find(mt6==mn6)+s_peak_pos_last(t)-1)=mn6;
end
t_peak_pos=find(t_peak~=sqrt(-1));

t_peak_last=hv*sqrt(-1);

for t=1:length(t_peak_pos)
    mt9=ecgSensorDatason(t_peak_pos(t)-10:t_peak_pos(t)+10);
    mn9=max(mt9);
    in = find(mt9==mn9);
    t_peak_last(in(1)+t_peak_pos(t)-11)=mn9;
end

t_peak_pos_last = find(t_peak_last~=sqrt(-1));

% for x=1:length(t_peak_pos_last)
% if(r_peak_pos_last(x+1)-t_peak_pos_last(x)<(ceil(fs/bpm)*6))
%     t_peak_last(t_peak_pos_last(x))=0;
% mt13=ecgSensorDatason(s_peak_pos_last(x):s_peak_pos_last(x)+(ceil(fs/bpm)*8));
% mn13=max(mt13);
% in=find(mt13==mn13);
% t_peak_last(in(1)+s_peak_pos_last(x)-1)=mn13;
% end
% end
% t_peak_pos_last = find(t_peak_last~=sqrt(-1));

p_peak=hv*sqrt(-1);

for p=1:length(q_peak_pos_last)
    mt10=e4(q_peak_pos_last(p)- 120:q_peak_pos_last(p));
    mn10=max(mt10);
    p_peak(find(mt10==mn10)+q_peak_pos_last(p)-120-1)=mn10;
end

p_peak_pos=find(p_peak~=sqrt(-1));

p_peak_last=hv*sqrt(-1);

for p=1:length(p_peak_pos)
    mt11=ecgSensorDatason(p_peak_pos(p)- 10:p_peak_pos(p)+10);
    mn11=max(mt11);
    in = find(mt11==mn11);
    p_peak_last(in(1)+p_peak_pos(p)-10+1)=mn11;
end

p_peak_pos_last=find(p_peak_last~=sqrt(-1));
%%
%%P ve T Start and Final Points
% P Start
p_peak_start=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt15=ecgSensorDatason(p_peak_pos_last(i)-80:p_peak_pos_last(i));
    mn15=min(mt15);
    in = find(mt15==mn15);
    p_peak_start(in(1)+p_peak_pos_last(i)-80-1)=mn15;
end
p_peak_start_pos=find(p_peak_start~=sqrt(-1));

%P Final
p_peak_final=hv*sqrt(-1);
for i=1:length(p_peak_pos_last)
    mt16=ecgSensorDatason(p_peak_pos_last(i):p_peak_pos_last(i)+45);
    mn16=min(mt16);
    in = find(mt16==mn16);
    p_peak_final(in(1)+p_peak_pos_last(i)-1)=mn16;
end
p_peak_final_pos=find(p_peak_final~=sqrt(-1));

for i=1:length(p_peak_pos_last)-1
    if(q_peak_pos_last(i+1)<=p_peak_final_pos(i))
        q_peak_pos_last(i+1)=q_peak_pos_last(i+1)+3;
        p_peak_final_pos(i)=p_peak_final_pos(i)-3;
    end
end

% T Final
t_peak_final=hv*sqrt(-1);
for i=1:length(t_peak_pos_last)
    mt18=ecgSensorDatason(t_peak_pos_last(i):t_peak_pos_last(i)+110);
    mn18=min(mt18);
    in = find(mt18==mn18);
    t_peak_final(in(1)+t_peak_pos_last(i)-1)=mn18;
end
t_peak_final_pos=find(t_peak_final~=sqrt(-1));

    t_peak_start=hv*sqrt(-1);   
for i=1:length(t_peak_pos_last)
    mt17=ecgSensorDatason(t_peak_pos_last(i)-120:t_peak_pos_last(i));
    mn17=min(mt17);
    in = find(mt17==mn17);
    t_peak_start(in(1)+t_peak_pos_last(i)-121)=mn17;
end
t_peak_start_pos=find(t_peak_start~=sqrt(-1));  


%kardineroData Characteristic Points
fs=1000;
%% Wavelet Transform
% Using Wavelet Transform to Decompose signal
kardineroData=kardineroData';
[c1,l1] = wavedec(kardineroData,8,'db6');
%% Wavelet Reconstruction with Coefficients
% DESCRIPTIVE TEXT
for t=1:8
    D_(:,t)=wrcoef('d',c1,l1,'db6',t);
end
a2=wrcoef('a',c1,l1,'db6',2);
%% R-Peak Detection
% Use selected coefficients to R Peak Detection
e1=D_(:,3)+D_(:,4)+D_(:,5);
e2=(D_(:,4).*(D_(:,3)+D_(:,5)))/2.^8;
r_peakK_Detect_Ecg_=e1.*e2;
r_peakK_Detect_Ecg__Positive = zeros(1,length(r_peakK_Detect_Ecg_));
for k=1:length(r_peakK_Detect_Ecg_)   
if r_peakK_Detect_Ecg_(k)>0
r_peakK_Detect_Ecg__Positive(k)=r_peakK_Detect_Ecg_(k);
end
end
last_ecg_=r_peakK_Detect_Ecg__Positive;
threshold=max(last_ecg_);
threshold=threshold*0.01;
r_peakK=0;
x=ceil(fs/9);

for k=x:(length(last_ecg_)-x)
        gecici=last_ecg_(k-x+1:k+x);
        if(gecici(ceil(length(gecici)/2))==max(gecici) & gecici(ceil(length(gecici)/2))>threshold)
        r_peakK(k) = last_ecg_(k);
    end
end
r_peakK_pos=find(r_peakK>0);
for j=1:length(r_peakK_pos)-1
    if(abs(r_peakK_pos(j)-r_peakK_pos(j+1))<=53) 
        if(r_peakK(r_peakK_pos(j))>r_peakK(r_peakK_pos(j+1)))
        r_peakK(r_peakK_pos(j+1))=0;
        elseif(r_peakK(r_peakK_pos(j))<r_peakK(r_peakK_pos(j+1)))
            r_peakK(r_peakK_pos(j))=0;
        end
    end
        
end
t=find(r_peakK>0);
r_peakK_pos=t;

r_peakK_last=zeros(1,length(last_ecg_));
for t=1:length(r_peakK_pos)
    mt3=kardineroData(r_peakK_pos(t)-60:r_peakK_pos(t)+60);
    mn3=max(mt3);
    r_peakK_last(find(mt3==mn3)+r_peakK_pos(t)-61)=mn3;
end
r_peakK_pos_last=find(r_peakK_last>0);

for i=1:length(r_peakK_pos_last)-5
if((r_peakK_pos_last(i+1)-r_peakK_pos_last(i))<500) & (kardineroData(r_peakK_pos_last(i+1))>kardineroData(r_peakK_pos_last(i)))
        r_peakK_pos_last(i)= [];
    else if(((r_peakK_pos_last(i+1)-r_peakK_pos_last(i))<500) & (kardineroData(r_peakK_pos_last(i+1))<kardineroData(r_peakK_pos_last(i))))
                r_peakK_pos_last(i+1)= [];
        end
    end
end
%% Q & S Detection
hv=ones(1,length(kardineroData));
s_peak_=sqrt(-1)*hv;
q_peak_=s_peak_;

for k=1:length(r_peakK_pos_last)
    mt=kardineroData(r_peakK_pos_last(k):r_peakK_pos_last(k)+40);
    mn=min(mt);
    s_peak_(find(mt==mn)+r_peakK_pos_last(k)-1)=mn;
end
s_peak__pos_last=find(s_peak_~=sqrt(-1));

for m=1:length(r_peakK_pos_last)
    mt2=kardineroData(r_peakK_pos_last(m)-50:r_peakK_pos_last(m));
    mn2=min(mt2);
    q_peak_(r_peakK_pos_last(m)-50-1+find(mt2==mn2))= mn2;
end
q_peak__pos_last=find(q_peak_~=sqrt(-1));

%% P and T Detection   ******* Burda kaldým Buraya kadar Süper 
e4=D_(:,4)+D_(:,5)+D_(:,6)+D_(:,7)+D_(:,8);

t_peak_=hv*sqrt(-1);
for t=1:length(s_peak__pos_last)-1
    mt6=e4(s_peak__pos_last(t):s_peak__pos_last(t)+250);
    mn6=max(mt6);
    t_peak_(find(mt6==mn6)+s_peak__pos_last(t)-1)=mn6;
end
t_peak__pos=find(t_peak_~=sqrt(-1));

t_peak__last=hv*sqrt(-1);

for t=1:length(t_peak__pos)
    mt9=kardineroData(t_peak__pos(t)-10:t_peak__pos(t)+10);
    mn9=max(mt9);
    in = find(mt9==mn9);
    t_peak__last(in(1)+t_peak__pos(t)-11)=mn9;
end

t_peak__pos_last = find(t_peak__last~=sqrt(-1));

% for x=1:length(t_peak__pos_last)
% if(r_peakK_pos_last(x+1)-t_peak__pos_last(x)<(ceil(fs/bpm)*6))
%     t_peak__last(t_peak__pos_last(x))=0;
% mt13=kardineroData(s_peak__pos_last(x):s_peak__pos_last(x)+(ceil(fs/bpm)*8));
% mn13=max(mt13);
% in=find(mt13==mn13);
% t_peak__last(in(1)+s_peak__pos_last(x)-1)=mn13;
% end
% end
% t_peak__pos_last = find(t_peak__last~=sqrt(-1));

p_peak_=hv*sqrt(-1);

for p=1:length(q_peak__pos_last)
    mt10=e4(q_peak__pos_last(p)- 120:q_peak__pos_last(p));
    mn10=max(mt10);
    p_peak_(find(mt10==mn10)+q_peak__pos_last(p)-120-1)=mn10;
end

p_peak__pos=find(p_peak_~=sqrt(-1));

p_peak__last=hv*sqrt(-1);

for p=1:length(p_peak__pos)
    mt11=kardineroData(p_peak__pos(p)- 10:p_peak__pos(p)+10);
    mn11=max(mt11);
    in = find(mt11==mn11);
    p_peak__last(in(1)+p_peak__pos(p)-10+1)=mn11;
end

p_peak__pos_last=find(p_peak__last~=sqrt(-1));
%%
%%P ve T Start and Final Points
% P Start
p_peak__start=hv*sqrt(-1);
for i=1:length(p_peak__pos_last)
    mt15=kardineroData(p_peak__pos_last(i)-80:p_peak__pos_last(i));
    mn15=min(mt15);
    in = find(mt15==mn15);
    p_peak__start(in(1)+p_peak__pos_last(i)-80-1)=mn15;
end
p_peak__start_pos=find(p_peak__start~=sqrt(-1));

%P Final
p_peak__final=hv*sqrt(-1);
for i=1:length(p_peak__pos_last)
    mt16=kardineroData(p_peak__pos_last(i):p_peak__pos_last(i)+45);
    mn16=min(mt16);
    in = find(mt16==mn16);
    p_peak__final(in(1)+p_peak__pos_last(i)-1)=mn16;
end
p_peak__final_pos=find(p_peak__final~=sqrt(-1));

for i=1:length(p_peak__pos_last)-1
    if(q_peak__pos_last(i+1)<=p_peak__final_pos(i))
        q_peak__pos_last(i+1)=q_peak__pos_last(i+1)+3;
        p_peak__final_pos(i)=p_peak__final_pos(i)-3;
    end
end

% T Final
t_peak__final=hv*sqrt(-1);
for i=1:length(t_peak__pos_last)
    mt18=kardineroData(t_peak__pos_last(i):t_peak__pos_last(i)+110);
    mn18=min(mt18);
    in = find(mt18==mn18);
    t_peak__final(in(1)+t_peak__pos_last(i)-1)=mn18;
end
t_peak__final_pos=find(t_peak__final~=sqrt(-1));

    t_peak__start=hv*sqrt(-1);   
for i=1:length(t_peak__pos_last)
    mt17=kardineroData(t_peak__pos_last(i)-120:t_peak__pos_last(i));
    mn17=min(mt17);
    in = find(mt17==mn17);
    t_peak__start(in(1)+t_peak__pos_last(i)-121)=mn17;
end
t_peak__start_pos=find(t_peak__start~=sqrt(-1));  


%% RESULTS & Findings
figure;
hold on;
plot(ecgSensorDatason)
plot(r_peak_pos_last,ecgSensorDatason(r_peak_pos_last),'r+','MarkerFaceColor','r')
plot(s_peak_pos_last,ecgSensorDatason(s_peak_pos_last),'r*','MarkerFaceColor','r')
plot(q_peak_pos_last,ecgSensorDatason(q_peak_pos_last),'r.','MarkerFaceColor','r')
plot(p_peak_start_pos,ecgSensorDatason(p_peak_start_pos),'k+','MarkerFaceColor','r')
plot(p_peak_pos_last,ecgSensorDatason(p_peak_pos_last),'k*','MarkerFaceColor','r')
plot(p_peak_final_pos,ecgSensorDatason(p_peak_final_pos),'k.','MarkerFaceColor','r')
plot(t_peak_start_pos,ecgSensorDatason(t_peak_start_pos),'b+','MarkerFaceColor','r')
plot(t_peak_pos_last,ecgSensorDatason(t_peak_pos_last),'b*','MarkerFaceColor','r')
plot(t_peak_final_pos,ecgSensorDatason(t_peak_final_pos),'b.','MarkerFaceColor','r')
title('Senseat Findings')
figure;
hold on;
plot(kardineroData)
plot(r_peakK_pos_last,kardineroData(r_peakK_pos_last),'r+','MarkerFaceColor','r')
plot(s_peak__pos_last,kardineroData(s_peak__pos_last),'r*','MarkerFaceColor','r')
plot(q_peak__pos_last,kardineroData(q_peak__pos_last),'r.','MarkerFaceColor','r')
plot(p_peak__start_pos,kardineroData(p_peak__start_pos),'k+','MarkerFaceColor','r')
plot(p_peak__pos_last,kardineroData(p_peak__pos_last),'k*','MarkerFaceColor','r')
plot(p_peak__final_pos,kardineroData(p_peak__final_pos),'k.','MarkerFaceColor','r')
plot(t_peak__start_pos,kardineroData(t_peak__start_pos),'b+','MarkerFaceColor','r')
plot(t_peak__pos_last,kardineroData(t_peak__pos_last),'b*','MarkerFaceColor','r')
plot(t_peak__final_pos,kardineroData(t_peak__final_pos),'b.','MarkerFaceColor','r')
title('Kardinero Findings')
% HRVs

for i=1:length(r_peak_pos_last)-1
    HRVsensEAT(i)=r_peak_pos_last(i+1)-r_peak_pos_last(i);
end

for i=1:length(r_peakK_pos_last)-1
    HRVkardinero(i)=r_peakK_pos_last(i+1)-r_peakK_pos_last(i);
end
HRVsensEAT=HRVsensEAT(1:127);
HRVkardinero=HRVkardinero(1:127);
% HRV is calculated as the standard deviation of all of the 
% RR intervals (the distance between each “R” of the QRS complex).

hrv_senseat=std(HRVsensEAT)
hrv_kardinero=std(HRVkardinero)
%% ST segment
st_interval_senseat= zeros(1,length(t_peak_start_pos));
% for i =1:length(st_interval_senseat)
%     st_interval_senseat(i)=t_peak_final_pos(i)-s_peak_pos_last(i); 
%     %t bitiþ - s
% end
% mean_st_interval_senseat= mean(st_interval_senseat(1:127))
% 
st_interval_kardinero= zeros(1,length(t_peak__start_pos));
% for i =1:length(st_interval_kardinero)
%     st_interval_kardinero(i)=t_peak__final_pos(i)-s_peak__pos_last(i); 
%     %t bitiþ - s
% end
% mean_st_interval_kardinero= mean(st_interval_kardinero(1:127))

st_segment_senseat=zeros(1,length(t_peak_pos_last));
for i =1:length(st_interval_senseat)
    st_segment_senseat(i)=t_peak_start_pos(i)-s_peak_pos_last(i);
    % t baþlangýç - s
end
mean_st_segment_senseat=mean(st_segment_senseat(1:127))

st_segment_kardinero=zeros(1,length(t_peak__pos_last));
for i =1:length(st_interval_kardinero)
    st_segment_kardinero(i)=t_peak__start_pos(i)-s_peak__pos_last(i);
    % t baþlangýç - s
end
mean_st_segment_kardinero=mean(st_segment_kardinero(1:127))

%% QRS Complex intervals

qrs_int_senseat=zeros(1,length(q_peak_pos_last));
for i=1:length(qrs_int_senseat)
    qrs_int_senseat(i)=s_peak_pos_last(i)-q_peak_pos_last(i);
end
qrs_interval_senseat=mean(qrs_int_senseat(1:127))

qrs_int_kardinero=zeros(1,length(q_peak__pos_last));
for i=1:length(qrs_int_kardinero)
    qrs_int_kardinero(i)=s_peak__pos_last(i)-q_peak__pos_last(i);
end
qrs_interval_kardinero=mean(qrs_int_kardinero(1:127))

%% RR interval
rr_int_senseat=zeros(1,length(r_peak_pos_last)-1);

for i=1:length(r_peak_pos_last)-1
    rr_int_senseat(i)=(r_peak_pos_last(i+1)-r_peak_pos_last(i));
end
mean_rr_senseat=mean(rr_int_senseat(1:127))

rr_int_kardinero=zeros(1,length(r_peak_pos_last)-1);

for i=1:length(r_peakK_pos_last)-1
    rr_int_kardinero(i)=(r_peakK_pos_last(i+1)-r_peakK_pos_last(i));
end
mean_rr_kardinero=mean(rr_int_kardinero(1:127))

