clc;
clear all;
close all;
 
JsfHead=gFJsfReadHeader('1000.jsf',1);
[Head,Data]=gFJsfRead0080(JsfHead,0,0);
Data1=abs(Data);
 
 i=1:999
 t=i*40*10.^(-6); % sampling time 40 us
%  figure(1);       %plot of envelop data(absolute of analytical data) 
%  plot(t,Data1(:,1));
%  hold on;
%  plot(t,Data1(:,2));
%  hold on;
%  plot(t,Data1(:,3));
%  hold on;
%  plot(t,Data1(:,5));
%  hold on;
%  plot(t,Data1(:,6));
%  hold on;
%  plot(t,Data1(:,7));
%  title('Envelope data');
%  xlabel('time in sec');
%  ylabel('amplitude');
 
 d=(1500*t)/2     % depth
% figure(2);       %plot depth vs the acoustic data
% plot(d,Data1(:,2160));
%  hold on;
%  plot(d,Data1(:,2));
%  hold on;
%  plot(d,Data1(:,3));
%  hold on;
%  plot(d,Data1(:,4));
%  title('Envelope data');
%  xlabel('Depth m');
%  ylabel('amplitude');
 fs=25000;
 
Data21=real(Data);
Data2=Data21./4;  % scaling factor of 4
 
 
 
 
[R,C]=size(Data2);
t=0:1/fs:(R-1)/fs;
 figure(40);
 plot(t,Data2(:,20));
%S=hilbert(Data);
%%%%%%%%%%%%
d=designfilt('bandpassiir','FilterOrder',20,'HalfPowerFrequency1',2000,'HalfPowerFrequency2',10000,'SampleRate',25000);
%[b a]=butter(10,[0.16 0.8],'bandpass');
 S=filtfilt(d,Data2(:,1:2162)); 
  figure(41);
  plot(t,S(:,20));
[R,C]=size(S);
 
t=0:1/fs:(R-1)/fs;
 hx = hilbert(S(:,2160));     % change here                                
 inst_amp = abs(hx);
 f_instphase= diff(unwrap(angle(hx)))/((1/fs)*2*pi);
  figure(3);
%plot(d,f1);
 d=(t(1:length(f_instphase))*750);
 d=d';
  plot((t(1:length(f_instphase))*750),f_instphase); %hilbert transform
  xlabel('Depth in m');
%  ylabel('frequency in Hz');
%  title('Instantaneous frequency');
%  
%%%%%%%%original chirp signal generation
 
%  fs=25000; %sampling frequency
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  f0=2000;% starting frequency of the chirp
%  f1=10000; %frequency of the chirp at t1=1 second
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  x = chirp([0:1/fs:0.02],2000,0.02,10000);
%  M = length(x);
%  n=(-(M-1)/2:(M-1)/2)';
%  w = exp(-n.*n./(2*sigma.*sigma));
%  xw = w(:) .* x(:);
%  t=0:1/fs:0.02;
%  figure(4);
%  subplot(121);
%  xw=20000*xw;
%  plot(t,(xw));
% [pxx,f] = pwelch(xw,50,25,1024,fs);
% figure(4);
% subplot(122);
% plot(f,pxx);
% xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot real pat of signal
[R,C]=size(Data2);
t=0:1/fs:(R-1)/fs;
% figure(5);
% subplot(121);
% for j=1:1
% plot(t,Data2(:,j));
% hold on;
% end;
% title('Real part of the signal');
% xlabel('time in sec');
% ylabel('amplitude');
%  
% i=1:999
% subplot(122);
% for j=1:1
% plot(i,Data2(:,j));
% hold on;
% end;
% title('Real part of the signal');
% xlabel('index');
% ylabel('amplitude');
%  
%%%%%%%%%%%%%%%%%%
NFFT=512;
S1fft=abs(fft(S(:,20),NFFT));
F=linspace(0,fs/2,(NFFT/2)+1);
figure(43);
subplot(121);
plot(F,S1fft(1:length(F)));title('Spectrum of filtered Signal');xlabel('Frequency in Hz');ylabel('FFT coeficient');
%%%%%%%%%%%%%%%%%%%%
NFFT=512;
S1fft=abs(fft(Data2(:,20),NFFT));
F=linspace(0,fs/2,(NFFT/2)+1);
figure(42);
subplot(121);
plot(F,S1fft(1:length(F)));title('Spectrum of raw Signal');xlabel('Frequency in Hz');ylabel('FFT coeficient');
 
% figure(6);
% subplot(122);
% plot(F,S1fft(1:length(F)));title('Spectrum of water sediment interface');xlabel('Frequency in Hz');ylabel('FFT coeficient');
%%%%%%%%%%%%%%%%%%%%%5
[pxx,f] = pwelch(Data2(:,666),100,50,1024,fs);
 figure(50);
 subplot(121);
 plot(f,20*log10(pxx/10^-6));
 title('PSD of transmitted Signal');
 xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 [pxx,f] = pwelch(S(:,666),100,50,1024,fs);
 
 subplot(122);
 plot(f,20*log10(pxx/10^-6));
 title('PSD of water sediment interface Signal');
 xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
%%%%%%%%%%%%%%
% [b a]=butter(2,[0.08 0.4],'bandpass');
%  S=filter(b,a,Data2(:,1:2162)); 
% [R,C]=size(S);
%%%%%%%%%%%%%%%%%%%%%%%%
[pxx,f] = pwelch(Data2(335:567,666),100,50,1024,fs);
figure(8);
subplot(121);
plot(f,20*log10(pxx/1));
title('PSD of original Signal after bandpass filter');
xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% NFFT=1024;
% S1fft=abs(fft(Data2(:,1),NFFT));
% F=linspace(0,fs/2,(NFFT/2)+1);
% subplot(122);
% plot(F,S1fft(1:length(F)));title('Spectrum of original Signal');xlabel('Frequency in Hz');ylabel('FFT coeficient');
% %%%%%%%%%%%%%%%%%%%%%%%%
%  [pxx,f] = pwelch(S(410:460,1),50,25,1024,fs);
%  figure(9);
%  subplot(121);
%  plot(f,pxx);
%  title('PSD of filtered water sediment interface Signal');
%  xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% % %%%%%%%%%%%%%%%%%%%%%%%
% NFFT=1024;
% S1fft=abs(fft(S,NFFT));
% F=linspace(0,fs/2,(NFFT/2)+1);
% figure(11);
% subplot(122);
% plot(F,S1fft(1:length(F)));title('Spectrum of filtered Signal');xlabel('Frequency in Hz');ylabel('FFT coeficient');
% %%%%%%%%%%%%%%%%%%%%%%%%%
% [pxx,f] = pwelch(S(480:525,1),45,22,1024,fs);
% figure(12);
% subplot(121);
% plot(f,20*log10(pxx/10^-6));
% title('PSD of original Signal');
% xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% % % hx = hilbert(S);                                     
% % % inst_amp = abs(hx);
% % % f1 = diff(unwrap(angle(hx)))/((1/fs)*2*pi);
% % % figure(2);
% % % plot(t(1:length(f1)),f1)
% % % xlabel('Time in sec');
% % % ylabel('frequency in Hz');title('Instantanious frequency');
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % [pxx,f] = pwelch(S,50,20,1024,fs);
% % % figure(3);
% % % subplot(122);
% % % plot(f,10*log10(pxx/10^-6));
% % % xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % corref=corr(Data2);
% % corref1=corrcoef(Data2);
% % imagesc(corref); % plot the matrix
% % % set(gca, 'XTick', 1:n); % center x-axis ticks on bins
% % % set(gca, 'YTick', 1:n); % center y-axis ticks on bins
% % % set(gca, 'XTickLabel', L); % set x-axis labels
% % % set(gca, 'YTickLabel', L); % set y-axis labels
% % title('Correlation Coefficient', 'FontSize', 14); % set title
% %colormapmap); % set the colorscheme
%  % enable colorbar
% % %Ensemble size selection using varience value
% % % ES=input('Entre value of ensemble size');
% ES=40;
% b = mod(C,ES)  consider for overall area
% Num=C-b;
% Data_new=S(:,1:Num);
% Dr=Num/ES;
% coldist = (size(Data_new,2)./Dr)*ones(1,Dr);
% D=mat2cell(Data_new,R,coldist);
% allined = cellfun(@(x) alline(x,ES),D,'UniformOutput',false);
 
%%%%%%%%%Plot allined signals
% plot(allined(1:2E-3*fs,1));
% hold on 
% plot(allined(1:2E-3*fs,2));
% hold on
%consider for overall area
% k=1;
% for i=1:54
%     a=allined{1,i};
%     a=mean(a');
%     b(k,1:999)=a;
%     k=k+1;
% end
%upto this
% for i=1:5
%     for j=1937:1:1941      
% final(:,i)=S(:,j);
%     
%     end;
% end;
%  
% Snew=S(:,542); % calculation of average of the 10 pings
%Savg= mean(Snew');
%Savg=Savg';
% Savg=Snew;
% %figure(13);
% %plot(1:999,Savg);
% %%%%%%%%%%%%%%%%%%%%%%%%
% t=0:1/fs:(R-1)/fs;
% hx = hilbert(Savg);     % hilbert of 1st ping                                
% inst_amp = abs(hx);
%  
% f_avginstphase= diff(unwrap(angle(hx)))/((1/fs)*2*pi);
% figure(333);
% % plot(d,f1);
% plot((t(1:length(f_instphase))*750),f_avginstphase); %hilbert transform
% xlabel('Depth in m');
% ylabel('frequency in Hz');
% title('Instantaneous frequency');
%  
%  %title('Instantaneous frequency of 20 averaged pings');
% tm=(t(1:length(f_avginstphase))*750);
% tm=tm';
% %%%%%%%%%%%%%%%%%
%figure(666)
% plot(tm(306:468,1),f_avginstphase(306:468,1));
 %%%%%%%%%%%%%%%%%
% x=tm(306:468,1);
% y=f_avginstphase(306:468,1);
% N=length(x);
% X=[ones(N,1),x];
% Y=y;
% phi=inv(X'*X)*X'*Y;
% %figure(5000);
% %plot(x,y,'bs',[0.5 5],phi(1)+phi(2)*[0.5 5],'-r');
% phi1=lsqcurvefit(@(x,xdata) myLinExample(x,xdata),[1;1],x,y);
 
%figure(5001);
%plot([0 20],phi1(1)+phi1(2)*[0.5 5],':g');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 
% [b,bint,r,rint,stats]= regress(f_avginstphase(306:468,1),tm(306:468,1));
% figure(21);
% plot(tm(306:468,1),f_avginstphase(306:468,1));
%%%%%%%%%%%%%%%% frequency analysis of average signal
% NFFT=1024;
% S1fft=abs(fft(Savg(437:450,1),NFFT)); 
% F=linspace(0,fs/2,(NFFT/2)+1);
% figure(14);
% subplot(121);
% plot(F,S1fft(1:length(F)));
% title(' Spectrum of averaged water sediment interface');xlabel('Frequency in Hz');ylabel('FFT coeficient');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  [pxxbu,fbu] = pwelch(Savg(437:450,1),10,5,1024,fs);% spetrun of top layer of sediment
%  [pxxbl,fbl] = pwelch(Savg(480:500,1),20,12,1024,fs); %spectrum of bottom layer of sediment
%   figure(15);
%  bldb=20*log10(pxxbl/10^-6)
%  budb=20*log10(pxxbu/10^-6);
%  plot(f,20*log10(pxxbu/10^-6));
%  hold on;
%  plot(f,20*log10(pxxbl/10^-6));
%  title('PSD of averaged upper water sediment interface and bottom sediment Signal');
%  xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
%  
 %%%%%%%%%%%%%%%%
% diff=bldu-bldb;
% A = cellfun(@(allined) mean(allined,2), allined,'UniformOutput',false);
% A_mat=cell2mat(A);
% [r,c]=size(A_mat);
% % 
% firstgr=Savg(17e-3*fs:20e-3*fs);
% figure(103);
% spectrogram(firstgr);
 
% % %%%%%%%%%Plot averaged signals
% figure
% for i=1:c
% plot(abs(A_mat(:,i)));
% hold on
% end
% % 
% % %%%%%%%%%Reflection coefficients
%Etx=1*10^9;
% Etx=31697863849222269704.04
% %%%%% francoise garrison eqn for calculation of absorption coefficient
% f=6;                                      %%%Centre frequency In khz 
% temp=25;                                  %%%%Average temperature
% W_depth=15;                               %%%%% SBES Data
% [ abs_coeff ] = fga(W_depth,temp,f);
% alpha=abs_coeff
% 
%      for  i=1:1                                   %%%C IS NUMBER OF COLUMNS i.e averaged pings
% rs1(:,i)=firstgr;
% Ers1(i)=sum(rs1(i).^2);
% Ers1_DB=db(Ers1,'power');
% R(i)=(2*W_depth*sqrt(Ers1(i)))./sqrt(Etx).*exp(-2*alpha*W_depth);
% R_Db(i)=20*log(R(i))
% end
% % %m=mean(R_Db)
% % 
% % x=1:ES:(ES*c);
% % t=0:0.5/92:0.5
% % figure(1)
% % g=plot(t,R_Db);
% % set(g,'linewidth',2);
% % xlabel('Time')
% % ylabel('Reflection Coefficients(Db)')
% %%%%%%%%%%%%%%%%%
% seg1=Savg(437:450,1); 
% seg2=Savg(480:500,1);
% [pxx,f] = pwelch(seg1,80,40,1024,fs);
% figure(5);
% subplot(121);
% plot(f,10*log10(pxx/10^-6));
% xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% figure(16);                      % spectrogram of top layer of sediment
%  spectrogram(seg1,13,8,1024,fs);
% figure(17);                      %spectrogram of bottom layer of sediment 
% spectrogram(seg2,20,10,1024,fs);
%  
 
% % %m=mean(R_Db)
% % 
% % x=1:ES:(ES*c);
% % t=0:0.5/92:0.5
% % figure(1)
% % g=plot(t,R_Db);
% % set(g,'linewidth',2);
% % xlabel('Time')
% % ylabel('Reflection Coefficients(Db)')
% %%%%%%%%%%%%%%%%%
% seg1=Savg(437:450,1); 
% seg2=Savg(480:500,1);
% % [pxx,f] = pwelch(seg1,80,40,1024,fs);
% % figure(5);
% % subplot(121);
% % plot(f,10*log10(pxx/10^-6));
% % xlabel('Frequency [Hz]');ylabel('PSD [dB re. 1 µPa ^2 /Hz]');
% figure(16);                      % spectrogram of top layer of sediment
%  spectrogram(seg1,13,8,1024,fs);
% figure(17);                      %spectrogram of bottom layer of sediment 
% spectrogram(seg2,20,10,1024,fs);
 
 
 
 

