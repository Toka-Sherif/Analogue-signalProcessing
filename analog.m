%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 1 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
[y,fs]=audioread('eric.wav');
ts=1/fs;
figure(1)
subplot(2,2,1);
plot(y);
title('Original signal in time Domain')
yf=fftshift(fft(y));
N_original = length(y);
l = linspace(-fs/2,fs/2,N_original);
subplot(2,2,2);
plot(l,abs(yf))
title('Original signal in frequency domain')

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 2 and 3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
fcut = 4000;
LPFsize=fix((N_original/fs)*(fcut*2) +1);
pass=ones(1,LPFsize);
stop=zeros(1,(N_original-LPFsize)/2); 
fil=[stop pass stop];
filtered_signal = fil'.*yf; %multiply the complex conj transpose of the filter with the spectrum of the signal
subplot(2,2,4);
plot(l,abs(filtered_signal));
title('Filtered signal in frequency domain')

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 4 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
filtered_signal_time=real(ifft(ifftshift(filtered_signal)));
subplot(2,2,3);
plot(l,filtered_signal_time);
title('Filtered signal in time domain')
%sound(filtered_signal_time,fs);

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 5 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
%XXXXXXXXXX resampling XXXXXXXXXXXXXX 
Fc=100000;
fc_new = 5*Fc ;
filtered_signal_resample_time = resample(filtered_signal_time,fc_new,fs);
f=linspace(-fc_new,fc_new,length(filtered_signal_resample_time));
figure;
subplot(3,1,1)
filted_signal_resample_freq = abs(fftshift(fft(filtered_signal_resample_time)));
plot(f,filted_signal_resample_freq);
title ('resampled filtered signal in frequency domain')
%XXXXXXXXXX implmenting DSB-SC modulation to resampled signal XXXXXXXXXXXXXX 
t_axis = length(filtered_signal_resample_time)/fc_new ;
t=linspace(0,t_axis,length(filtered_signal_resample_time));
F_carrirer=cos(2*pi*Fc*t);

DSB_SC_time = filtered_signal_resample_time'.*F_carrirer;

subplot(3,1,3)
plot(DSB_SC_time);
title ('Modulated  resampled filtered signal (DSB-SC)in time domain')

subplot(3,1,2)
DSB_SC_freq = abs(fftshift(fft(DSB_SC_time)));
fvec=linspace(-fs/2,fs/2,length(DSB_SC_freq));
plot(fvec,DSB_SC_freq );
title ('Modulated  resampled filtered signal (DSB-SC)in frequency domain')


%XXXXXXXXXX implmenting DSB-tC modulation to resampled  XXXXXXXXXXXXXX 

DSB_TC_time =((2*max(filtered_signal_resample))*(1+0.5*filtered_signal_resample))'.*F_carrirer(1:length(filtered_signal_resample)); 
figure;
subplot(3,1,1)
plot(f,filted_signal_resample_mag);
title ('resampled filtered signal in frequency domain')

subplot(3,1,3)
plot(DSB_TC_time);
title ('Modulated  resampled filtered signal (DSB-TC)in time domain')
% f=linspace((-fc_new/3.5)-Fc,fc_new/3.5+Fc,length(DSB_TC_time));

subplot(3,1,2)
DSB_TC_freq = abs(fftshift(fft(DSB_TC_time)));
f = linspace(-fs/2,fs/2,length(DSB_TC_freq));
plot(f,DSB_TC_freq);
title ('Modulated  resampled filtered signal (DSB-TC)in frequency domain')


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 6 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
envelope_DSB_SC = abs(hilbert(DSB_SC_time));
envelope_DSB_TC = abs(hilbert(DSB_TC_time));


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 7 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
figure;
subplot(2,1,1)
plot(envelope_DSB_SC);
title ('Modulated  resampled filtered signal (DSB-SC) envelope in time domain')
envelope_DSB_SC_resampled = resample(envelope_DSB_SC,fs,fc_new);
% sound(envelope_DSB_SC_resampled,fs); 
% audiowrite('envelope_DSB_SC_resampled.wav',envelope_DSB_SC_resampled,fs);

%-------comment on sound----------
%the sound is unclear (interfed sound)
%envelope detector can't be used with DSB-SC modulated signal
%------------------------

subplot(2,1,2)
plot(envelope_DSB_TC);
title ('Modulated  resampled filtered signal (DSB-TC) envelope in time domain')
envelope_DSB_TC_resampled=resample(envelope_DSB_TC,fs,fc_new);
% sound(envelope_DSB_TC_resampled,fs);
% audiowrite('envelope_DSB_TC_resampled.wav',envelope_DSB_TC_resampled,fs);

%-------comment on sound----------
%the sound is low added to it a long wistle song 
%envelope detector can be used with DSB-SC modulated signal
%------------------------



%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 8 applied on DSB_TC XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
%add noise to DSB_TC signal
snr_0=awgn( DSB_TC_time ,0); %for SNR=0 %
snr_10=awgn(DSB_TC_time ,10);  %for SNR=10 %
snr_30=awgn(DSB_TC_time ,30);  %for SNR=30 %

%recieve them by envelope detctor and play back each 
envelope_snr_0=abs(hilbert(snr_0));
envelope_snr_0_resample=resample(envelope_snr_0,fs,fc_new);
% sound(envelope_snr_0_resample,fs); 
audiowrite('envelope_snr_0_resample.wav',envelope_snr_0_resample,fs);

envelope_snr_10=abs(hilbert(snr_10));
envelope_snr_10_resample=resample(envelope_snr_10,fs,fc_new);
% sound(envelope_snr_10_resample,fs); 
audiowrite('envelope_snr_10_resample.wav',envelope_snr_10_resample,fs);

envelope_snr_30=abs(hilbert(snr_30));
envelope_snr_30_resample=resample(envelope_snr_30,fs,fc_new);
% sound(envelope_snr_30_resample,fs); 
audiowrite('envelope_snr_30_resample.wav',envelope_snr_30_resample,fs);

%sketch it each in time domain
figure;
subplot(3,1,1)
plot(envelope_snr_0_resample);
title ('SNR = 0')
subplot(3,1,2)
plot(envelope_snr_10_resample);
title ('SNR = 10')
subplot(3,1,3)
plot(envelope_snr_30_resample);
title ('SNR = 30')


%----------conclusion----------
% As signal to noise ratio increase 
% the output signal is returning to the original signal(where was no interfernce)  
%------------------------------

% at 100 db the ouput signal is nearly equal to the orignal signal
% snr_100 = awgn(DSB_TC_time ,100);
% envelope_snr_100=abs(hilbert(snr_100));
% envelope_snr_100_resample = resample(envelope_snr_100,fs,fc_new);
% sound(envelope_snr_100_resample,fs); 




%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 9 applied on DSB_SC XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
%XXXXXXXXXXXXXXXXXXXXXXXXXcoherent detection with no noises XXXXXXXXXXXXXXXX
figure;
subplot(2,1,1)
f=linspace(-fs/2,fs/2,length(filtered_signal_resample));
plot(f,filted_signal_resample_mag);
title('the original signal in frquency domain');


DSBC_SC_coh_time = DSB_SC_time.*F_carrirer;
DSBC_SC_coh_freq = fftshift(fft(DSBC_SC_coh_time));
N = length(DSBC_SC_coh_freq);
LPFsize=fix((N/fs)*(fcut*2));
pass=ones(1,LPFsize);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
length(filter)
% applying LPF
DSB_SC_filtered_freq=filter.*DSBC_SC_coh_freq;
DSBC_SC_coh_filtered_freq_resampled = resample(DSB_SC_filtered_freq,fs,fc_new);
subplot(2,1,2)
y=length(DSBC_SC_coh_filtered_freq_resampled);
f_vec = linspace(-fs/2,fs/2,y);
plot(f_vec , abs(DSBC_SC_coh_filtered_freq_resampled));
title('Message after using coherent detector for (DSBSC) freq-domain')
% subplot(2,1,2)
% DSB_SC_filtered_time = real(ifft(ifftshift(DSB_SC_filtered_freq)));
% DSB_SC_filtered_time_resampled = resample(DSB_SC_filtered_time,fs,fc_new);
% plot(DSB_SC_filtered_time_resampled);
% title('Message after using coherent detector for (DSBSC) time-domain')


%XXXXXXXXXXXXXXXXXXXXXcoherent detection with noises XXXXXXXXXXXXXXXX
snr_0=awgn( DSB_SC_time ,0); %for SNR=0 %
snr_10=awgn(DSB_SC_time ,10);  %for SNR=10 %
snr_30=awgn(DSB_SC_time ,30);  %for SNR=30 %
%for clarifying the observation that as power increase the resulted signal
%as similar to original signal
snr_100=awgn(DSB_SC_time ,100);  %for SNR=100 %


S_snr_0=snr_0.*F_carrirer;
S_snr_0_freq = fftshift(fft(S_snr_0)) ; 
N = length(S_snr_0_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
S_snr_0_filtered_freq=filter.*S_snr_0_freq;
S_snr_0_filtered_resample_freq =resample(S_snr_0_filtered_freq,fs,fc_new);
%sound(S_snr_0_filtered_resample,fs);
audiowrite('S_snr_0_filtered_resample.wav',S_snr_0_filtered_resample_freq,fs);
figure;
subplot(4,2,1)
y=length(S_snr_0_filtered_resample_freq);
fvec =linspace(-fs/2,fs/2,y) ;  
plot(fvec,abs(S_snr_0_filtered_resample_freq));
title ('SNR = 0 (DSBSC)coherent detection freq domain')
subplot(4,2,2)
S_snr_0_filtered_time =real(ifft(ifftshift(S_snr_0_filtered_freq)));
S_snr_0_filtered_resambled_time = resample(S_snr_0_filtered_time,fs,fc_new);
plot(S_snr_0_filtered_resambled_time);
title('SNR = 0 (DSBSC)coherent detection time domain')

S_snr_10=snr_10.*F_carrirer;
S_snr_10_freq = fftshift(fft(S_snr_10)) ; 
N = length(S_snr_10_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
S_snr_10_filtered_freq=filter.*S_snr_10_freq;
S_snr_10_filtered_resample_freq =resample(S_snr_10_filtered_freq,fs,fc_new);
%sound(S_snr_10_filtered_resample,fs);
audiowrite('S_snr_10_filtered_resample.wav',S_snr_10_filtered_resample_freq,fs);
subplot(4,2,3)
y=length(S_snr_10_filtered_resample_freq);
fvec =linspace(-fs/2,fs/2,y) ;  
plot(fvec,abs(S_snr_10_filtered_resample_freq));
title ('SNR = 10 (DSBSC)coherent detection freq domain')
subplot(4,2,4)
S_snr_10_filtered_time =real(ifft(ifftshift(S_snr_10_filtered_freq)));
S_snr_10_filtered_resambled_time = resample(S_snr_10_filtered_time,fs,fc_new);
plot(S_snr_10_filtered_resambled_time);
title('SNR = 10 (DSBSC)coherent detection time domain')


S_snr_30=snr_30.*F_carrirer;
S_snr_30_freq = fftshift(fft(S_snr_30)) ; 
N = length(S_snr_30_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
S_snr_30_filtered_freq=filter.*S_snr_30_freq;
S_snr_30_filtered_resample_freq =resample(S_snr_30_filtered_freq,fs,fc_new);
sound(S_snr_30_filtered_resample,fs);
audiowrite('S_snr_30_filtered_resample.wav',S_snr_30_filtered_resample_freq,fs);
subplot(4,2,5)
y=length(S_snr_30_filtered_resample_freq);
fvec =linspace(-fs/2,fs/2,y) ;  
plot(fvec,abs(S_snr_30_filtered_resample_freq));
title ('SNR = 30 (DSBSC)coherent detection freq domain')
subplot(4,2,6)
S_snr_30_filtered_time =real(ifft(ifftshift(S_snr_30_filtered_freq)));
S_snr_30_filtered_resambled_time = resample(S_snr_30_filtered_time,fs,fc_new);
plot(S_snr_30_filtered_resambled_time);
title('SNR = 30 (DSBSC)coherent detection time domain')


%--------------------------------------------------------------
%for clarifying the observation that as power increase the resulted signal
%as similar to original signal
S_snr_100=snr_100.*F_carrirer;
S_snr_100_freq = fftshift(fft(S_snr_100)) ; 
N = length(S_snr_100_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
S_snr_100_filtered_freq=filter.*S_snr_100_freq;
S_snr_100_filtered_resample_freq =resample(S_snr_100_filtered_freq,fs,fc_new);
%sound(S_snr_100_filtered_resample,fs);
audiowrite('S_snr_100_filtered_resample.wav',S_snr_100_filtered_resample_freq,fs);
subplot(4,2,7)
y=length(S_snr_100_filtered_resample_freq);
fvec =linspace(-fs/2,fs/2,y) ;  
plot(fvec,abs(S_snr_100_filtered_resample_freq),'r');
title ('SNR = 100 (DSBSC)coherent detection freq domain')
subplot(4,2,8)
S_snr_100_filtered_time =real(ifft(ifftshift(S_snr_100_filtered_freq)));
S_snr_100_filtered_resambled_time = resample(S_snr_100_filtered_time,fs,fc_new);
plot(S_snr_100_filtered_resambled_time,'r');
title('SNR = 100 (DSBSC)coherent detection time domain')
%--------------------------------------------------------------



%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 10 applied on DSB_SC XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% coherent detection with frequency error
Carrier_with_freqerror=cos((Fc+100)*2*pi*t);
DSBSC_with_freqerror_time= DSB_SC_time.*Carrier_with_freqerror;
DSBSC_with_freqerror_freq = fftshift(fft(DSBSC_with_freqerror_time));
N = length(DSBSC_with_freqerror_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
DSBSC_with_freqerror_filtered=filter.*DSBSC_with_freqerror_freq;
DSBSC_with_freqerror_filtered_time = real(ifft(ifftshift(DSBSC_with_freqerror_filtered)));
DSBSC_with_freqerror_filtered_resampled_time =resample(DSBSC_with_freqerror_filtered_time,fs,fc_new);
% sound(DSBSC_with_freqerror_filtered_resampled,fs);

figure; 
subplot(2,1,1)
plot(DSBSC_with_freqerror_filtered_resampled_time);
title('Message after using coherent detector with frequency error in Time domain')

subplot(2,1,2)
DSBSC_with_freqerror_filtered_resambled_freq = resample(DSBSC_with_freqerror_filtered , fs ,fc_new);
y = length(DSBSC_with_freqerror_filtered_resambled_freq);
fvec = linspace(-fs/2,fs/2,y);
plot(fvec ,abs(DSBSC_with_freqerror_filtered_resambled_freq));
title('Message after using coherent detector with frequency error in Frequency domain')


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 11 applied on DSB_SC XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% coherent detection with phase error
Carrier_with_phase_error=cos((20*pi/180)+(Fc*2*pi*t));
DSBSC_with_PhaseError= DSB_SC_time.*Carrier_with_phase_error;
DSBSC_with_PhaseError_freq = fftshift(fft(DSBSC_with_PhaseError));
N = length(DSBSC_with_PhaseError_freq);
LPFsize=fix((N/fs)*(fcut*2) +1);
pass=ones(1,LPFsize+1);
stop=zeros(1,(N-LPFsize)/2); 
filter=[stop pass stop];
DSBSC_with_PhaseError_filtered=DSBSC_with_PhaseError_freq.*filter;
DSBSC_with_freqerror_filtered_time = real(ifft(ifftshift(DSBSC_with_PhaseError_filtered)));
%sound(DSBSC_with_PhaseError_filtered_resampled,fs);

figure; 
subplot(2,1,1)
DSBSC_with_freqerror_filtered_resampled_time =resample(DSBSC_with_freqerror_filtered_time,fs,fc_new);
plot(DSBSC_with_freqerror_filtered_resampled_time);
title('Message after using coherent detector with phase error in Time domain')

subplot(2,1,2)
DSBSC_with_freqerror_filtered_resampled_freq =resample(DSBSC_with_PhaseError_filtered,fs,fc_new);
y = length(DSBSC_with_freqerror_filtered_resampled_freq);
fvec = linspace(-fs/2,fs/2,y);
plot(fvec,abs(DSBSC_with_freqerror_filtered_resampled_freq));
title('Message after using coherent detector with phase error in Frequency domain')

