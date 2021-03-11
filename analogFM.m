%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 1 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
[y,fs]=audioread('eric.wav');  %read audio file
figure(1)
subplot(2,2,1);
plot(y); %plot the signal of the audio file
title('Original signal in time Domain')
yf=fftshift(fft(y)); %use fourier transform to convert the signal from time domain to frequency domain
N = length(y);  %measure the length of the original signal
l = linspace(-fs/2,fs/2,N); %create a linspace by the size of the frequency of the original signal
subplot(2,2,2);
plot(abs(yf))  %plot the original signal in frequency domain
title('Original signal in frequency domain')

fcut = 4000;  %givin the cutting frequency by 4KHz
LPFsize=fix((N/fs)*(fcut*2) +1);   %the size of the low pass filter equal to the length of the original signal divided the frequency of the original signal times twice cutting freq 
pass=ones(1,LPFsize);
stop=zeros(1,(N-LPFsize)/2); 
fil=[stop pass stop];
filteredSignal = fil'.*yf; %multiply the complex conj transpose of the filter with the spectrum of the signal
subplot(2,2,4);
plot(l,abs(filteredSignal));
title('Filtered signal in frequency domain')

z=real(ifft(ifftshift(filteredSignal)));
subplot(2,2,3);
plot(z);
title('Filtered signal in time domain')

% sound(z,fs);

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 2 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 

fc=100000;  
fsnew=5*fc;  %new fs
x=resample(z',fsnew,fs);  %resample the filtered signal with the new fs
% lnew = linspace((-fsnew/2)-fc,(fsnew/2)+fc,length(x));
kf=60; %assumption for kf, B = kf/fsnew <<<1
t = linspace(0, length(x)/fsnew, length(x));
phase_dev = kf*cumsum(x);              % Integrating signal
NBFM=cos(2*pi*fc*t)-phase_dev.*sin(2*pi*fc*t);  
figure(2)
subplot(2,1,1);
plot(NBFM);
title('NBFM in time domain')
lnew = linspace((-fsnew/2)-fc,(fsnew/2)+fc,length(NBFM));
subplot(2,1,2);
plot(lnew,abs(fftshift(fft(NBFM))));
title('NBFM in frequency domain')

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% B<<<1 -> (delta f will be small, BW=2*fm, phase deviation will be small)

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 4 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% envelope=abs(hilbert(diff(NBFM)));
% figure(3);
% plot(envelope);
% title('Demodulation of the NBFM signal using a differentiator and an ED');

t = linspace(0, length(x)/fsnew, length(x));
c=(cumsum(x))/fsnew;

NBFM=cos((2*pi*fc*t)+(73.1*2*pi*c));  

dem = fsnew*diff(NBFM);

% dem=[NBFM(1) dem ];

envelope=abs(hilbert(dem));
envelope=0.01*(envelope - abs(mean(envelope)));
down=resample(envelope,fs,fsnew);
% sound(down,fs);
figure(4);
fvec = linspace(-fs/2,fs/2,length(down));  
plot(fvec,down);
ylim([-1 1])
title('NBFM demodulated signal');

% figure
% spectrum=fftshift(fft(NBFM)/fsnew);
% freq=linspace(-fsnew/2,fsnew/2,length(spectrum));
% plot(freq,abs(spectrum));
% title('spectrum of NBFM')


% 
% figure(2)
% subplot(2,1,1);
% plot(NBFM);
% title('NBFM in time domain')
% lnew = linspace((-fsnew/2)-fc,(fsnew/2)+fc,length(NBFM));
% subplot(2,1,2);
% plot(lnew,abs(fftshift(fft(NBFM))));
% title('NBFM in frequency domain')
% 
% %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 3 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% % B<<<1 -> (delta f will be small, BW=2*fm, phase deviation will be small)
% 
% 
% %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Part 4 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
% envelope=abs(hilbert(diff(NBFM)));
% figure(3);
% plot(envelope);
% title('Demodulation of the NBFM signal using a differentiator and an ED');
