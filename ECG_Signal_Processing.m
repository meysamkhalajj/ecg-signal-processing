clear
clc

%%% Plotting the ecg signal %%%
ecg = (load('Exercise_ecg.mat').ecg)';
Fs = 500;
Ts = 1/Fs;
L = numel(ecg);
t  = 0:Ts:L*Ts - Ts;
plot(t, ecg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FFT calculation of the ecg signal %%%
ECG = fft(ecg);
abs_ECG = abs(ECG);
L = numel(ECG);

Halved_ECG = ECG(1:L/2);
abs_Halved_ECG = abs(Halved_ECG);
Halved_f = 0:Fs/L:(Fs/2) - Fs/L;
f = 0:Fs/L:Fs-Fs/L;

subplot(2,1,1)
plot(Halved_f, abs_Halved_ECG)  % One-sided FFT

subplot(2,1,2)
plot(f, abs_ECG)                % Two_sided FFT                                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setting frequencies below 0.5 Hz to zero %%%
Filter = ones(1, L);
Filter(f<=0.5) = 0;
Filter(f>=(Fs-0.5)) = 0;
im_filter = ifft(Filter);  % Filter has a real impulse response

Filtered_ECG = ECG .* Filter;
filtered_ecg1 = ifft(Filtered_ECG);

subplot(2,1,1)
plot(t, ecg)

subplot(2,1,2)
plot(t, filtered_ecg1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Removing 50 Hz interference %%%
Fc = 50;
f1 = (Fc-1)/(Fs/2);
f2 = (Fc+1)/(Fs/2);
[b, a] = butter(4,[f1 f2],'stop');
[h, w]=freqz(b,a, [], Fs);
[gd, wd]=grpdelay(b, a, [], 500);

subplot(2,1,1)
plot(w, abs(h))

subplot(2,1,2)
plot(wd, gd)

filtered_ecg2 = filter(b, a, filtered_ecg1);

subplot(2,1,1)
plot(t, filtered_ecg1)

subplot(2,1,2)
plot(t, filtered_ecg2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Increasing the signal-to-noise ratio %%%
Fc = 35;  % 35 25 45 Comparision
f3 = Fc/(Fs/2);
[b, a] = butter(7,f3,"low");
[h, w] = freqz(b,a,[],500);
plot(w, abs(h))
filtered_ecg3 = filter(b,a,filtered_ecg2);

subplot(2,1,1)
% plot(t(t>=5 & t<=6), filtered_ecg2(t>=5 & t<=6))
plot(t, filtered_ecg2)
subplot(1,1,1)
% plot(t(t>=5 & t<=6), filtered_ecg3(t>=5 & t<=6))
plot(t, filtered_ecg3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the heart rate using autocorrelation %%%
[corr,lags] = xcorr(filtered_ecg3,filtered_ecg3,'normalized');
Corr_length = numel(corr);
plot(lags, corr)
xlabel('Lag (samples)')
ylabel('ACF')
ylim([-0.4 1.1])

Halved_corr = corr(Corr_length/2:end);

Max_HR = 90;
Min_HR = 60;
Search_min = floor(60*Fs/(Max_HR));
Search_max = floor(60*Fs/(Min_HR));

Ns = find(Halved_corr == max(Halved_corr(Search_min:Search_max))) - 1;
HR = 60*Fs/Ns;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fs = 500;
% fc = 0.5;
% 
% [b,a] = butter(12,fc/(Fs/2), "high");
% 
% [h,w] = freqz(b,a,'whole');
% 
% plot(w/pi*Fs/2, abs(h))
